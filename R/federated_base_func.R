
## return a table of summary statistics for federated results 
## input args: coef<vector>: a vector of pooled coefficients 
##              cov<matrix>: a matrix of pooled covariance 
## 
## return: a summary table of pooled Estimate, Std. Error, t stats and p value
make_summary_table <- function(coef, cov) {
  sd <- diag(cov)**0.5
  z.value <- coef/sd
  p.value <- 2*(1-pnorm(abs(z.value)))
  summary.table <- cbind(coef, sd, z.value, p.value)
  colnames(summary.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|z|)")
  return(summary.table)
}

## reformat the summary table 
## input args: summary_table<matrix>: a summary table of pooled results; return by make_summary_table function 
## 
## return: kable stylized summary table 
Stylize <- function(summary_table){
  kable(summary_table) %>% kable_styling(full_width = F) #position = "left"
}


## compute estimated probabilites with estimated coefficients 
## input args: X<matrix>: a design matrix X
##          beta<vector>: an estimated coefficients 
##       subset<boolean>: indicator if only a subset of coefficients is used; default is FALSE 
## return: estimated probabilites, i.e. expit(X %*% beta) if subset = FALSE 
eval_prob <- function(X, beta, subset=FALSE) {
  if (subset) {
    covars <- colnames(X)
    return(1/(1+exp(-X %*% beta[covars])))
  } else {
    return(1/(1+exp(-X %*% beta)))
  }
}

## get the inverse propensity weights
##    input args: treat<vector>: a vector of treatment assignments  
##                    e<vector>: a vector of estimated propensity scores 
##               estimand<char>: "ATE" or "ATT"
##  propensity_threshold<float>: enforced lower bound for propensity scores  
## 
## return: weights that will be used in estimating IPW-MLE variance
propensity_weight <- function(treat, e, estimand="ATE", propensity_threshold = 0, normalize_propensity_weights) {
  e <- pmax(e, propensity_threshold)
  if (estimand == "ATT") {
    w <- treat + (1 - treat) * e / (1 - e) 
  } else {
    w <- treat / e + (1 - treat) / (1 - e)
  }
  return (w)
}

## get the alternative inverse propensity weights (only for ATT weights)
##    input args: treat<vector>: a vector of treatment assignments  
##                    e<vector>: a vector of estimated propensity scores 
##               estimand<char>: "ATE" or "ATT"
## 
## return: weights that will be used in estimating IPW-MLE variance (for ATT)
alternative_propensity_weight <- function(treat, e, estimand="ATE") {
  if (estimand == "ATT") {
    w <- (1 - treat) / (1 - e)
  } else {
    w <- treat / e + (1 - treat) / (1 - e)
  }
  return (w)
}

## evaluate the log likelihood function
##    input args: beta<vector>: a vector of estimated coefficents 
##                   y<vector>: a vector of binary outcome values 
##                   X<matrix>: a design matrix  
##                   w<vector>: a vector of weights 
## 
## return: log likelihood evaluated at the estimated coefficients 
log_lik <- function(beta, y, X, w) {
  prob <- eval_prob(X, beta)
  ll_vec <- y * log(prob) + (1 - y) * log(1 - prob)
  ll <- sum(w * ll_vec)
  ll
}

## evaluate the gradient of log likelihood function
##    input args: beta<vector>: a vector of estimated coefficents 
##                   y<vector>: a vector of binary outcome values 
##                   X<matrix>: a design matrix  
##                   w<vector>: a vector of weights 
## 
## return: gradient of log likelihood function evaluated at the estimated coefficients 
grad_log_lik <- function(beta, y, X, w) {
  prob <- eval_prob(X, beta)
  grad_ll <- diag(as.vector(w) * as.vector(y - prob)) %*% X 
  grad_ll
}

## evaluate the outer product of gradient of log likelihood function
##    input args: beta<vector>: a vector of estimated coefficents 
##                   y<vector>: a vector of binary outcome values 
##                   X<matrix>: a design matrix  
##                   w<vector>: a vector of weights 
##                wsq<boolean>: indicator if weighted log likelihood is used; default is FALSE
## 
## return: the outer product of gradient of log likelihood function evaluated at the estimated coefficients 
outer_grad_log_lik <- function(beta, y, X, w, wsq=FALSE) {
  prob <- eval_prob(X, beta)
  grad_ll <- diag(as.vector(y - prob)) %*% X
  
  # if weighted log likelihood is used
  if (wsq) {
    outer_grad_ll <- t(grad_ll) %*% diag(as.vector(w)**2) %*% grad_ll
  } else {
    outer_grad_ll <- t(grad_ll) %*% diag(as.vector(w)) %*% grad_ll
  }
  outer_grad_ll
}

## evaluate the hessian of log likelihood function
##    input args: beta<vector>: a vector of estimated coefficents 
##                   y<vector>: a vector of binary outcome values 
##                   X<matrix>: a design matrix  
##                   w<vector>: a vector of weights 
## 
## return: the hessian matrix of log likelihood function evaluated at the estimated coefficients 
hessian_log_lik <- function(beta, y, X, w) {
  prob <- eval_prob(X, beta)
  hessian_ll <- -t(X) %*%  diag(as.vector(w) * as.vector((1 - prob) * prob)) %*%  X 
  hessian_ll
}

## evaluate the gradient and hessian of log likelihood function
##    input args: est<vector>: a vector of estimated coefficents 
##             this.y<vector>: a vector of binary outcome values 
##             this.X<matrix>: a design matrix  
##             this.w<vector>: a vector of weights 
##               wsq<boolean>: indicator if weighted log likelihood is used; default is FALSE
##          return_V<boolean>: indicator if covariance matrix is returned; default is FALSE 
##            robust<boolean>: indicator if robust Huber (sandwich) variance estimator is used; default is TRUE
##           dataset_idx<int>: index of this dataset 
##          unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##    unstable_covars<vector>: if models are unstable, a vector of unstable covariate names 
##          num_datasets<int>: number of datasets
## 
## return: the gradient, the outer product of gradient, hessian of log likelihood and covariance matrix 
eval_grad_hessian <- function(est, this.y, this.X, this.w, wsq=FALSE, return_V=FALSE, robust=TRUE, dataset_idx=NULL, 
                              unstable=FALSE, unstable_covars=NULL, num_datasets=NULL) {

  if(unstable) {
    this.X <- AddInteractCols(this.X, dataset_idx, unstable_covars, num_datasets)
    this.X <- this.X[, names(est)]
  }
  
  grad <- grad_log_lik(est, y=this.y, X=this.X, w=this.w)    
  outer.grad <- outer_grad_log_lik(est, y=this.y, X=this.X, w=this.w, wsq=wsq)
  hessian <- hessian_log_lik(est, y=this.y, X=this.X, w=this.w)
  
  if (return_V) {
    if (robust) {
      V <- solve(hessian) %*% outer.grad %*% solve(hessian)
    } else {
      V <- solve(-hessian)
    }
    return(list(grad=grad, outer.grad=outer.grad, hessian=hessian, V=V))
  } else {
    return(list(grad=grad, outer.grad=outer.grad, hessian=hessian))
  }
}


## add interaction columns to the design matrix for the unstable covariates
##    input args: this.X<matrix>: a design matrix  
##              dataset_idx<int>: index of this dataset 
##       unstable_covars<vector>: if models are unstable, a vector of unstable covariate names; if stable, set to NULL
##          num_datasets<int>: number of datasets
## 
## return: a zero-padded design matrix 
AddInteractCols <- function(this.X, dataset_idx, unstable_covars, num_datasets) {
  this.stable_cols <- as.matrix(this.X[,!colnames(this.X) %in% c("intercept", unstable_covars)]) 
  colnames(this.stable_cols) <- colnames(this.X)[!colnames(this.X) %in% c("intercept", unstable_covars)]
  
  # if idx is not given, append unstable columns corresponding to all datasets
  if (is.null(dataset_idx)) {
    intercept.unstable_covars <- c("intercept", unstable_covars)
    this.unstable_datasets <- NULL
    for (c in 1:length(intercept.unstable_covars)) {
      this.unstable_col <- as.matrix(this.X[, intercept.unstable_covars[c]])
      rownames(this.unstable_col) <- NULL
      this.unstable_cols <- this.unstable_col[, rep(1, each=num_datasets)]
      colnames(this.unstable_cols) <- paste("I", 1:num_datasets, "_", intercept.unstable_covars[c], sep ="")
      this.unstable_datasets <- cbind(this.unstable_datasets, this.unstable_cols)
    }
    interact.X <- cbind(this.stable_cols, this.unstable_datasets)
  
  # if idx is given, append unstable columns corresponding to specific dataset 
  } else {
    this.unstable_cols <- as.matrix(this.X[, c("intercept", unstable_covars)])
    rownames(this.unstable_cols) <- NULL
    colnames(this.unstable_cols) <- paste("I", dataset_idx, "_", c("intercept", unstable_covars), sep ="")
    
    # fill with 0 for other interation columns corresponding to other datasets 
    filled_zero_names <- apply(expand.grid(paste("I", c(1:num_datasets)[c(1:num_datasets) != dataset_idx], sep = ""),
                                           c("intercept", unstable_covars)), 1, function(x) paste0(x, collapse="_"))
    filled_zero_cols <- matrix(0, ncol = length(filled_zero_names), nrow = nrow(this.X))
    colnames(filled_zero_cols) <- filled_zero_names
    interact.X <- cbind(this.stable_cols, this.unstable_cols, filled_zero_cols)
  }
  
  # reorder the columns: I1_covar1, I2_covar1, I3_covar1, I1_covar2, I2_covar2, I3_covar2.  
  unstable_coeff_names <- apply(expand.grid(paste("I", 1:num_datasets, sep = ""),
                                            c("intercept", unstable_covars)), 1, function(x) paste0(x, collapse="_"))
  coeff_names <- c(setdiff(colnames(this.X), c("intercept", unstable_covars)), unstable_coeff_names)
  interact.X <- interact.X[, coeff_names]
  return(interact.X)
}


## compute the IPW-MLE estimator variance if propensity scores are estimated
##    input args: e.grad<vector>: gradient of log likelihood function for treatment model 
##          e.outer.grad<matrix>: a matrix of outer product of gradient of log likelihood function for treatment model 
##             e.hessian<matrix>: hessian of log likelihood function for treatment model 
##                y.grad<vector>: weights * gradient of log likelihood function for outcome model (first two component in C1)
##            y.grad.alt<vector>: alternative weights * gradient of log likelihood function for outcome model (first two component in C2)
##          y.outer.grad<matrix>: a matrix of outer product of gradient of log likelihood function for outcome model 
##             y.hessian<matrix>: hessian of log likelihood function for outcome model 
##                    n_obs<int>: number of observations 
##             unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##
## return: matrices for computing IPW-MLE variance
calc_ht_cov <- function(e.grad, e.outer.grad, e.hessian, y.grad, y.grad.alt, y.outer.grad, y.hessian, n_obs, unstable=FALSE) {
  
  A0 <- y.hessian/n_obs
  D0 <- y.outer.grad/n_obs
  C1 <- t(y.grad.alt) %*% e.grad/n_obs
  C2 <- t(y.grad) %*% e.grad/n_obs
  
  if (!unstable) {
    M0.bread <- solve(e.hessian/n_obs)
    M0.meat <- e.outer.grad/n_obs
    M0.inv <- M0.bread %*% M0.meat %*% M0.bread
    
    B0 <- D0 - C1 %*% M0.inv %*% t(C2) - C2 %*% M0.inv %*% t(C1) + C2 %*% M0.inv %*% t(C2)
    A0.inv <- solve(A0)
    V <- A0.inv %*% B0 %*% A0.inv/n_obs 
    return(list(V=V, A0=A0, B0=B0, C1=C1, C2=C2, D0=D0, M0.inv=M0.inv, M0.bread=M0.bread, M0.meat=M0.meat, n_obs=n_obs))
  } else {
    return(list(A0=A0, C1=C1, C2=C2, D0=D0, n_obs=n_obs))
  }
}

## compute the HT estimator variance if propensity scores are known
##    input args: y.outer.grad<matrix>: a matrix of outer product of gradient of log likelihood function for outcome model 
##                   y.hessian<matrix>: hessian of log likelihood function for outcome model 
##                          n_obs<int>: number of observations 
##
## return: matrices for computing IPW-MLE variance
calc_ht_cov_simple <- function(y.outer.grad, y.hessian, n_obs) {
  
  A0 <- y.hessian/n_obs
  D0 <- y.outer.grad/n_obs
  B0 <- D0 
  
  A0.inv <- solve(A0)
  V <- A0.inv %*% B0 %*% A0.inv/n_obs
  return(list(V=V, A0=A0, B0=B0, D0=D0, n_obs=n_obs))
}

############# Base Estimation Functions (MLE, IPW-MLE, AIPW) #############

## obtain data set-specific MLE estimation
##            input args: this.y<vector>: a vector of binary outcome values 
##                        this.X<matrix>: a design matrix  
##                        this.w<vector>: a vector of weights 
##                          wsq<boolean>: indicator if weighted log likelihood is used; default is FALSE
##                       robust<boolean>: indicator if robust Huber (sandwich) variance estimator is used; default is TRUE
##                    treat.reg<boolean>: indicator if estimate treatment model; default is FALSE
##                      dataset_idx<int>: index of this dataset 
##                     unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##               unstable_covars<vector>: if models are unstable, a vector of unstable covariate names 
##                  reg.formula<formula>: a outcome regression formula 
##                full_reg.formula<list>: list of all outcome regression formulas
##         reg.formula.unstable<formula>: a outcome regression formula (if unstable models)
##       full_reg.formula.unstable<list>: list of all outcome regression formulas (if unstable models)
##            treat.reg.formula<formula>: a treatment regression formula 
##          full_treat.reg.formula<list>: list of all treatment regression formulas
##   treat.reg.formula.unstable<formula>: a treatment regression formula (if unstable models)
## full_treat.reg.formula.unstable<list>: list of all treatment regression formulas (if unstable models)
##
## return: the MLE coefficients, zero-padded coefficients, gradient, outer product of gradient, hessian, 
##         estimated probabilities, variance covariance matrix 
est_mle <- function(this.y, this.X, this.w=NULL, wsq=FALSE, robust=TRUE, treat.reg=FALSE, 
                    dataset_idx=NULL, num_datasets=NULL, unstable=FALSE, unstable_covars=NULL, 
                    reg.formula=NULL, full_reg.formula=NULL, 
                    reg.formula.unstable=NULL, full_reg.formula.unstable=NULL,
                    treat.reg.formula=NULL, full_treat.reg.formula=NULL, 
                    treat.reg.formula.unstable=NULL, full_treat.reg.formula.unstable=NULL) {
  
  if (is.null(this.w)) {
    this.w <- rep(1, length(this.y))
  }
  if(!unstable) {
    data <- data.frame(this.y, this.w, this.X)
  } else {
    X.interact <- AddInteractCols(this.X, dataset_idx, unstable_covars, num_datasets)
    data <- data.frame(this.y, this.w, X.interact)
  }
  if (!unstable) {
    if (treat.reg) {
      if (is.null(dataset_idx)) {
        this.reg.formula <- treat.reg.formula
      } else {
        this.reg.formula <- full_treat.reg.formula[[dataset_idx]]
      }
      colnames(data)[1] <- "treat"
    } else {
      if (is.null(dataset_idx)) {
        this.reg.formula <- reg.formula
      } else {
        this.reg.formula <- full_reg.formula[[dataset_idx]]
      }
      colnames(data)[1] <- "y"
    }
  } else {
    if (treat.reg) {
      if (is.null(dataset_idx)) {
        this.reg.formula <- treat.reg.formula.unstable
      } else {
        this.reg.formula <- full_treat.reg.formula.unstable[[dataset_idx]]
      }
      colnames(data)[1] <- "treat"
    } else {
      if (is.null(dataset_idx)) {
        this.reg.formula <- reg.formula.unstable
      } else {
        this.reg.formula <- full_reg.formula.unstable[[dataset_idx]]
      }
      colnames(data)[1] <- "y"
    }
  }
  dstrat <- svydesign(id=~1, weights=this.w, data=data, variables=this.reg.formula)
  this.fit <- svyglm(formula=this.reg.formula, design=dstrat, family=binomial)
  est <- coef(this.fit)
  est.prob <- eval_prob(as.matrix(data[, names(est)]), est)
  
  if (!unstable) {
    out_grad_hessian <- eval_grad_hessian(est, this.y, this.X, this.w, wsq=wsq, return_V=TRUE, robust=robust)
  } else {
    out_grad_hessian <- eval_grad_hessian(est, this.y, this.X, this.w, wsq=wsq, return_V=FALSE,
                                          dataset_idx=dataset_idx, robust=robust,
                                          unstable=TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets)
  }
  return(list(est=est, grad=out_grad_hessian$grad, outer.grad=out_grad_hessian$outer.grad,
              hessian=out_grad_hessian$hessian, est.prob=est.prob, V=out_grad_hessian$V))
}


## obtain data set-specific IPW-MLE estimation 
## (*) this function will be used in generating dataset-specific results and oracle results 
## (*) this function will not be used in generating federated results
## 
##            input args: this.y<vector>: a vector of binary outcome values 
##                  this.X.tilde<matrix>: a design matrix (include column of treatment)
##                        this.X<matrix>: a design matrix (baseline covariates only)
##                    this.treat<vector>: a vector of treatment assignments
##                        this.w<vector>: a vector of weights 
##                          wsq<boolean>: indicator if weighted log likelihood is used; default is FALSE
##                       robust<boolean>: indicator if robust Huber (sandwich) variance estimator is used; default is TRUE
##                    treat.reg<boolean>: indicator if estimate treatment model; default is FALSE
##                      dataset_idx<int>: index of this dataset 
##                     unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##               unstable_covars<vector>: if models are unstable, a vector of unstable covariate names 
##                  reg.formula<formula>: a outcome regression formula 
##            treat.reg.formula<formula>: a treatment regression formula 
##                        estimand<char>: "ATE" or "ATT"; default is "ATE"
##         estimated_propensity<boolean>: if propensity scores are estimated or known; default is TRUE
## 
## return: the IPW-MLE coefficients, variance and estimated probabilities
est_ht <- function(this.y, this.X.tilde, this.X, this.treat, reg.formula, treat.reg.formula, 
                   estimand="ATE", estimated_propensity=TRUE) {
  
  n_obs <- length(this.y)
  this.w <- rep(1, length(this.y))
  
  e_mle_out <- est_mle(this.treat, this.X, this.w, treat.reg=TRUE, treat.reg.formula=treat.reg.formula)
  e.grad <- e_mle_out$grad; e.outer.grad <- e_mle_out$outer.grad; e.hessian <- e_mle_out$hessian; est.e.prob <- e_mle_out$est.prob
  
  this.w.y <- propensity_weight(this.treat, est.e.prob, estimand = estimand)
  this.w.alt.y <- alternative_propensity_weight(this.treat, est.e.prob, estimand = estimand)
  
  y_mle_out <- est_mle(this.y, this.X.tilde, this.w.y, wsq = TRUE, reg.formula=reg.formula)
  est.y.coef <- y_mle_out$est
  y.grad <- y_mle_out$grad; y.outer.grad <- y_mle_out$outer.grad; y.hessian <- y_mle_out$hessian; est.y.prob <- y_mle_out$est.prob
  y.grad.alt <- diag(as.vector(this.w.alt.y) / as.vector(this.w.y)) %*% y_mle_out$grad
  
  if (estimated_propensity) {
    ht_cov_out <- calc_ht_cov(e.grad, e.outer.grad, e.hessian, y.grad, y.grad.alt, y.outer.grad, y.hessian, n_obs)
  } else {
    ht_cov_out <- calc_ht_cov_simple(y.grad, y.outer.grad, y.hessian, n_obs)
  }
  return(list(est=est.y.coef, V=ht_cov_out$V, est.prob=est.y.prob))
}


## obtain data set-specific AIPW estimation
##            input args: this.y<vector>: a vector of binary outcome values 
##                  this.X.tilde<matrix>: a design matrix (include column of treatment)
##                        this.X<matrix>: a design matrix (baseline covariates only)
##                    this.treat<vector>: a vector of treatment assignments
##                        estimand<char>: "ATE" or "ATT"; default is "ATE"
##                      dataset_idx<int>: index of this dataset 
##                 pooled.y.coef<vector>: federated outcome model coefficients 
##                 pooled.e.coef<vector>: federated treatment model coefficients
##                  reg.formula<formula>: a outcome regression formula 
##                full_reg.formula<list>: list of all outcome regression formulas
##            treat.reg.formula<formula>: a treatment regression formula 
##          full_treat.reg.formula<list>: list of all treatment regression formulas
##
## return: the AIPW estimates, variance and influence functions 
est_aipw <- function(this.y, this.X.tilde, this.X, this.treat, estimand="ATE", dataset_idx=NULL, 
                     pooled.y.coef=NULL, pooled.e.coef=NULL, 
                     reg.formula=NULL, full_reg.formula=NULL, 
                     treat.reg.formula=NULL, full_treat.reg.formula=NULL) {
  
  n_obs <- length(this.y); this.w <- rep(1, length(this.y))
  
  if (is.null(pooled.e.coef)) {
    e_mle_out <- est_mle(this.treat, this.X, this.w, dataset_idx=dataset_idx, treat.reg=TRUE, 
                         treat.reg.formula=treat.reg.formula, full_treat.reg.formula=full_treat.reg.formula)
    est.e.prob <- e_mle_out$est.prob
    y_mle_out <- est_mle(this.y, this.X.tilde, this.w, wsq=TRUE, dataset_idx=dataset_idx, 
                         reg.formula=reg.formula, full_reg.formula=full_reg.formula)
    est.y.coef <- y_mle_out$est    
  } else {
    this.pool_e_coef <- pooled.e.coef[rownames(pooled.e.coef) %in% colnames(this.X), ]
    est.e.prob <- eval_prob(this.X, this.pool_e_coef)
    est.y.coef <- pooled.y.coef[rownames(pooled.y.coef) %in% colnames(this.X.tilde), ]
  }
  
  this.X.tilde.treat <- this.X.tilde; this.X.tilde.treat[,'treat'] <- 1
  this.X.tilde.control <- this.X.tilde; this.X.tilde.control[,'treat'] <- 0
  this.y.treat.pred <- eval_prob(this.X.tilde.treat, est.y.coef)
  this.y.control.pred <- eval_prob(this.X.tilde.control, est.y.coef)
  
  if (estimand == "ATT") {
    influence.function <- this.treat*(this.y-this.y.control.pred)-est.e.prob*(1-this.treat)/(1-est.e.prob)*(this.y-this.y.control.pred)
  } else {
    influence.function <- (this.y.treat.pred-this.y.control.pred)+this.treat/est.e.prob*(this.y-this.y.treat.pred)-(1-this.treat)/(1-est.e.prob)*(this.y-this.y.control.pred)
  }
  
  est.tau <- mean(influence.function)
  V <- mean((influence.function - est.tau)**2)/n_obs
  
  return(list(est=est.tau, V=V, influence.function=influence.function))
}
