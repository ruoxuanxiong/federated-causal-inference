
## get the all covariates names (union of dataset-specific covariates names)
## input args: full_list<list>: a list of dataset-specific estimated coefficients 
##           full_matrix<list>: a list of dataset-specific matrices; these can be hessian, outer product of gradient of log likelihood function, etc. 
##  
## return: full_covar_names<vector>: a vector of covariates names; a union of dataset-specific covariates  
##         full_covar<list<vector>>: each element in this list is a vector of dataset-specific covariate names  
get_full_covar <- function(full_list=NULL, full_matrix=NULL) {
  full_covar <- list()
  if (is.null(full_matrix)){
    for (i in c(1:length(full_list))) {
      full_covar[[i]] <- names(full_list[[i]])
    }
  } else {
    for (i in c(1:length(full_matrix))) {
      full_covar[[i]] <- colnames(full_matrix[[i]])
    }
  }
  full_covar_names <- full_covar[[1]]
  if (length(full_covar) > 1) {
    for (i in c(2:length(full_covar))) {
      full_covar_names <- union(full_covar_names, full_covar[[i]])
    }
  }
  return(list(full_covar_names=full_covar_names, full_covar=full_covar))
}

## get the federated MLE point estimates 
## input args: full_coef<list>: a list of dataset-specific MLE coefficients 
##          full_hessian<list>: a list of dataset-specific hessian matrices 
## 
## return: a vector of federated MLE point estimates 
pool_mle_coef <- function(full_coef, full_hessian) {
  
  full_covar_out <- get_full_covar(full_list = full_coef)
  full_covar_names <- full_covar_out$full_covar_names; full_covar <- full_covar_out$full_covar
  
  hess_sum <- matrix(0, nrow = length(full_covar_names), ncol = length(full_covar_names))
  row.names(hess_sum) <- full_covar_names
  colnames(hess_sum) <- full_covar_names
  
  hess_coeff_sum <- matrix(0, nrow = length(full_covar_names), ncol=1)
  row.names(hess_coeff_sum) <- full_covar_names
  
  for (i in c(1:length(full_hessian))) {
    hess_sum[full_covar[[i]], full_covar[[i]]] <- hess_sum[full_covar[[i]], full_covar[[i]]] + full_hessian[[i]]
    hess_coeff_sum[full_covar[[i]],] <- hess_coeff_sum[full_covar[[i]],] + full_hessian[[i]] %*% full_coef[[i]]
  }
  pool_coef <- solve(hess_sum) %*% hess_coeff_sum
  return(pool_coef)
}

## get the federated MLE covariance estimates 
## input args: full_bread<list>: a list of dataset-specific Fisher information matrices (bread) 
##              full_meat<list>: a list of dataset-specific outer product of gradient of log likelihood function matrices (meat)
##       model_misspec<boolean>: indicator if robust variance formulation is used; default is TRUE
## 
## return: a matrix of federated MLE covariance estimates 
pool_mle_cov <- function(full_bread, full_meat, model_misspec=TRUE) {
  full_covar_out <- get_full_covar(full_matrix = full_bread)
  full_covar_names <- full_covar_out$full_covar_names; full_covar <- full_covar_out$full_covar
  
  pool_bread <- matrix(0, nrow = length(full_covar_names), ncol = length(full_covar_names))
  pool_meat <- matrix(0, nrow = length(full_covar_names), ncol = length(full_covar_names))
  row.names(pool_bread) <- full_covar_names; colnames(pool_bread) <- full_covar_names
  row.names(pool_meat) <- full_covar_names; colnames(pool_meat) <- full_covar_names
  
  for (i in c(1:length(full_bread))) {
    pool_bread[full_covar[[i]], full_covar[[i]]] <- pool_bread[full_covar[[i]], full_covar[[i]]] + full_bread[[i]]
    pool_meat[full_covar[[i]], full_covar[[i]]] <- pool_meat[full_covar[[i]], full_covar[[i]]] + full_meat[[i]]
  }
  
  if (model_misspec) {
    pool_cov <- solve(pool_bread) %*% pool_meat %*% solve(pool_bread)
  } else {
    pool_cov <- solve(-pool_bread)
  }
  return(pool_cov)
}

## get the federated IPW-MLE covariance estimates 
## input args: full_ht_cov_out<list>: a list of IPW-MLE output by calc_ht_cov/calc_ht_cov_simple function 
##                pool_e_cov<matrix>: the variance covariance matrix of propensity score estimation; NULL if propensity score is known
##     estimated_propensity<boolean>: if propensity scores are estimated or known; default is TRUE
## 
## return: a matrix of federated IPW-MLE covariance estimates 
pool_ht_cov <- function(full_ht_cov_out, pool_e_cov=NULL, estimated_propensity=TRUE) {
  n_obs_total <- 0
  full_A0 <- list(); 
  
  for (i in c(1:length(full_ht_cov_out))) {
    full_A0[[i]] <- full_ht_cov_out[[i]]$A0
    n_obs_total <- n_obs_total + full_ht_cov_out[[i]]$n_obs
  }
  full_covar_out <- get_full_covar(full_matrix = full_A0)
  full_covar_names <- full_covar_out$full_covar_names; full_covar <- full_covar_out$full_covar

  pool_A0 <- matrix(0, nrow = length(full_covar_names), ncol = length(full_covar_names))
  pool_D0 <- matrix(0, nrow = length(full_covar_names), ncol = length(full_covar_names))
  row.names(pool_A0) <- full_covar_names; colnames(pool_A0) <- full_covar_names
  row.names(pool_D0) <- full_covar_names; colnames(pool_D0) <- full_covar_names
  
  for (i in c(1:length(full_ht_cov_out))) {
    this.n_obs <- full_ht_cov_out[[i]]$n_obs
    pool_A0[full_covar[[i]], full_covar[[i]]] <- pool_A0[full_covar[[i]], full_covar[[i]]] + full_ht_cov_out[[i]]$A0 * this.n_obs
    pool_D0[full_covar[[i]], full_covar[[i]]] <- pool_D0[full_covar[[i]], full_covar[[i]]] + full_ht_cov_out[[i]]$D0 * this.n_obs
  }
  
  pool_A0 <- pool_A0/n_obs_total
  pool_D0 <- pool_D0/n_obs_total
  pool_B0 <- pool_D0
  
  if (estimated_propensity) {
    full_C1 <- list()
    full_C2 <- list()
    for (i in c(1:length(full_ht_cov_out))) {
      full_C1[[i]] <- full_ht_cov_out[[i]]$C1
      full_C2[[i]] <- full_ht_cov_out[[i]]$C2
    }
    full_e_covar_out_1 <- get_full_covar(full_matrix = full_C1)
    full_e_covar_out_2 <- get_full_covar(full_matrix = full_C2)
    full_e_covar_names <- full_e_covar_out_1$full_covar_names; full_e_covar <- full_e_covar_out_1$full_covar
    pool_C1 <- matrix(0, nrow = length(full_covar_names), ncol = length(full_e_covar_names))
    row.names(pool_C1) <- full_covar_names; colnames(pool_C1) <- full_e_covar_names
    pool_C2 <- matrix(0, nrow = length(full_covar_names), ncol = length(full_e_covar_names))
    row.names(pool_C2) <- full_covar_names; colnames(pool_C2) <- full_e_covar_names
    
    for (i in c(1:length(full_ht_cov_out))) {
      pool_C1[full_covar[[i]], full_e_covar[[i]]] <- pool_C1[full_covar[[i]], full_e_covar[[i]]] + full_ht_cov_out[[i]]$C1 * this.n_obs
      pool_C2[full_covar[[i]], full_e_covar[[i]]] <- pool_C2[full_covar[[i]], full_e_covar[[i]]] + full_ht_cov_out[[i]]$C2 * this.n_obs
    }
    pool_C1 <- pool_C1/n_obs_total
    pool_C2 <- pool_C2/n_obs_total
    pool_M0.inv <- n_obs_total * pool_e_cov
    pool_B0 <- pool_D0 - pool_C1 %*% pool_M0.inv %*% t(pool_C2) - pool_C2 %*% pool_M0.inv %*% t(pool_C1) + pool_C2 %*% pool_M0.inv %*% t(pool_C2)
  }
  pool_A0.inv <- solve(pool_A0)
  pool_cov <- pool_A0.inv %*% pool_B0 %*% pool_A0.inv / n_obs_total
  
  return(pool_cov)
}

## get the federated MLE estimation results
##        input args: full_y<list>: a list of dataset-specific binary outcomes 
##                    full_X<list>: a list of dataset-specific design matrices 
##               unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##         unstable_covars<vector>: if models are unstable, a vector of unstable covariate names 
##               num_datasets<int>: number of datasets
##          full_reg.formula<list>: list of all outcome regression formulas
##          model_misspec<boolean>: indicator if robust variance formulation is used; default is TRUE
##
## return: a summary table of federated MLE estimation results  
poolMLE <- function(full_y, full_X, unstable=FALSE, unstable_covars=NULL, num_datasets=NULL, 
                    full_reg.formula=NULL, model_misspec=TRUE) {
  
  full_coef <- list(); full_outer_grad <- list(); full_hessian <- list(); full_V <- list()
  
  if (unstable) {
    full_reg.formula.unstable <- UnstableRegModels(full_reg.formula, unstable_covars, num_datasets)
  }
  
  for (i in c(1:length(full_y)) ) {
    this.y <- full_y[[i]]; this.X <- full_X[[i]]
    
    if (!unstable) {
      out_mle <- est_mle(this.y, this.X, robust=FALSE, dataset_idx=i, full_reg.formula=full_reg.formula) 
    } else {
      out_mle <- est_mle(this.y, this.X, robust=FALSE, dataset_idx=i, 
                         unstable=TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets, 
                         full_reg.formula.unstable=full_reg.formula.unstable)  
    }
    full_coef[[i]] <- out_mle$est
    full_hessian[[i]] <- out_mle$hessian
    full_V[[i]] <- out_mle$V
  }
  
  pool_coef <- pool_mle_coef(full_coef, full_hessian)
  
  for (i in c(1:length(full_y))) {
    this.y <- full_y[[i]]
    this.X <- full_X[[i]]
    this.w <- rep(1, length(this.y))
    
    if (!unstable) {
      this.pool_coef <- pool_coef[colnames(this.X), ]
      out_grad_hessian <- eval_grad_hessian(this.pool_coef, this.y, this.X, this.w, robust=model_misspec)
    } else {
      this.pool_coef <- pool_coef[rownames(pool_coef) %in% names(full_coef[[i]]), ]
      out_grad_hessian <- eval_grad_hessian(this.pool_coef, this.y, this.X, this.w, dataset_idx = i, robust=model_misspec,
                                            unstable = TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets)
    }
    full_outer_grad[[i]] <- out_grad_hessian$outer.grad; full_hessian[[i]] <- out_grad_hessian$hessian
  }
  pool_cov <- pool_mle_cov(full_hessian, full_outer_grad, model_misspec=model_misspec)
  summary_table <- make_summary_table(pool_coef, pool_cov)
  return(summary_table)
}

## get the federated IPW-MLE estimation results
##              input args: full_y<list>: a list of dataset-specific binary outcomes 
##                    full_X.tilde<list>: a list of dataset-specific design matrices (include column of treatment)
##                          full_X<list>: a list of dataset-specific design matrices 
##                      full_treat<list>: a list of dataset-specific treatment assignments
##                     unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##               unstable_covars<vector>: if models are unstable, a vector of unstable covariate names 
##                     num_datasets<int>: number of datasets
##                full_reg.formula<list>: list of all outcome regression formulas
##          full_treat.reg.formula<list>: list of all treatment regression formulas
##                        estimand<char>: "ATE" or "ATT"; default is "ATE"
##         estimated_propensity<boolean>: if propensity scores are estimated or known; default is TRUE
##                model_misspec<boolean>: indicator if robust variance formulation is used; default is TRUE
##
## return: a summary table of federated IPW-MLE estimation results  
poolHT <- function(full_y, full_X.tilde, full_X, full_treat, 
                   unstable=FALSE, unstable_covars=NULL,num_datasets=NULL, 
                   full_reg.formula=NULL, full_treat.reg.formula=NULL, 
                   estimand="ATE", estimated_propensity=TRUE, model_misspec=TRUE) {
  
  full_e_coef <- list(); full_e_grad <- list(); full_e_outer_grad <- list(); full_e_hessian <- list()
  
  if (unstable) {
    full_reg.formula.unstable <- UnstableRegModels(full_reg.formula, unstable_covars, num_datasets)
    full_treat.reg.formula.unstable <- UnstableRegModels(full_treat.reg.formula, unstable_covars, num_datasets)
  }
  
  for (i in c(1:length(full_y)) ) {
    this.treat <- full_treat[[i]]; this.X <- full_X[[i]]
    if (!unstable) {
      e_out_mle <- est_mle(this.treat, this.X, robust=FALSE, dataset_idx=i, treat.reg=TRUE, full_treat.reg.formula=full_treat.reg.formula)
    } else {
      e_out_mle <- est_mle(this.treat, this.X, robust=FALSE, dataset_idx=i, treat.reg=TRUE, unstable=TRUE, 
                           unstable_covars=unstable_covars, num_datasets=num_datasets, 
                           full_treat.reg.formula.unstable=full_treat.reg.formula.unstable)  
    }
    full_e_coef[[i]] <- e_out_mle$est
    full_e_hessian[[i]] <- e_out_mle$hessian
  }
  pool_e_coef <- pool_mle_coef(full_e_coef, full_e_hessian)
  
  for (i in c(1:length(full_y))) {
    this.treat <- full_treat[[i]]; this.X <- full_X[[i]]; this.w <- rep(1, length(this.treat))
    
    if(!unstable) {
      this.pool_e_coef <- pool_e_coef[colnames(this.X),]
      e_out_grad_hessian <- eval_grad_hessian(this.pool_e_coef, this.treat, this.X, this.w)
    } else {
      this.pool_e_coef <- pool_e_coef[rownames(pool_e_coef) %in% names(full_e_coef[[i]]), ]
      e_out_grad_hessian <- eval_grad_hessian(this.pool_e_coef, this.treat, this.X, this.w, dataset_idx=i,
                                              unstable=TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets)
    }
    full_e_grad[[i]] <- e_out_grad_hessian$grad
    full_e_outer_grad[[i]] <- e_out_grad_hessian$outer.grad
    full_e_hessian[[i]] <- e_out_grad_hessian$hessian
  }
  
  pool_e_cov <- pool_mle_cov(full_e_hessian, full_e_outer_grad, model_misspec=model_misspec)
  full_w.y <- list(); full_est.e.prob <- list()
  full_w.alt.y <- list()
  for (i in c(1:length(full_y))) {
    this.treat <- full_treat[[i]]; this.X <- full_X[[i]]
    
    if (!unstable) {
      this.pool_e_coef <- pool_e_coef[colnames(this.X),]
      this.est.e.prob <- eval_prob(this.X, this.pool_e_coef)
    } else {
      this.pool_e_coef <- pool_e_coef[rownames(pool_e_coef) %in% names(full_e_coef[[i]]), ]
      this.X.interact <- AddInteractCols(this.X, dataset_idx = i, unstable_covars=unstable_covars, num_datasets=num_datasets)
      this.X.interact <- this.X.interact[,names(this.pool_e_coef)]
      this.est.e.prob <- eval_prob(this.X.interact, this.pool_e_coef)
    }
    full_est.e.prob[[i]] <- this.est.e.prob
    full_w.y[[i]] <- propensity_weight(this.treat, this.est.e.prob, estimand = estimand)
    full_w.alt.y[[i]] <- alternative_propensity_weight(this.treat, this.est.e.prob, estimand = estimand)
  }
  
  full_y_coef <- list(); full_y_grad <- list(); full_y_outer_grad <- list(); full_y_hessian <- list(); full_y_grad.alt <- list()
  for (i in c(1:length(full_y)) ) {
    this.y <- full_y[[i]]; this.X.tilde <- full_X.tilde[[i]]; this.w.y <- full_w.y[[i]]
    if (!unstable) {
      y_mle_out <- est_mle(this.y, this.X.tilde, this.w.y, robust=FALSE,  wsq=TRUE, dataset_idx=i, full_reg.formula=full_reg.formula)
    } else {
      y_mle_out <- est_mle(this.y, this.X.tilde, this.w.y, robust=FALSE, wsq=TRUE, 
                           dataset_idx=i, unstable=TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets, 
                           full_reg.formula.unstable=full_reg.formula.unstable)  
    }

    full_y_coef[[i]] <- y_mle_out$est
    full_y_grad[[i]] <- y_mle_out$grad
    full_y_outer_grad[[i]] <- y_mle_out$outer.grad
    full_y_hessian[[i]] <- y_mle_out$hessian
  }
  pool_y_coef <- pool_mle_coef(full_y_coef, full_y_hessian)
  
  for (i in c(1:length(full_y))) {
    this.y <- full_y[[i]]; this.X.tilde <- full_X.tilde[[i]]; this.w.y <- full_w.y[[i]]
    if (!unstable) {
      this.pool_y_coef <- pool_y_coef[colnames(this.X.tilde),]
      y_out_grad_hessian <- eval_grad_hessian(this.pool_y_coef, this.y, this.X.tilde, this.w.y, wsq=TRUE)
    } else {
      this.pool_y_coef <- pool_y_coef[rownames(pool_y_coef) %in% names(full_y_coef[[i]]), ]
      y_out_grad_hessian <- eval_grad_hessian(this.pool_y_coef, this.y, this.X.tilde, this.w.y, wsq=TRUE, dataset_idx=i, 
                                              unstable=TRUE, unstable_covars=unstable_covars, num_datasets=num_datasets)
    }
    full_y_grad[[i]] <- y_out_grad_hessian$grad
    full_y_outer_grad[[i]] <- y_out_grad_hessian$outer.grad
    full_y_hessian[[i]] <- y_out_grad_hessian$hessian
    full_y_grad.alt[[i]] <- diag(as.vector(full_w.alt.y[[i]]) / as.vector(full_w.y[[i]])) %*% y_out_grad_hessian$grad
  }
  
  full_ht_cov_out <- list()
  for (i in c(1:length(full_y))) {
    n_obs <- length(full_y[[i]])
    if (estimated_propensity) {
      full_ht_cov_out[[i]] <- calc_ht_cov(full_e_grad[[i]], full_e_outer_grad[[i]], full_e_hessian[[i]], 
                                          full_y_grad[[i]], full_y_grad.alt[[i]], full_y_outer_grad[[i]], full_y_hessian[[i]], n_obs, unstable = unstable)
    } else {
      full_ht_cov_out[[i]] <- calc_ht_cov_simple(full_y_grad[[i]], full_y_outer_grad[[i]], full_y_hessian[[i]], n_obs)
    }
  }
  pool_y_cov <- pool_ht_cov(full_ht_cov_out, pool_e_cov, estimated_propensity=estimated_propensity)
  summary_table <- make_summary_table(pool_y_coef, pool_y_cov)
  return(summary_table)
}

## get the federated AIPW estimation results
##              input args: full_y<list>: a list of dataset-specific binary outcomes 
##                    full_X.tilde<list>: a list of dataset-specific design matrices (include column of treatment)
##                          full_X<list>: a list of dataset-specific design matrices 
##                      full_treat<list>: a list of dataset-specific treatment assignments
##                     unstable<boolean>: indicator if the models are unstable (Condition 4 and/or 5); default is FALSE
##                full_reg.formula<list>: list of all outcome regression formulas
##          full_treat.reg.formula<list>: list of all treatment regression formulas
##                        estimand<char>: "ATE" or "ATT"; default is "ATE"
##
## return: a summary table of federated AIPW estimation results  
poolAIPW <- function(full_y, full_X.tilde, full_X, full_treat, unstable=FALSE, 
                     full_reg.formula=NULL, full_treat.reg.formula=NULL, estimand="ATE") {
  
  full_aipw_out <- list(); full_n_obs <- list(); n_obs_total <- 0
  if (!unstable) {
    full_e_coef <- list(); full_e_hessian <- list()
    full_y_coef <- list(); full_y_hessian <- list()
    
    for (ix in c(1:length(full_y)) ) {
      this.treat <- full_treat[[ix]]; this.X <- full_X[[ix]]; this.y <- full_y[[ix]]; this.X.tilde <- full_X.tilde[[ix]]
      e_out_mle <- est_mle(this.treat, this.X, robust=FALSE, dataset_idx=ix, treat.reg=TRUE, full_treat.reg.formula=full_treat.reg.formula)
      y_out_mle <- est_mle(this.y, this.X.tilde, robust=FALSE, dataset_idx=ix, full_reg.formula=full_reg.formula)
      full_e_coef[[ix]] <- e_out_mle$est
      full_e_hessian[[ix]] <- e_out_mle$hessian
      full_y_coef[[ix]] <- y_out_mle$est
      full_y_hessian[[ix]] <- y_out_mle$hessian
    }
    pooled.y.coef <- pool_mle_coef(full_y_coef, full_y_hessian)
    pooled.e.coef <- pool_mle_coef(full_e_coef, full_e_hessian)
  }
  
  for (i in c(1:length(full_y))) {
    this.y <- full_y[[i]]; this.X.tilde <- full_X.tilde[[i]]; this.X <- full_X[[i]]; this.treat <- full_treat[[i]]
    treat.reg.formula <- full_treat.reg.formula[[i]]; reg.formula <- full_reg.formula[[i]]
    if (unstable) {
      aipw_out <- est_aipw(this.y, this.X.tilde, this.X, this.treat, estimand=estimand, 
                           treat.reg.formula=treat.reg.formula, reg.formula=reg.formula)
    } else {
      aipw_out <- est_aipw(this.y, this.X.tilde, this.X, this.treat, estimand=estimand, 
                           pooled.y.coef=pooled.y.coef, pooled.e.coef=pooled.e.coef)
    }
    full_aipw_out[[i]] <- aipw_out
    full_n_obs[[i]] <- length(this.y)
    n_obs_total <- n_obs_total + full_n_obs[[i]]
  }
  
  pool_coef <- 0
  pool_cov <- 0
  
  if (unstable) {
    for (i in c(1:length(full_y))) {
      pool_coef <- pool_coef + full_n_obs[[i]]/n_obs_total * full_aipw_out[[i]]$est
      pool_cov <- pool_cov + full_n_obs[[i]]**2/n_obs_total * full_aipw_out[[i]]$V
    }
    pool_cov <- pool_cov / n_obs_total
  } else {
    for (i in c(1:length(full_y))) {
      pool_cov <- pool_cov + 1/full_aipw_out[[i]]$V
    }
    pool_cov <- 1/pool_cov
    for (i in c(1:length(full_y))) {
      pool_coef <- pool_coef + pool_cov * full_aipw_out[[i]]$est / full_aipw_out[[i]]$V
    }
    pool_cov <- pool_cov 
  }
  pool_cov <- as.matrix(pool_cov)
  summary_table <- make_summary_table(pool_coef, pool_cov)
  return(summary_table)
}



