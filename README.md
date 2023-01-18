## Federated Causal Inference

This repository contains the code for "Federated causal inference in heterogeneous observational data".

We are interested in estimating the effect of a treatment applied to individuals at multiple sites, where data is stored locally for each site. Due to privacy constraints, individual-level data cannot be shared across sites; the sites may also have heterogeneous populations and treatment assignment mechanisms. 


See "code/Federation Documentation.Rmd" for the example to use our federated methods to draw inference on the average treatment effects of combined data across sites. Our methods first compute summary statistics locally using propensity scores and then aggregate these statistics across sites to obtain point and variance estimators of average treatment effects.


### Reference

```
@article{xiong2021federated,
  title={Federated causal inference in heterogeneous observational data},
  author={Xiong, Ruoxuan and Koenecke, Allison and Powell, Michael and Shen, Zhu and Vogelstein, Joshua T and Athey, Susan},
  journal={arXiv preprint arXiv:2107.11732},
  year={2021}
}
```
