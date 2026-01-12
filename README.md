# On the Hierarchical Bayes justification of Empirical Bayes Confidence Intervals

## Overview 

Multi-level normal hierarchical models play an important role in developing statistical theory in multiparameter estimation for a wide range of applications. In this article, we propose a new reconciliation framework of the empirical Bayes (EB) and hierarchical Bayes approaches for interval estimation of random effects under a two-level normal model. Our framework shows that a second-order efficient EB confidence interval, with  EB coverage error of order $O(m^{-3/2})$, can also be viewed as a credible interval whose posterior coverage is close to the nominal level, provided a carefully chosen prior—referred to as a ‘matching prior’—is placed on the hyperparameters. While existing literature has examined matching priors that reconcile frequentist and Bayesian inference in various settings, this paper is the first to study matching priors with the goal of interval estimation of random effects in a two-level model. We obtain an area-dependent matching prior on the variance component that achieves a proper posterior under mild regularity conditions. The theoretical results in the paper are corroborated through a Monte Carlo simulation study and a real data analysis.


## Main Article

The full article is available at https://arxiv.org/abs/2511.13037.


## Repository Structure

```
Repo/
│
├── data_analysis/                          # Scripts for real baseball data analysis
│   ├── baseball.txt                        # Baseball dataset
│   ├── baseball M3 PC cm.R                 # Model M3 with common mean, producing length and PC (Posterior Coverage)
│   ├── baseball M4 PC.R                    # Model M4 with covariates, producing length and PC
│   ├── baseball S4 S5 PC cm sim.R          # Setups S4 and S5 with simulated response, producing PC
│   │
│   └── plots/                             # Plotting and visualization scripts for baseball data
│       ├── baseball_plot_ebc.R            # Plots for EBC (Empirical Bayes Coverage)
│       └── baseball_plot_pc.R             # Plots for PC (Posterior Coverage)
│
├── simulation/                             # Monte Carlo simulation scripts
│   ├── sim S1 AS unbal EBC.R                # Setup S1, AS and YL method, unbalanced case, producing EBC
│   ├── sim S1 AS unbal PC.R                 # Setup S1, AS and YL method, unbalanced case, producing PC
│   ├── sim S2 YL unbal cm EBC.R             # Setup S2, YL method, unbalanced and common mean case, producing EBC
│   ├── sim S2 YL unbal cm PC.R              # Setup S2, YL method, unbalanced and common mean case, producing PC
│   ├── sim S3 AS bal EBC.R                  # Setup S3, AS and YL method, balanced case, producing EBC
│   └── sim S3 AS bal PC.R                   # Setup S3, AS and YL method, balanced case, producing PC
│
├── supporting_functions/                  # Core statistical routines
│   ├── AS_bal_fun.R                        # AS method for balanced case
│   ├── AS_unbal_fun.R                      # AS method for unbalanced case
│   ├── YL_bal_fun.R                        # YL method for balanced case
│   ├── YL_unbal_fun.R                      # YL method for unbalanced case
│   └── YL_unbal_fun_cm.R                   # YL unbalanced with common mean case
│
└── README.md                              # This file
```

## Requirements

The code is written in R and requires the following packages:

```
install.packages(c(
  "psych",
  "MASS",
  "mvtnorm",
  "rstan",
  "ggplot2"
  "dplyr",
  "foreach",
  "doParallel",
  "tidyverse",
  "parallel",
  "mcprogress"
))
```

## Reproducing the Results

1. Clone the repository:

```bash
git clone https://github.com/asen123/matching_prior.git
```

## Run simulations 

Run codes provided in the "simulation" folder. Source codes are available in the folder "supporting_functions".


## Run real data analysis

Run codes provided in the "data_analysis" folder.


## Data

The real data, "baseball.txt", is provided in the "data_analysis" folder. Source codes are available in the folder "supporting_functions".


## Citation

If you use this code or data, please cite:

```mathematica
@article{,
  title={{On the Hierarchical Bayes justification of Empirical Bayes Confidence Intervals}},
  author={Sen, Aditi and Hirose, Masayo Y. and Lahiri, Partha},
  journal={arXiv preprint},
  volume={arXiv:2511.13037},
  year={2025},
  month={Nov},
  note={\url{https://doi.org/10.48550/arXiv.2511.13037}}
}
```

