# tas: Forbidden Knowledge and Specialized Training

This repository contains the code used to perform the analyses in the paper:

Rohlfs, Chris, 2023. "Forbidden Knowledge and Specialized Training: A Versatile Solution for the Two Main Sources of Overfitting in Linear Regression," *The American Statistician*, 77(2): 160-8.  
[https://doi.org/10.1080/00031305.2022.2128874](https://doi.org/10.1080/00031305.2022.2128874)

## Overview

This repository includes scripts for both the empirical application and the Monte Carlo simulations presented in the paper. Below are descriptions of the key scripts, specifications, inputs, and outputs.

### Empirical Analysis

- **`baseline.empirical.R`**  
  Computes the baseline specification for the empirical application.

  **Key Specifications:**
  - **line 7**: `cores <- 20` (for multicore processing).
    - Parallel processing is used, so exact replication of the paper's results using `set.seed` is not guaranteed.
  - **line 177**: `pct <- 0.50` (determines fraction of the sample to be used for training).

  **Inputs:**
  - `mri.RDS`: Cleaned predictors and outcomes from the Neurocognitive Aging Data.

  **Outputs:**
  - `components.RDS`: Large (~0.6 GB) file containing a `data.table` of results by observation, iteration, and specification.
  - `rss.RDS`: Smaller `data.table` with summary statistics for each specification and iteration.
  - `MRI_50pct_rss.png`: Sample output graph; additional graphs are in `doc.graphs.R`.

### Monte Carlo Simulation

- **`baseline.simulation.R`**  
  Computes the baseline specification for the Monte Carlo simulation.

  **Specifications:**  
  Specification details are in lines 6 to 13.  
  For homoskedastic specifications, data generation and specification details differ slightly.  
  A sample version of this calculation is in `homosked45_125.R`.

  **Outputs:**
  - `simulation.data.hetero.txt`: Text version of components (similar to the empirical analysis). Saved as a text file to use the "append" function, given the large size (some versions can be as large as 25-30 GB).
  - `sims.RDS`: Smaller `data.table` with summary statistics for each specification and iteration.

- **`homosked45_125.R`**  
  Computes an alternative homoskedastic specification for the Monte Carlo simulation.

## Citing This Work

If you use this code or the results from the paper in your research, please cite the following publication:

```bibtex
@article{rohlfs2023forbidden,
  title={Forbidden Knowledge and Specialized Training: A Versatile Solution for the Two Main Sources of Overfitting in Linear Regression},
  author={Rohlfs, Chris},
  journal={The American Statistician},
  volume={77},
  number={2},
  pages={160--168},
  year={2023},
  publisher={Taylor \& Francis},
  doi={10.1080/00031305.2022.2128874}
}
```

- Note: also calculates some extra statistics related to MAPE (used in figures) and sigmas

doc.graphs - Generates all the graphs for the paper. Requires that the empirical & simulation scripts have been run with various parameter configurations.
