# Master Thesis – Simulation Codes

This repository contains the MATLAB codes used for the numerical simulations and bifurcation analysis presented in the Master's thesis together with the summary table for the models review:

**From Negative Plant--Soil Feedback to Vegetation Patterns: A Cross-Diffusion Approach**

by **Sara Bonacina**
MSc in Mathematics, University of Trento
2026

## Contents

The repository contains the following MATLAB codes:

* `1Dsim_ode15s.m` – code for the one-dimensional numerical simulations.
* `evol2d_4scenarios.m` – code for the two-dimensional numerical simulations.
* `simtheta_1d_video.m` – implementation of the different transition functions $\theta$ considered in the one-dimensional simulations.
* `continuation/` – files used for the numerical continuation and bifurcation analysis of stationary solutions with `pde2path`.

In addition, the repository includes:

* `modelsNPSF_review.pdf` – table summarizing all the models considered in the thesis with a link to the corresponding paper for each model.
* `modelsNPSF_review.xlsx` 

## Software

The simulations were performed using **MATLAB R2024b**.

The bifurcation analysis was carried out using **pde2path**.

## Reproducibility

The codes contain the model parameters, initial conditions, and numerical settings used to obtain the results presented in the thesis.

The numerical results obtained from these codes may depend on the MATLAB version and, for the bifurcation analysis, on the version of `pde2path` used.

## Citation

If you use these codes, please cite the associated Master's thesis:

> S. Bonacina, *From Negative Plant--Soil Feedback to Vegetation Patterns: A Cross-Diffusion Approach*, Master's thesis, University of Trento, 2026.

Repository:

> https://github.com/Sara-Bonacina/MasterThesis-Materials
