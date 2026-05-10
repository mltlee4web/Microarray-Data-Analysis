# Microarray-Data-Analysis

A collection of R scripts and functions for performing permutation-based differential expression analysis on microarray data for the paper: "Lee ML, Gray RJ, Björkbacka H, Freeman MW. Generalized rank tests for replicated microarray data. Stat Appl Genet Mol Biol. 2005;4:Article3. doi:10.2202/1544-6115.1093". This repository provides non-parametric statistical tools to handle independent, paired, and stratified experimental designs using exact and Monte Carlo permutations.

## Overview

This repository contains modernized R implementations of rank-analysis algorithms. It is designed to be a lightweight, easily accessible tool for researchers needing robust permutation tests, false discovery rate (FDR) estimations, and handling of complex batch effects without relying on theoretical distribution assumptions.

## Repository Structure

* `R/permax.R`: The core statistical engine containing the `permax()`, `perm_test()`, and `plot.expord()` functions.
* `R/Rank_Analysis_of_Microarrays.R`: Higher-level wrapper functions and analysis pipelines dependent on the `permax` engine.
* `data/`: Example datasets for demonstrating paired and stratified test configurations.
* `examples/`: Sample scripts demonstrating how to execute various experimental designs (e.g., independent groups, repeated measures).

## Getting Started

Because this is a research repository and not a CRAN package, you do not need to install it via standard package managers. You can load the functions directly into your R session by sourcing the scripts from this repository.

### Prerequisites
* R (version 3.5.0 or higher recommended)