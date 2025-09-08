# concurrentFLMM: Code for generating data and figures for Xin et al. (2025)

## Overview

Each `[NUMBER]-[NAME].Rmd` file in the base directory contains the code needed to generate the corresponding figure. We also include the directories

- `data/`: datasets for generating the figures,
- `fun/`: supporting functions, 
- `img/`: image outputs, with subfolders named according to the `.Rmd` files,
- `pdf/`: knitted PDFs of each R markdown file, and
- `svg/`: formatted figures as SVG files.

To load package dependencies, see the file `0-setup.Rmd`. 

## File contents

- `0-setup.Rmd`: Loading dependencies with `renv` and `devtools`.
- `1-cflmm.Rmd`: Visualization of an example concurrent functional linear mixed model on synthetic data.
- `2-lick.Rmd`: Comparison of concurrent models to non-concurrent models when summarizing the behavioral covariate of licking data in mice. 
- `3-d2pvt.Rmd`: Comparisons of coefficient estimates for an experiment where mice receive different rewards at variable times.
