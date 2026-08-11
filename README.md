# Supplementary files for "Microbial public goods dynamics create a competition-colonization trade-off that leads to the evolution of coexisting strategies"

## Contents

### mathematica/
Wolfram Mathematica notebooks (compiled in version 14.0 on Windows 11).

Contents:

- 1-in-patch-model.nb: In-patch model solution, Figure 1B-E
- 2-metapopulation-model.nb: Figure 3, Figure S1, system solution for D6 and D7
- AppendixB.nb: Analysis of the unreduced in-patch model
- AppendixE.nb: Analysis of the lactamase model
- AppendixF.nb: Analysis of the invertase model

### simulations/
R and Rust source code used to run simulations and generate figures.

Contents:

R/patch_occupancy_sim/
- rust/: simulation source code; see readme therein
- figure3.rmd: Figure 3
- lactamase.rmd: Figure 5D
- invertase.rmd: Figure 5H
- stats.rmd: Figure D2 and D3
- sigma_sweep.rmd: Figure D4
- mu_sweep.rmd: Figure D5


R/gillespie_sim/
- rust/: simulation source code; see readme therein
- nomut.rmd: Figures D6 and D7
- example.rmd: Figure D8
- stats.rmd: Figure D9

R/kick_flow_sim/
- rust/: simulation source code; see readme therein
- sweeps.rmd: Figures 4A and D10

### isolate_comparison/
R code used to generate Figure 5B and 5C.