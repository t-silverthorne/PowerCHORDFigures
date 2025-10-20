# Overview
This repository contains the code necessary for generating the figures in our manuscript. The actual PowerCHORD repository, which contains the optimization methods and power analysis tools can be found [here](https://github.com/t-silverthorne/PowerCHORD). By default, it is loaded as a git submodule when you clone this repo. 


To recreate the figures from the manuscript, please consult the `vector_figures/` directory.


## Dependencies 
This project uses `renv` to manage external dependencies. Depending on your `R` configuration, it may be necessary to call `renv::load()` when you first clone this repository. You can run the following to confirm that your `R` environment is aware of your project library.
```R
.libPaths()
```
If you see a library path that is contained within your current working directory (i.e. your clone of this repository), then `renv` is aware of your local library. You can install external dependencies using the following command.
```R
renv::hydrate()
```


## PowerCHORD setup
Only the basic `R` functions in PowerCHORD are required for recreating the figures. To load these functions, simply ensure that the git submodule is initialized
```bash
git submodule init
git submodule update --remote
```
To make sure `renv` is aware of PowerCHORD, run the following within an R session:
```R
renv::install("./PowerCHORD")
```

## Git-lfs setup
Some of the figures require data from large simulations. You can download this data as follows:
```bash
git lfs install
git lfs pull
```
 
If you want to explore the code on a deeper level and solve your own optimization problems, you will need to follow the full installation instructions in the [PowerCHORD repo](https://github.com/t-silverthorne/PowerCHORD).
