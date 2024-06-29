# Overview
This repository contains revised code for the Power-CHORD project. When we release an updated preprint, the contents of this repository will either be merged with the public repository or replace the public repository.

## Initial setup of dependencies 
This project uses `renv` to manage dependencies, depending on your `R` configuration, it may be necessary to call `renv::load()` when you first clone this repository. You can run the following to confirm that your `R` environment is aware of your project library.
```R
.libPaths()
```
If you see a library path that is contained within your repo, this likely means that `renv` is working properly. 