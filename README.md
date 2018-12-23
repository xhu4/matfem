# matfem - a MATlab Finite Element Method package

[![GitHub license](https://img.shields.io/apm/l/vim-mode.svg)](https://github.com/xhu4/matfem/blob/master/LICENSE)
![matlab](https://img.shields.io/badge/language-Matlab-blue.svg)
[![Github all releases](https://img.shields.io/github/downloads/xhu4/matfem/total.svg)](https://GitHub.com/xhu4/matfem/releases/)

An object-oriented highly vectorized matlab package for solving finite element models. Only support 2D problems.

This package is a (object oriented & vectorized) reimplementation of my homework solution for a finite element class. This code is developed to help me understand FEM, and is
* **good for**: numerical analysis, try ideas fast and fail fast in small model,
* **not good for**: production code, or any analysis that FEM algorithm itself is not of interest.


## Benefits

* Easier to code than Deal.II
* Lower level than FEniCS
* Faster than coding yourself
* Cheaper than Matlab PDE toolbox

## Limitations

* 2D only
* At most degree 2 Lagrange element
* Even though it's highly vectorized, it's not very fast
* Does not run in parallel

## Examples

All `.m` files in this root folder are example codes.
Need to write examples in matlab notebook format.

## ToDo List

- [ ] Add documentation for package code :scream:
- [ ] `.m` -> `.mlx`
