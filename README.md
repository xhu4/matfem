# matfem - a MATlab Finite Element Method package

[![GitHub license](https://img.shields.io/apm/l/vim-mode.svg)](https://github.com/xhu4/matfem/blob/master/LICENSE)
![matlab](https://img.shields.io/badge/language-Matlab-blue.svg)
[![Github all releases](https://img.shields.io/github/downloads/xhu4/matfem/total.svg)](https://GitHub.com/xhu4/matfem/releases/)

An object-oriented highly vectorized matlab package for solving finite element
models. Only support 2D problems. Developed in 2016b. Have not tested in any
other Matlab version.

![](demos/coupled.gif)

This package is a (object oriented & vectorized) reimplementation of my homework
solution for a finite element class. This code is developed to help me
understand FEM, and is

- **good for**: numerical analysis, try ideas fast and fail fast in small model,
- **not good for**: production code, or any analysis that FEM algorithm itself
  is not of interest.

## Benefits

- Easier to code than Deal.II
- Lower level than FEniCS
- Faster than coding yourself
- Cheaper than Matlab PDE toolbox

## Limitations

- 2D only
- At most degree 2 Lagrange element
- Even though it's highly vectorized, it's not very fast
- Does not run in parallel

## Examples (Under Construction :construction:)

All `.xlm` files in the `demos` folder are example codes.
The hello world `a_poisson.xlm` is a good place to start.

Matlab Live Codes:
1. [Poisson Equation](demos/a_poisson.mlx)
2. [Poisson Equation with Neumann Boundaries](demos/b_poisson_neumann.mlx)

Matlab Script (Not documented yet):
1. [Poisson Equation with Robin Boundaries](demos/exmp3.m)
2. [Poisson Equation with Matrix-valued Coefficient](demos/final1.m)
3. [Elasticity Equation](demos/elasticity.m)
4. [Steady-State Stokes Equation](demos/stokesSteady.m)
5. [Stokes Equation (time dependent)](demos/stokes.m)
6. [Coupled Dual-Porosity-Stokes Model](demos/coupled.m)

## Roadmap

- [ ] Add more detailed documentation for package code :memo:
- [ ] Translate matlab script into well documented live code.
