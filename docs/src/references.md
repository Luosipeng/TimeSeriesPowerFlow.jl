# References

## How to Cite HyDistFlow.jl

If you use HyDistFlow.jl in your research, please use the following citation format:

```bibtex
@software{HyDistFlow2025,
  author = {Sipeng Luo,Tianyang Zhao,Zhaohong Bie},
  title = {HyDistFlow.jl: a Julia package for time series power flow analysis},
  year = {2025},
  url = {https://github.com/Luosipeng/HyDistFlow.jl}
}
```
## Referenced Code and Inspirations

HyDistFlow.jl was developed with inspiration from and reference to the following open source projects:

### MATPOWER
An open-source MATLAB package for power system simulation and optimization.

```bibtex
@article{zimmerman2011matpower,
  title = {MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education},
  author = {Zimmerman, Ray D and Murillo-S{\'a}nchez, Carlos E and Thomas, Robert J},
  journal = {IEEE Transactions on Power Systems},
  volume = {26},
  number = {1},
  pages = {12--19},
  year = {2011},
  doi = {10.1109/TPWRS.2010.2051168}
}
```
```bibtex
@software{zimmerman2024matpower,
  author = {Zimmerman, Ray D and Murillo-S{\'a}nchez, Carlos E},
  title = {MATPOWER},
  version = {8.0},
  year = {2024},
  url = {https://matpower.org},
  doi = {10.5281/zenodo.11212330}
}
```

### PandaPower
An easy to use open source tool for power system modeling, analysis and optimization with a high degree of automation.

```bibtex
@ARTICLE{pandapower.2018,
    author={L. Thurner and A. Scheidler and F. Sch{\"a}fer and J. Menke and J. Dollichon and F. Meier and S. Meinecke and M. Braun},
    journal={IEEE Transactions on Power Systems},
    title={pandapower — An Open-Source Python Tool for Convenient Modeling, Analysis, and Optimization of Electric Power Systems},
    year={2018},
    month={Nov},
    volume={33},
    number={6},
    pages={6510-6521},
    doi={10.1109/TPWRS.2018.2829021},
    ISSN={0885-8950}}
```

## External Packages

HyDistFlow.jl relies on several external packages. We acknowledge and thank the contributors of these packages:

### AMD.jl

A Julia interface to the AMD library of Amestoy, Davis and Duff.

```bibtex
@software{montoison2020amd,
  author = {Montoison, Alexis and Orban, Dominique and Siqueira, Abel S. and contributors},
  title = {AMD.jl: A Julia interface to the AMD library of Amestoy, Davis and Duff},
  year = {2020},
  month = {5},
  day = {3},
  url = {https://github.com/JuliaSmoothOptimizers/AMD.jl},
  doi = {10.5281/zenodo.3381898},
  version = {0.4.0}
}
```

### DataFrames.jl

Tools for working with tabular data in Julia.

```bibtex
@article{JSSv107i04,
 title={DataFrames.jl: Flexible and Fast Tabular Data in Julia},
 volume={107},
 url={https://www.jstatsoft.org/index.php/jss/article/view/v107i04},
 doi={10.18637/jss.v107.i04},
 abstract={DataFrames.jl is a package written for and in the Julia language offering flexible and efficient handling of tabular data sets in memory. Thanks to Julia’s unique strengths, it provides an appealing set of features: Rich support for standard data processing tasks and excellent flexibility and efficiency for more advanced and non-standard operations. We present the fundamental design of the package and how it compares with implementations of data frames in other languages, its main features, performance, and possible extensions. We conclude with a practical illustration of typical data processing operations.},
 number={4},
 journal={Journal of Statistical Software},
 author={Bouchet-Valat, Milan and Kamiński, Bogumił},
 year={2023},
 pages={1--32}
}
```

### Plots.jl

A powerful visualization interface and toolkit in Julia.

```bibtex
@article{PlotsJL,
  doi = {https://doi.org/10.5334/jors.431},
  url = {https://openresearchsoftware.metajnl.com/articles/10.5334/jors.431/},
  author = {Christ, Simon and Schwabeneder, Daniel and Rackauckas, Christopher and Borregaard, Michael Krabbe and Breloff, Thomas},
  keywords = {Graphics (cs.GR), FOS: Computer and information sciences, FOS: Computer and information sciences, I.3.3},
  title = {Plots.jl -- a user extendable plotting API for the julia programming language},
  publisher = {Journal of Open Research Software},
  year = {2023},
  copyright = {Creative Commons Attribution 4.0 International}
}
```

### JuMP.jl
A Powerful optimization interface and toolkit in Julia

```bibtex
@article{Lubin2023,
    author = {Miles Lubin and Oscar Dowson and Joaquim {Dias Garcia} and Joey Huchette and Beno{\^i}t Legat and Juan Pablo Vielma},
    title = {{JuMP} 1.0: {R}ecent improvements to a modeling language for mathematical optimization},
    journal = {Mathematical Programming Computation},
    year = {2023},
    doi = {10.1007/s12532-023-00239-3}
}
```

### Graphs.jl

A powerful graph analysis tool in Julia

```bibtex
@misc{Graphs2021,
  author       = {Fairbanks, James and Besan{\c{c}}on, Mathieu and Simon, Sch{\"o}lly and Hoffiman, J{\'u}lio and Eubank, Nick and Karpinski, Stefan},
  title        = {JuliaGraphs/Graphs.jl: an optimized graphs package for the Julia programming language},
  year         = 2021,
  url = {https://github.com/JuliaGraphs/Graphs.jl/}
}
```
### LinearOperator

A Julia Linear Operator Package

```bibtex
@software{leconte2023linearoperators,
  title = {LinearOperators.jl: Linear Operators for Julia},
  author = {Leconte, Geoffroy and Orban, Dominique and Soares Siqueira, Abel and contributors},
  year = {2023},
  month = {12},
  day = {11},
  version = {2.6.0},
  doi = {10.5281/zenodo.2559295},
  url = {https://github.com/JuliaSmoothOptimizers/LinearOperators.jl},
  keywords = {Linear Algebra, Julia, Matrix Factorization},
  license = {MPL-2.0}
}
```

### Krylov

A Julia basket of hand-picked Krylov methods

```bibtex
@article{montoison2023krylov,
  title = {Krylov.jl: A Julia basket of hand-picked {Krylov} methods},
  author = {Montoison, Alexis and Orban, Dominique},
  journal = {Journal of Open Source Software},
  volume = {8},
  number = {89},
  pages = {5187},
  year = {2023},
  month = {9},
  day = {26},
  publisher = {Open Journals},
  doi = {10.21105/joss.05187},
  url = {https://joss.theoj.org/papers/10.21105/joss.05187},
  issn = {2475-9066}
}
```

### KrylovKit

A Julia package collecting a number of Krylov-based algorithms for linear problems, singular value and eigenvalue problems and the application of functions of linear maps or operators to vectors.

```bibtex
@software{haegeman2024krylovkit,
  author = {Haegeman, Jutho},
  title = {KrylovKit},
  version = {0.7.0},
  year = {2024},
  month = {3},
  day = {14},
  doi = {10.5281/zenodo.10622234},
  url = {https://github.com/Jutho/KrylovKit.jl}
}
```

## Acknowledgments

We would like to thank all the contributors to the Julia ecosystem whose work has made this package possible. Special thanks to the developers of the packages listed above for their valuable tools that form the foundation of HyDistFlow.jl.