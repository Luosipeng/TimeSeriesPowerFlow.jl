# References

## How to Cite TSPF.jl

If you use TSPF.jl in your research, please use the following citation format:

```bibtex
@misc{luo2025tspf,
  author = {Luo, Sipeng},
  title = {TSPF.jl: A Julia Package for Time Series Power Flow Analysis},
  year = {2025},
  url = {https://github.com/Luosipeng/TSPF.jl}
}
```

## External Packages

TSPF.jl relies on several external packages. We acknowledge and thank the contributors of these packages:

### AMD.jl

A Julia interface to the AMD library of Amestoy, Davis and Duff.

```bibtex
@misc{amd2023,
  author = {Montoison},
  title = {AMD.jl: A Julia interface to the AMD library of Amestoy, Davis and Duff},
  year = {2023},
  url = {https://github.com/JuliaSmoothOptimizers/AMD.jl}
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


## Note on Julia Standard Libraries

This package also uses the following standard libraries that are part of Julia:
- Dates: Julia's standard library for working with dates and time

## Acknowledgments

We would like to thank all the contributors to the Julia ecosystem whose work has made this package possible. Special thanks to the developers of the packages listed above for their valuable tools that form the foundation of TSPF.jl.