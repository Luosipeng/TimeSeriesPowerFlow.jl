# TSPflow.jl

[![Documentation](https://github.com/Luosipeng/TSPflow.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/Luosipeng/TSPflow.jl/actions/workflows/documentation.yml)
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://luosipeng.github.io/TSPflow.jl/stable)
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://luosipeng.github.io/TSPflow.jl/dev)

## Description

TSPflow.jl is a Julia package developed by the HR-PES team of Xi'an Jiaotong University, which provides a dynamic power flow simulation framework for distribution systems incorporating various renewable resources.

## Installation

You can install the package via the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add TSPflow
```

Or, you can install it directly in Julia code:

```julia
using Pkg
Pkg.add("TSPflow")
```

## Documentation

Comprehensive documentation is available online:

- **Stable version**: [https://luosipeng.github.io/TSPflow.jl/stable](https://luosipeng.github.io/TSPflow.jl/stable)
- **Development version**: [https://luosipeng.github.io/TSPflow.jl/dev](https://luosipeng.github.io/TSPflow.jl/dev)


## Features

- **Integration of PF and OPF**: Using relaxed OPF to bridge time-varying loads and generation profiles with ESS status, integrating VSC power allocation into Power Flow calculations
- **Customizable Simulation Environment**: Comprehensive parameter customization including irradiance data, electricity prices, load profiles, and flexible network topology configuration
- **Comprehensive Renewable Energy Models**: Validated mathematical models for energy storage, photovoltaic generation, and VSCs for reliable distribution network analysis
- **Advanced VSC Control Framework**: Implementation of seven control modes (voltage, reactive power, droop control, etc.) with corresponding solution algorithms
- **Flexible Case Management**: Support for both externally imported and internally generated case studies with complete data structures

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use TSPflow.jl in your research, please cite:

```bibtex
@software{TSPflow2025,
  author = {Sipeng Luo,Tianyang Zhao,Zhaohong Bie},
  title = {TSPflow.jl: a Julia package for time series power flow analysis},
  year = {2025},
  url = {https://github.com/Luosipeng/TSPflow.jl}
}
```

## Contact

For questions and feedback, please open an issue on the [GitHub repository](https://github.com/Luosipeng/TSPflow.jl/issues) or contact the maintainer.