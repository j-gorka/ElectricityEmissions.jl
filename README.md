# ElectricityEmissions.jl

[![Build Status](https://github.com/j-gorka/ElectricityEmissions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/j-gorka/ElectricityEmissions.jl/actions/workflows/CI.yml?query=branch%3Amain)



ElectricityEmissions.jl is a Julia package for the calculation of carbon emissions intensity metrics for power system test cases. It also includes visualization functionality via the [PowerPlots.jl](https://github.com/WISPO-POP/PowerPlots.jl) package.

## Installation
The latest development version of `ElectricityEmissions.jl` can be added via:
```Pkg> add https://github.com/WISPO-POP/ElectricityEmissions.jl.git```


## Basic Use
To use `ElectricityEmissions.jl`, simply load a test case in MATPOWER format, assign emissions intensity values to all generators, and calculate nodal carbon intensity values. Optionally, it is then possible to visualize these nodal intensity values.

```
using ElectricityEmissions, PowerModels, HiGHS

case = PowerModels.parse_file("./test_cases/case30pwl.m")

#randomly assign generator emissions factors between 0 and 1
for (id,gen) in case["gen"]
	gen["emissions"] = rand()
end


#calculate LMCE
lmce = calculate_LMCE(case,HiGHS.Optimizer)

#(Optional) Update test case with nodal emissions info and plot
update_emissions_intensity!(case,lmce)
plot_emissions(case)
```

`ElectricityEmissions.jl` currently supports the following carbon intensity calculation methods:
- `calculate_LMCE`: Locational Marginal Carbon Emissions
- `calculate_ACE`: Average Carbon Emissions
- `calculate_LACE`: Locational Average Carbon Emissions (as proposed in [Chen et. al.](https://arxiv.org/abs/2311.03712))
- `calculate ALMCE:` Adjusted Locational Marginal Carbon Emissions (as we propose in [ElectricityEmissions.jl: A Framework for the Comparison of Carbon Intensity Signals](https://arxiv.org/abs/2411.06560)


## Generator Cost Function Support
Currently, `ElectricityEmissions.jl` supports the use of test cases with linear or piece-wise linear generator cost functions. At present, quadratic and higher terms within polynomial cost functions will be ignored. However, An upcoming update will a add support for quadratic terms.

## Troubleshooting
The calculation of marginal metrics, such as LMCE and ALMCE, relies on the inversion of a constraint matrix and hence requires a unique solution to the generation dispatch problem (OPF). Situations in which there is more than one generator with the same cost function at a particular bus may result in a `SingularException`. To remedy this, we have included the function `add_gen_cost_noise!(case)`which adds a small amount of noise to all generator cost functions.

## Acknowledgements
This code is primarily developed by Joe Gorka \<jgorka@wisc.edu\>, Noah Rhodes \<nrhodes@lanl.gov\>, and Line Roald \<roald@wisc.edu\>.


## Citation Information
If you use ElectricityEmissions.jl in your project, we would appreciate if you would cite the following publication:

```
@misc{gorka2024electricityemissionsjl,
      title={ElectricityEmissions.jl: A Framework for the Comparison of Carbon Intensity Signals}, 
      author={Joe Gorka and Noah Rhodes and Line Roald},
      year={2024},
      eprint={2411.06560},
      archivePrefix={arXiv},
      primaryClass={eess.SY},
      url={https://arxiv.org/abs/2411.06560}, 
}
```


