module ElectricityEmissions

    import PowerModels
    import JuMP
    import MathOptInterface

    _PM = PowerModels
    _MOI = MathOptInterface

    include("types.jl")

    include("accounting.jl")
    include("calculate_emissions.jl")
    include("utils.jl")
    include("plot_prep.jl")

export LMCE, calculate_LMCE, ALMCE, calculate_ALMCE, ACE, calculate_ACE, LACE, calculate_LACE, plot_emissions,update_emissions_intensity!
end
