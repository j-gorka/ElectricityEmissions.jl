module ElectricityEmissionsExt

using ElectricityEmissions
using PowerPlots
# using ColorSchemes
# using Colors
using Setfield

const co2_intensity_scheme = ["#33a663", "#dcde4f","#c28838","#672e14","#2a1602"]

function ElectricityEmissions.plot_emissions(data::Dict{String,<:Any};
    components = ["bus","branch","gen"],
    bus_data = :emissions_intensity,
    bus_color = co2_intensity_scheme,
    bus_data_type = :quantitative,
    gen_data = :emissions,
    gen_color = co2_intensity_scheme,
    gen_data_type = :quantitative,
    branch_color = ["#C6DBEF", "#d73027"],
    gen_size_data = :pg,
    gen_size_range = [50,500],
    bus_size_data = :pd,
    bus_size_range = [50,500],
    show_branch_utilization=false,
    title="Emissions Intensity",
    power_domain = nothing,
    intensity_domain = nothing,
    kwargs...)

    # assert contains intensity data
    @assert ElectricityEmissions._contains_emissions_intensity(data)

    data_copy = deepcopy(data)
    for (id,bus) in data_copy["bus"]
        bus["pd"] = 0.0
    end
    for (id,load) in data_copy["load"]
        data_copy["bus"][string(load["load_bus"])]["pd"] += load["pd"]
    end

    # size domain
    if isnothing(power_domain)
        min_power, max_power = extrema_power_values(data_copy)
        power_domain = [min_power, max_power]
    end

    # color domain
    if isnothing(intensity_domain)
        min_intensity, max_intensity = extrema_emissions_intensity_values(data_copy)
        intensity_domain = [min_intensity, max_intensity]
    end

    # Create base plot
    p=powerplot(data_copy,
        components=components,
        bus_data=bus_data,
        bus_color=bus_color,
        bus_data_type=bus_data_type,
        gen_data=gen_data,
        gen_color=gen_color,
        gen_data_type=gen_data_type,
        branch_color=branch_color,
        ;
        kwargs...
    )

    # Adjust color
    #one is gen, one is bus
    p.layer[3]["encoding"]["color"]["title"] = "Emissions Intensity"
    p.layer[3]["encoding"]["color"]["scale"]["domain"] = intensity_domain
    p.layer[3]["encoding"]["color"]["legend"] = Dict("gradientLength"=>150, "orient"=>"bottom")
    p.layer[4]["encoding"]["color"]["title"] = "Emissions Intensity"
    p.layer[4]["encoding"]["color"]["scale"]["domain"] = intensity_domain
    p.layer[4]["encoding"]["color"]["legend"] = Dict("gradientLength"=>150, "orient"=>"bottom")

    # size is load/pg
    # p.layer[3]["encoding"]["size"]=Dict(
    #     "field"=>bus_size_data, "title"=>"Power",
    #     "type"=>"quantitative",
    #     "scale"=>Dict{String,Any}("range"=>bus_size_range, "domain"=>power_domain),
    #     "legend"=>Dict("orient"=>"right")
    # )
    # p.layer[4]["encoding"]["size"]=Dict(
    #     "field"=>gen_size_data, "title"=>"Power",
    #     "type"=>"quantitative",
    #     "scale"=>Dict{String,Any}("range"=>gen_size_range, "domain"=>power_domain),
    #     "legend"=>Dict("orient"=>"right")
    # )


    # set branch color values
    if show_branch_utilization
        p.layer[1]["transform"] = Dict{String, Any}[
            Dict("calculate"=>"abs(datum.pt)/datum.rate_a*100", "as"=>"branch_Percent_Loading"),
            Dict("calculate"=>"abs(datum.pt)", "as"=>"BranchPower")
        ]
        p.layer[1]["layer"][1]["encoding"]["color"]["field"]="branch_Percent_Loading"
        p.layer[1]["layer"][1]["encoding"]["color"]["type"]="quantitative"
        p.layer[1]["layer"][1]["encoding"]["color"]["title"]="Branch Utilization %"
        p.layer[1]["layer"][1]["encoding"]["color"]["scale"]["domain"]=[0,100]
        p.layer[1]["layer"][1]["encoding"]["color"]["legend"] = Dict("gradientLength"=>150, "orient"=>"bottom")

    end


    # Set gen shape
    p.layer[4]["mark"]["type"]=:square

    # set legend
    @set! p.resolve.scale.size=:independent
    @set! p.resolve.scale.color=:independent

    p.layer[4]["encoding"]["size"]["legend"] = Dict("disable"=>true)
    p.layer[4]["encoding"]["color"]["legend"]["disable"]=true
    p.layer[2]["encoding"]["color"]["legend"] = Dict("disable"=>true)

    # set title
    @set! p.title=title

    return p
end


function extrema_power_values(data::Dict{String,<:Any})
    min_load, max_load = extrema(bus["pd"] for (id,bus) in data["bus"])
    min_pg, max_pg = extrema(gen["pg"] for (id,gen) in data["gen"])
    min_power, max_power = min(min_load, min_pg), max(max_load, max_pg)
    return min_power, max_power
end

function extrema_emissions_intensity_values(data::Dict{String,<:Any})
    min_node_intensity, max_node_intensity = extrema(bus["emissions_intensity"] for (id,bus) in data["bus"])
    min_gen_intensity, max_gen_intensity = extrema(gen["emissions"] for (id,gen) in data["gen"])
    min_intensity, max_intensity = min(min_node_intensity, min_gen_intensity), max(max_node_intensity, max_gen_intensity)
    return min_intensity, max_intensity
end


end