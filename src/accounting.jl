# basic accounting

# Update pomwermodels data dicts

function update_emissions_intensity!(data::Dict{String,<:Any}, lmce::LMCE)
    # check that data bus and lmce bus are the same
    @assert _contains_matching_bus_ids(data, lmce)

    for (id,bus) in data["bus"]
        bus["emissions_intensity"] = lmce.emissions_intensity[findfirst(==(parse(Int,id)),lmce.node_ids)]
    end
end

function update_emissions_intensity!(data::Dict{String,<:Any}, almce::ALMCE)
    # check that data bus and lmce bus are the same
    @assert _contains_matching_bus_ids(data, almce)

    for (id,bus) in data["bus"]
        bus["emissions_intensity"] = almce.emissions_intensity[findfirst(==(parse(Int,id)),almce.node_ids)]
    end
end

function update_emissions_intensity!(data::Dict{String,<:Any}, lace::LACE)
    # check that data bus and lmce bus are the same
    @assert _contains_matching_bus_ids(data, lace)

    for (id,bus) in data["bus"]
        bus["emissions_intensity"] = lace.emissions_intensity[findfirst(==(parse(Int,id)),lace.node_ids)]
    end
end

function update_emissions_intensity!(data::Dict{String,<:Any}, ace::ACE)

    # check that data bus and lmce bus are the same
    for (id,bus) in data["bus"]
        bus["emissions_intensity"] = ace.emissions_intensity[findfirst(==(parse(Int,id)),ace.node_ids)]
    end
end



function _contains_matching_bus_ids(nodal_intensity::AbstractNodalCarbonIntensity, data::Dict{String,<:Any})
    _contains_matching_bus_ids(data, nodal_intensity)
end

function _contains_matching_bus_ids(data::Dict{String,<:Any}, nodal_intensity::AbstractNodalCarbonIntensity)
    assertion = true

    data_bus_ids = parse.(Int,keys(data["bus"]))
    for id in data_bus_ids
        if !in(id, nodal_intensity.node_ids)
            println("Bus $id is not in the carbon intensity metric")
            assertion = false
        end
    end

    for id in nodal_intensity.node_ids
        if !in(id, data_bus_ids)
            println("Bus $id is not in power grid data")
            assertion =  false
        end
    end

    return assertion
end

function _contains_emissions_intensity(data::Dict{String,<:Any})
    for (id,bus) in data["bus"]
        if !haskey(bus, "emissions_intensity")
            println("Bus $id does not have key `emissions_intensity`")
            return false
        end
    end
    return true
end
