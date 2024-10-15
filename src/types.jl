
abstract type AbstractCarbonIntensity end

abstract type AbstractNodalCarbonIntensity <: AbstractCarbonIntensity end

abstract type AbstractSystemCarbonIntensity <: AbstractCarbonIntensity end

abstract type AbstractCarbonAccount end

mutable struct LMCE <: AbstractNodalCarbonIntensity
    node_ids::Vector{Int}
    emissions_intensity::Vector{Float32}
end


mutable struct ALMCE <: AbstractNodalCarbonIntensity
    node_ids::Vector{Int}
    emissions_intensity::Vector{Float32}
end


mutable struct ACE <: AbstractSystemCarbonIntensity
    node_ids::Vector{Int}
    emissions_intensity::Vector{Float32}
end



mutable struct Node
    bus_id::Int
    inflows::Vector{Any} #should be type Line, just can't define yet
    outflows::Vector{Any} #should be type Line...
    gen_cont::Union{Vector{Float64},Nothing}
end
#outer constructor to make in-flows and out-flows empty lists to start, gen_cont allowed to be nothing. Will sort that out later
Node(bus_id::Int) = Node(bus_id::Int,Any[],Any[],nothing)
Node(bus_id::Int,gen_cont::Vector{Float64}) = Node(bus_id::Int,Any[],Any[],gen_cont::Vector{Float64})


mutable struct Line
    from_node::Int
    to_node::Int
    total_flow::Float64
    gen_cont::Union{Vector{Vector{Float64}},Nothing}
end
#outer constructor for default values (might be needed?)
Line(from_node::Int,to_node::Int,total_flow::Float64) = Line(from_node::Int,to_node::Int,total_flow::Float64,nothing)


#graph is made up of Nodes and Lines, useful for storing Node/Line objects (or modifying them)
mutable struct Graph
    lines::Vector{Line}
    nodes::Dict{Int,Node}
end

#outer constructor for default values:
Graph() = Graph(Line[],Dict{Int,Node}())


mutable struct LACE <: AbstractNodalCarbonIntensity
    node_ids::Vector{Int}
    emissions_intensity::Vector{Float32}
    gen_cont::Union{Vector{Vector{Float64}},Nothing}
end
LACE(node_ids::Vector{Int}, emissions_intensity::Vector{Float32}) = LACE(node_ids::Vector{Int}, emissions_intensity::Vector{Float32},nothing)