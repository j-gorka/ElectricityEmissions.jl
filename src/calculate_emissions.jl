

function calculate_ACE(file::String, optimizer)
    data = PowerModels.parse_file(file)
    calculate_ACE(data, optimizer)
end


# function calc_LMCE(data;return_delta_pg=false,adjusted_LMCE = false,account=false)
function calculate_ACE(data::Dict{String,<:Any}, optimizer)
    # data checks
    @assert _contains_generator_emissions(data)

    #PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = _solve_dcopf(ref, optimizer)

    total_system_emissions = sum(JuMP.value(model[:pg][id])*gen["emissions"] for (id,gen) in ref[:gen]) #get sum of all system emissions
    avg_system_emissions = total_system_emissions/sum(JuMP.value(model[:pg][id]) for (id,gen) in ref[:gen]) #get average system emissions

    #grab bus bus_ids and create vector of ACE emissions (same for each bus) to match format of other emissions metrics
    bus_ids = sort(collect(keys(ref[:bus])))
    bus_ace = ones(length(bus_ids)).*avg_system_emissions

    return ACE(bus_ids,bus_ace)
end


function calculate_system_emissions(data::Dict{String,<:Any}, optimizer)
    # data checks
    @assert _contains_generator_emissions(data)

    #PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = _solve_dcopf(ref, optimizer)

    total_system_emissions = sum(JuMP.value(model[:pg][id])*gen["emissions"] for (id,gen) in ref[:gen]) #get sum of all system emissions

    return total_system_emissions

end




#recursive function for graph traversal portion of LACE calculation, modifies input data dictionary
function _recurse_LACE!(curr_node::Node,gen::Int,system_graph::Graph, ref, parent_line::Union{Line,Nothing})
    #current node calculations
    if isnothing(parent_line) #equivalent to checking whether root node
        #in-flow contributed by gen is known (just gen set point), no need to store anything
        in_flow = ref[:gen][gen]["pg"]
    else
        #1. grab gen contribution of parent in-flow
        #in_flow = parent_line.gen_cont[gen]
        in_flow = parent_line.gen_cont[gen][end] #get the current path value for the gen cont on the inflow line (set on prev level outflow calc)
    end

    #2. calculate gen contribution on loads, store this
    total_load = get_total_bus_load(curr_node.bus_id,ref)
    total_outflow = get_total_outflow(curr_node)
    load_gen_cont = in_flow * (total_load/(total_load + total_outflow)) #directly from paper


    for outflow_line in curr_node.outflows #these are Line objects
        #calculate gen contribution on particular outflow that you're going to recurse on...
        #...store this by incrementing the gen contribution of the line in question

        #calculate gen contribution to outflow, based on eq in paper
        outflow_cont = in_flow * (outflow_line.total_flow/(total_load + total_outflow))
        #then add outflow gen contribution (along current path) to outflow Line
        #outflow_line.gen_cont[gen]+=outflow_cont
        append!(outflow_line.gen_cont[gen], outflow_cont)


        outflow_node = system_graph.nodes[outflow_line.to_node]
        #note that outflow_node is the node we're going to, and outflow_line will be it's 'parent line' at the next level of recursion
        _recurse_LACE!(outflow_node,gen,system_graph,ref,outflow_line)
    end


end


function calculate_LACE(file::String, optimizer)
    data = PowerModels.parse_file(file)
    calculate_LACE(data, optimizer)
end


function calculate_LACE(data::Dict{String,<:Any}, optimizer)
    @assert _contains_generator_emissions(data)
    @assert _no_duplicate_costs_at_bus(data)

    #PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = _solve_dcopf(ref, optimizer)
    update_ref!(ref,model)


    #grab number of generators
    num_gens = length(ref[:gen])

    #instantiate graph
    system_graph = Graph()

    # will loop through lines, adding them as objects as well as their endpoint Node objects
    # both will be added to the overall Graph object
    for line in values(ref[:branch])
        #get f and t bus based on power flow direction, flip if pf <0
        if line["pf"]>=0
            line_flow = line["pf"]
            f_bus = line["f_bus"]
            t_bus = line["t_bus"]
        else
            line_flow = line["pf"]*-1
            t_bus = line["f_bus"]
            f_bus = line["t_bus"]
        end

        #create Node objects for endpoints of line
        f_node = Node(f_bus,zeros(num_gens))


        if !(f_bus in collect(keys(system_graph.nodes)))
            f_node = Node(f_bus,zeros(num_gens))
            system_graph.nodes[f_bus] = f_node
        else
            f_node = system_graph.nodes[f_bus]
        end

        if !(t_bus in collect(keys(system_graph.nodes)))
            t_node = Node(t_bus,zeros(num_gens))
            system_graph.nodes[t_bus] = t_node
        else
            t_node = system_graph.nodes[t_bus]
        end


        #and create the Line object itself, using total flow stored before. Gen contributions are left until later, as with Nodes
        #line_obj = Line(f_bus,t_bus,line_flow,zeros(num_gens))
        line_obj = Line(f_bus,t_bus,line_flow,[Float64[] for i=1:num_gens])


        #and finally add this line object to the list of outflows for the f_node, and the list of in-flows for the t_node
        #as well as the master line list stored by the system graph
        push!(f_node.outflows,line_obj)
        push!(t_node.inflows,line_obj)
        push!(system_graph.lines,line_obj)

    end

    #calculate generator contributions across all nodes, all generators

    for gen in sort(collect(keys(ref[:gen])))
        gen_bus = ref[:gen][gen]["gen_bus"]
        gen_node = system_graph.nodes[gen_bus]
        _recurse_LACE!(gen_node, gen, system_graph,ref, nothing)

        for node in collect(values(system_graph.nodes))
            total_load = get_total_bus_load(node.bus_id,ref)
            total_outflow = get_total_outflow(node)
            total_load_gen_cont = 0 #initialize gen cont to zero
            for inflow_line in node.inflows #go through all in-flows, calculating load cont related to each inflow
                #note below now sums across line gen cont across all paths for calculation:
                inflow_load_gen_cont = sum(inflow_line.gen_cont[gen]) * (total_load/(total_load + total_outflow))
                total_load_gen_cont += inflow_load_gen_cont #increment total
            end
            #also calculate contribution generator local to current node
            if node.bus_id == ref[:gen][gen]["gen_bus"]
                local_gen_cont = ref[:gen][gen]["pg"] * (total_load/(total_load + total_outflow))
                total_load_gen_cont += local_gen_cont
            end

            node.gen_cont[gen] = total_load_gen_cont #assign to node
        end
    end

    #finally, loop through all nodes to calculate LACE at each (using gen contributions and gen emissions intensities)
    bus_ids = sort(collect(keys(ref[:bus])))
    bus_lace = []
    bus_gen_cont = []
    sorted_gen_emission_factors = [ref[:gen][id]["emissions"] for id in sort(collect(keys(ref[:gen])))]

    for bus in bus_ids
        bus_load = get_total_bus_load(bus,ref)
        node = system_graph.nodes[bus] #grab node object
        if bus_load >0
            node_LACE = sum(node.gen_cont.*sorted_gen_emission_factors)/bus_load #inner product of gen intensities and gen contributions, divided by nodal load to get carbon intensity
        else
            node_LACE = 0 #set intensity to zero in case of zero-load bus
        end
        push!(bus_lace,node_LACE)
        push!(bus_gen_cont, node.gen_cont)

    end

    # for node in collect(values(system_graph.nodes))
    #     node_LACE = sum(node.gen_cont.*sorted_gen_emission_factors) #inner product of gen intensities, gen contributions
    #     push!(bus_ids,node.bus_id)
    #     push!(bus_lace,node_LACE)

    # end

    lace = LACE(bus_ids,bus_lace,bus_gen_cont)

    return lace


end



function calculate_LMCE(file::String, optimizer)
    data = PowerModels.parse_file(file)
    calculate_LMCE(data, optimizer)
end


# function calc_LMCE(data;return_delta_pg=false,adjusted_LMCE = false,account=false)
function calculate_LMCE(data::Dict{String,<:Any}, optimizer)
    # data checks
    @assert _contains_generator_emissions(data)
    @assert _no_duplicate_costs_at_bus(data)

    #PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = _solve_dcopf(ref, optimizer)

    ## Calculate ΔPg ##
    #initialize A matrix to be NxN where n = N_b + N_g (number of buses plus number of generators)
    A = assemble_A(model)
    A_inv = inv(A)

    #grab relevant portion of A_inv matrix (Last num_gens rows, first num_bus columns of A_inv)
    Nb = length(model[:va])
    Ng = length(model[:pg])
    B = A_inv[Nb+1:Ng+Nb,1:Nb]

    ## Calculate lmce by multiplying B (generator response) by generator emissions intensity (emissions/MWh)
    sorted_gen_emission_factors = [ref[:gen][id]["emissions"] for id in sort(collect(keys(ref[:gen])))]
    for r in eachcol(B)
        r .=  r .* sorted_gen_emission_factors
    end

    bus_ids = sort(collect(keys(ref[:bus])))
    bus_lmce = vec(sum(B,dims=1))
    return LMCE(bus_ids, bus_lmce)
end

function calculate_ALMCE(file::String, optimizer)
    data = PowerModels.parse_file(file)
    calculate_ALMCE(data, optimizer)
end


# function calc_LMCE(data;return_delta_pg=false,adjusted_LMCE = false,account=false)
function calculate_ALMCE(data::Dict{String,<:Any}, optimizer)
    # data checks
    @assert _contains_generator_emissions(data)
    @assert _no_duplicate_costs_at_bus(data)

    #PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    model = _solve_dcopf(ref, optimizer)

    ## Calculate ΔPg ##
    #initialize A matrix to be NxN where n = N_b + N_g (number of buses plus number of generators)
    A = assemble_A(model)
    A_inv = inv(A)

    #grab relevant portion of A_inv matrix (Last num_gens rows, first num_bus columns of A_inv)
    Nb = length(model[:va])
    Ng = length(model[:pg])
    B = B = A_inv[Nb+1:Ng+Nb,1:Nb]

    ## Calculate emissions/MWh
    sorted_gen_emission_factors = [ref[:gen][id]["emissions"] for id in sort(collect(keys(ref[:gen])))]
    for r in eachcol(B)
        r .=  r .* sorted_gen_emission_factors
    end

    bus_ids = sort(collect(keys(ref[:bus])))
    bus_lmce = vec(sum(B,dims=1))
    lmce = LMCE(bus_ids, bus_lmce)

    # adusted LMCE
    sorted_bus_list = sort(collect(keys(ref[:bus_loads]))) #grab list of buses
    sorted_bus_loads = zeros(length(sorted_bus_list)) #initialize zeros for total load at each bus

    #loop through buses and their constituent loads to get total system load
    for (bus_idx,bus) in enumerate(sorted_bus_list)
        bus_loads = ref[:bus_loads][bus] #grab list of loads at bus
        for load in bus_loads
            sorted_bus_loads[bus_idx] += ref[:load][load]["pd"] #add real power demand for load to corresponding bus
        end
    end

    #grab lmce-assigned emissions and total system load:
    lmce_bus_emissions = sorted_bus_loads.*lmce.emissions_intensity
    total_system_load = sum(sorted_bus_loads)

    #now need to get current generator set points and calculate total system emissions
    sorted_pg_idx = sort(axes(model[:pg])[1])
    sorted_pg_vars = Array(model[:pg][sorted_pg_idx])
    gen_set_points = JuMP.value.(sorted_pg_vars)
    total_system_emissions = sum(gen_set_points.*sorted_gen_emission_factors) #get sum of all system emissions

    #calculate emissions discrepancy with LMCE-assigned emissions
    emissions_discrep = total_system_emissions - sum(lmce_bus_emissions)

    LMCE_adjustment_factor = emissions_discrep/total_system_load

    adjusted_LMCE = lmce.emissions_intensity .+LMCE_adjustment_factor
    return ALMCE(sorted_bus_list, adjusted_LMCE)
end