
function _contains_generator_emissions(data)
    contains_emissions = true
    for (id,gen) in data["gen"]
        if !haskey(gen, "emissions")
            println("Gen $id does not have key `emissions`")
            contains_emissions = false
        end
    end
    return contains_emissions
end

function _no_duplicate_costs_at_bus(data)
    no_duplicate = true
    bus_costs = Dict()
    for (id,gen) in data["gen"]
        if haskey(bus_costs, gen["gen_bus"])
            if gen["cost"] in bus_costs[gen["gen_bus"]]
                println("Gen $id has duplicate cost at bus $(gen["gen_bus"])")
                no_duplicate = false
                break
            else
                push!(bus_costs[gen["gen_bus"]], gen["cost"])
            end
        else
            bus_costs[gen["gen_bus"]] = [gen["cost"]]
        end
    end
    return no_duplicate
end

function _all_gen_costs_single_type(data)
    same_type = true
    first_type = collect(values(data["gen"]))[1]["model"]
    for (id,gen) in data["gen"]
        if !(gen["model"]==first_type)
            print("Generators have multiple cost function types, not supported.")
            same_type=false
        end
    end

    return same_type

end


function _build_dcopf(ref::Dict{Symbol,<:Any}, optimizer;cost_type="PWL")
    model = JuMP.Model()
    JuMP.set_optimizer(model, optimizer)

    # Add Optimization and State Variables
    # ------------------------------------
    JuMP.@variable(model, va[i in keys(ref[:bus])])
    JuMP.@variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])

    # Add Objective Function
    # ----------------------
    # JuMP.@objective(model, Min,
    #     sum(gen["cost"][2]*pg[i] for (i,gen) in ref[:gen])
    # )

    if cost_type=="PWL"
        #create dummy cost variable for each generator. Will be constrained such that cg is higher than all lines formed by extending PWL segments (assumes convexity)
        JuMP.@variable(model,cg[i in keys(ref[:gen])])
        #constrain cg variables
        for gen_idx in keys(ref[:gen])

            gen_cost = ref[:gen][gen_idx]["cost"]
            cost_points = []

            #create easier to deal with cost points array
            for i in 1:length(gen_cost)-0 #loop through gen["cost] array and grab points (-0 silences annoying julia warning)
                if i%2 != 0 #grab x points (odd indices)
                    push!(cost_points,[gen_cost[i],gen_cost[i+1]]) #add point to list
                end

            end

            #create constraints for the example generator
            for i in 1:length(cost_points)-1
                pt_1 = cost_points[i]
                pt_2 = cost_points[i+1]

                delta = pt_2.-pt_1
                slope = delta[2]/delta[1]

                JuMP.@constraint(model, cg[gen_idx] >= slope*(pg[gen_idx]-pt_1[1]) + pt_1[2]) #constrain cg>=m(pg-x_1) + y_1
            end

        end

        #and finally, set objective to be minimization of sum of cg variables
        JuMP.@objective(model, Min, sum(cg))


    elseif cost_type=="POLY"
        #no support for quadratic terms, use only linear and constant
        #this assumes that these are the final two terms in a polynomial cost fn
        JuMP.@objective(model, Min,
        sum(gen["cost"][end-1]*pg[i] + gen["cost"][end] for (i,gen) in ref[:gen])
        )

    else
        error("Unsupported cost function type $(cost_type) provided")
    end

    # Add Constraints
    # ---------------

    # Fix the voltage angle to zero at the reference bus
    for (i,bus) in ref[:ref_buses]
        JuMP.@constraint(model, va[i] == 0)
    end

    for (i,branch) in ref[:branch]
        # Compute the branch parameters and transformer ratios from the data
        g, b = PowerModels.calc_branch_y(branch)
        branch["g"] = g
        branch["b"] = b

        k = branch["f_bus"]
        j = branch["t_bus"]

        #constrain the thermal limit
        JuMP.@constraint(model, -branch["b"]*(va[k] - va[j]) <= branch["rate_a"])
        JuMP.@constraint(model, -branch["b"]*(va[k] - va[j]) >= -branch["rate_a"])

    end

    # Nodal power balance constraints
    for i in sort(collect(keys(ref[:bus])))
        ref[:bus][i]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        # Flows out of bus i are equal to generation at bus i minus demand at bus i
        JuMP.@constraint(model,
            sum(-ref[:branch][l]["b"]*(va[k] - va[j]) for (l,k,j) in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -                 # sum of active power generation at bus i -
            sum(load["pd"] for load in bus_loads) -                 # sum of active load consumption at bus i -
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2          # sum of active shunt element injections at bus i
        )
    end
    return model
end


function _solve_dcopf(data, optimizer;cost_type="PWL")
    model = _build_dcopf(data, optimizer,cost_type=cost_type)
    JuMP.optimize!(model)
    return model
end


#this function updates ref with the solved line flows (not contained in JuMP model but calculated from voltage angles)
function _update_solved_line_flow!(ref,model::JuMP.AbstractModel)
    for branch in values(ref[:branch])
        k = branch["f_bus"]
        j = branch["t_bus"]
        b = branch["b"]

        line_flow = -b*(JuMP.value(model[:va][k]) - JuMP.value(model[:va][j]))

        branch["pf"] = line_flow
        branch["pt"] = -line_flow
    end

    return
end

#this function does the same thing for voltage angles on all buses
function _update_solved_va!(ref,model::JuMP.AbstractModel)
    for bus in values(ref[:bus])
        bus_id = bus["index"]
        bus["va"] = JuMP.value(model[:va][bus_id])
    end
end


#and this one does the same for all generator set points
function _update_solved_pg!(ref,model::JuMP.AbstractModel)
    for gen in values(ref[:gen])
        gen_id = gen["index"]
        gen["pg"] = JuMP.value(model[:pg][gen_id])

    end

end

function update_ref!(ref,model::JuMP.AbstractModel)
    _update_solved_line_flow!(ref,model)
    _update_solved_va!(ref,model)
    _update_solved_pg!(ref,model)
end




function collect_binding_constraints(model::JuMP.AbstractModel, N::Int)
    eq_cons = collect_equality_constraints(model, N)
    binding_ineq_cons = collect_binding_inequality_constraints(model)

    binding_cons = [eq_cons;binding_ineq_cons]
    return binding_cons
end


function collect_equality_constraints(model::JuMP.AbstractModel, N::Int)
    equality_constraints_ref = JuMP.all_constraints(model, JuMP.AffExpr,_MOI.EqualTo{Float64})
    eq_cons = [JuMP.constraint_object(x).func for x in equality_constraints_ref[2:N+1]]     # select equality constraints EXCEPT ref_bus angle constraint
    ref_bus_eq_con = JuMP.constraint_object(equality_constraints_ref[1]).func     #select ref_bus angle=0 constraint
    eq_cons = [eq_cons;ref_bus_eq_con].*-1     #stack and flip sign
    return eq_cons
end


function collect_binding_inequality_constraints(model::JuMP.AbstractModel)
    var_gt = JuMP.all_constraints(model, JuMP.VariableRef,_MOI.GreaterThan{Float64})
    var_lt = JuMP.all_constraints(model, JuMP.VariableRef,_MOI.LessThan{Float64})
    con_gt = JuMP.all_constraints(model, JuMP.AffExpr,_MOI.GreaterThan{Float64})
    con_lt = JuMP.all_constraints(model, JuMP.AffExpr,_MOI.LessThan{Float64})

    all_ineq = [var_gt;var_lt;con_gt;con_lt]

    binding_ineq_refs = []
    for ineq in all_ineq
        if abs(JuMP.dual(ineq))>0
            push!(binding_ineq_refs,ineq)
        end
    end

    return [JuMP.constraint_object(con).func for con in binding_ineq_refs] #convert to function representation (from object)
end


function assemble_A_poly(model::JuMP.AbstractModel)
    Nb = length(model[:va])
    Ng = length(model[:pg])
    n = Nb + Ng #voltage angle and pg variables
    A = zeros(n,n)

    all_binding_constraints = collect_binding_constraints(model, Nb) #Nb passed in in order to identify the nodal power balance constraints

    #println("There are $(length(all_binding_constraints)) total binding constraints")

    i = 1
    for con in all_binding_constraints
        #initialize row
        temp_row = zeros(n)

        #grab the voltage angle and pg vars in order
        sorted_va_idx = sort(axes(model[:va])[1])
        sorted_pg_idx = sort(axes(model[:pg])[1])

        va_vars = Array(model[:va][sorted_va_idx])
        pg_vars = Array(model[:pg][sorted_pg_idx])

        #and combine them into one list
        va_pg_vars = [va_vars;pg_vars]

        if typeof(con)==JuMP.VariableRef #if constraint is just variable ref, then row in A should just have a 1 where that variable is
            var_idx = findfirst(==(con),va_pg_vars)
            temp_row[var_idx]=1
        else #otherwise if constraint is proper funtion (equality constraint or non-variable-definition inequality)
            #now grab the constraint function that we're encoding
            con_func_terms = con.terms
            #and loop through the variables in va_pg_vars, checking whether they're in the con_func (and if so, grabbing their coefficient)
            for (j,var) in enumerate(va_pg_vars)
                if var in keys(con_func_terms)
                    temp_row[j] = con_func_terms[var]
                end
            end
        end

        #make temp row permanent
        A[i,:] = temp_row
        i = i+1
    end
    return A
end


function assemble_A_pwl(model::JuMP.AbstractModel)
    Nb = length(model[:va])
    Ng = length(model[:pg])
    n = Nb + Ng + Ng #2x Ng because one cg variable for each gen as well as one pg variable
    A = zeros(n,n)

    all_binding_constraints = collect_binding_constraints(model, Nb) #Nb passed in in order to identify the nodal power balance constraints

    #println("There are $(length(all_binding_constraints)) total binding constraints")

    i = 1
    for con in all_binding_constraints
        #initialize row
        temp_row = zeros(n)

        #grab the voltage angle and pg vars in order
        sorted_va_idx = sort(axes(model[:va])[1])
        sorted_pg_idx = sort(axes(model[:pg])[1])
        sorted_cg_idx = sort(axes(model[:cg])[1])

        va_vars = Array(model[:va][sorted_va_idx])
        pg_vars = Array(model[:pg][sorted_pg_idx])
        cg_vars = Array(model[:cg][sorted_cg_idx])

        #and combine them into one list
        va_pg_cg_vars = [va_vars;pg_vars;cg_vars]

        if typeof(con)==JuMP.VariableRef #if constraint is just variable ref, then row in A should just have a 1 where that variable is
            var_idx = findfirst(==(con),va_pg_cg_vars)
            temp_row[var_idx]=1
        else #otherwise if constraint is proper funtion (equality constraint or non-variable-definition inequality)
            #now grab the constraint function that we're encoding
            con_func_terms = con.terms
            #and loop through the variables in va_pg_vars, checking whether they're in the con_func (and if so, grabbing their coefficient)
            for (j,var) in enumerate(va_pg_cg_vars)
                if var in keys(con_func_terms)
                    temp_row[j] = con_func_terms[var]
                end
            end
        end

        #make temp row permanent
        A[i,:] = temp_row
        i = i+1
    end
    return A
end



function get_total_bus_load(bus::Int,ref)
    total_load = 0
    bus_loads = ref[:bus_loads][bus] #grab list of loads at bus
    #and total them up
    for l in bus_loads
        total_load += ref[:load][l]["pd"]
    end

    return total_load
end

function get_total_outflow(node::Node)
    total_outflow = 0
    for outflow in node.outflows
        total_outflow+=outflow.total_flow
    end
    return total_outflow
end