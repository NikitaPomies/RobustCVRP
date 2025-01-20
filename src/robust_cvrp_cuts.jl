using JuMP, CPLEX, LinearAlgebra


include("instance.jl")

instance = read_instance("../data/n_10-euclidean_true")

function find_one_route_clients(pairs::Vector{Tuple{Int,Int}}, n::Int64)
    point_to_depot = Dict()
    depot_point_to = Dict()
    for pair in pairs
        if pair[1] == 1
            depot_point_to[pair[2]] = true
        end
        if pair[2] == 1
            point_to_depot[pair[1]] = true
        end
    end
    return [i for i in 1:n if get(point_to_depot, i, false) && get(depot_point_to, i, false)]
end


function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

function fractional_knapsack(weights::Dict{Tuple{Int, Int}, Float64}, 
                             values::Dict{Tuple{Int, Int}, Float64}, 
                             capacity::Float64, max_var::Float64 = 1.0)
    """
    Résoud le problème du sac à dos continu avec max_var la borne supérieure de la variable de choix continue
    """
    # Calcul des ratios valeur/poids pour chaque (i, j)
    items = [(values[k] / weights[k], k, weights[k], values[k]) for k in keys(weights)]
    
    # Tri par ratio décroissant
    sort!(items, by = x -> x[1], rev = true)
    
    total_value = 0.0
    remaining_capacity = capacity
    solution = Dict((i, j) => 0.0 for (i,j) in keys(weights))  

    for (_, (i, j), weight, value) in items
        if remaining_capacity == 0
            break
        end

        if weight*max_var <= remaining_capacity
            total_value += value*max_var
            remaining_capacity -= weight*max_var
            solution[(i, j)] = max_var  
        else
            fraction = remaining_capacity/weight
            total_value += value * fraction
            solution[(i, j)] = fraction 
            remaining_capacity = 0
        end
    end

    return total_value, solution
end


function solve_subproblem(I::Instance, x_bar)
    n = length(I.demands)

    selected_e = selected_edges(x_bar, n)
    
    values_1 = Dict((i, j) => 0. for i in 1:n, j in 1:n)
    for (i,j) in selected_e
        values_1[i,j] = Float64(I.th[i]+I.th[j]) 
    end
    
    weights_1 = Dict((i, j) => 1.0 for i in 1:n, j in 1:n)
    capacity_1 = Float64(I.T)
    
    val_1, delta_1 = fractional_knapsack(weights_1, values_1, capacity_1)

    values_2 = Dict((i, j) => 0. for i in 1:n, j in 1:n)
    for (i,j) in selected_e
        values_2[i,j] = Float64(I.th[i]*I.th[j]) 
    end
    
    weights_2 = Dict((i, j) => 1.0 for i in 1:n, j in 1:n)
    capacity_2 = Float64(I.T*I.T)
    
    val_2, delta_2 = fractional_knapsack(weights_2, values_2, capacity_2, 2.0)
    
    #obj = sum((I.distances[i, j] + delta_1[i, j] * (I.th[i] + I.th[j]) + delta_2[i, j] * (I.th[i] * I.th[j])) * x_bar[i, j] for i in 1:n, j in 1:n if i != j)
    
    obj = sum((I.distances[i, j]) * x_bar[i, j] for i in 1:n, j in 1:n if i != j)+val_1+val_2

    return (obj=obj, delta_1_opt=delta_1, delta_2_opt=delta_2)

end

function build_RCVRP_cuts_model(I::Instance)

    # Create the model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "TimeLimit", 5)


    # Decision variables
    @variable(model, z >= 0)
    @variable(model, x[i=1:n, j=1:n], Bin)  # Binary variable: 1 if arc (i, j) is used
    @variable(model, u[1:n] >= 0)           # MTZ variables: auxiliary variables for subtour elimination

    @constraint(model, [j in 2:n], sum(x[:, j]) == 1)  # Each node is visited once
    @constraint(model, [i in 2:n], sum(x[i, :]) == 1)  # Each node is visited once

    # MTZ subtour elimination constraints
    for i in 2:n, j in 2:n
        if i != j
            @constraint(model, u[i] - u[j] <= (I.capacity) * (1 - x[i, j]) - I.demands[j])
        end
    end

    for i in 2:n, j in 2:n
        if i != j
            @constraint(model, x[i, j] + x[j, i] <= 1)
        end
    end

    @constraint(model, [i in 2:n], I.demands[i] <= u[i] <= I.capacity)

    @constraint(model, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n))

    # Constraints
    #@constraint(model, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
    #@constraint(model, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot
    @constraint(model, sum(x[i, 1] for i in 2:n) >= ceil(sum(I.demands)/I.capacity) ) # 8 vehicles return to the depot

    #z constraint 
    @constraint(model, obj_cstr, sum(I.distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j) <= z)

    for i in 1:n
        @constraint(model, x[i, i] == 0)
    end

    @objective(model, Min, z)

    return model
end

model = build_RCVRP_cuts_model(instance)

function solve_RCVRP_iterative_cuts(I::Instance, model)

    MAXIMUM_ITERATIONS = 100
    #set_silent(model)

    optimize!(model)
    println(objective_value(model))

    @assert is_solved_and_feasible(model)

    for k in 1:MAXIMUM_ITERATIONS
        println("Iteration $(k)")
        #set_silent(model)


        @assert is_solved_and_feasible(model)
        lower_bound = objective_value(model)
        x_k = value.(model[:x])
        ret = solve_subproblem(I, x_k)

        println(lower_bound - ret.obj)
        for i in 1:n
            for j in 1:n
                if ret.delta_1_opt[i, j] > 0
                    println(x_k[i, j])
                end
            end
        end

        if lower_bound >= ret.obj
            break
        end

        one_route = find_one_route_clients(selected_edges(x_k, n), n)
        println(one_route)
        new_distances = zeros(n, n)
        for i in 1:n
            for j in 1:n
                new_distances[i, j] = I.distances[i, j] + ret.delta_1_opt[i, j] * (I.th[i] + I.th[j]) + ret.delta_2_opt[i, j] * (I.th[i] * I.th[j])

            end
        end

        new_distances_symetric = zeros(n, n)
        for i in 1:n
            for j in 1:n
                new_distances_symetric[i, j] = I.distances[j, i] + ret.delta_1_opt[j, i] * (I.th[i] + I.th[j]) + ret.delta_2_opt[j, i] * (I.th[i] * I.th[j])

            end
        end

        z = value.(model[:z])
        println("sum before : ", z - sum(new_distances[i, j] * x_k[i, j] for i in 1:n, j in 1:n if i != j))
        cut = @constraint(model, sum(new_distances[i, j] * model[:x][i, j] for i in 1:n, j in 1:n if i != j) <= model[:z])
        cut2 = @constraint(model, sum(new_distances_symetric[i, j] * model[:x][i, j] for i in 1:n, j in 1:n if i != j) <= model[:z])

        set_silent(model)
        optimize!(model)
        x_k = value.(model[:x])
        z = value.(model[:z])


        println("sum after : ", z - sum(new_distances[i, j] * x_k[i, j] for i in 1:n, j in 1:n if i != j))
        println("VRP value :", z)

    end
    return model
end

function solve_RCVRP_callback_cuts(I::Instance, model)

    function cut_callback(cb_data)
        status = callback_node_status(cb_data, model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        x_k = callback_value.(cb_data, model[:x])
        lower_bound = callback_value.(cb_data, model[:z])



        ret = solve_subproblem(I, x_k)
        if lower_bound >= ret.obj
            return
        end



        new_distances = zeros(n, n)
        for i in 1:n
            for j in 1:n
                new_distances[i, j] = I.distances[i, j] + ret.delta_1_opt[i, j] * (I.th[i] + I.th[j]) + ret.delta_2_opt[i, j] * (I.th[i] * I.th[j])
                #=                if (new_distances[j, i] < new_distances[i, j])
                                   new_distances[j, i] = new_distances[i, j] # symmetric 
                               end
                               if (new_distances[j, i] > new_distances[i, j])
                                   new_distances[i, j] = new_distances[j, i] # symmetric 
                               end =#
                #=   if ((x_k[i, j] > 0) && (ret.delta_1_opt[i, j] > 0 || ret.delta_2_opt[i, j] > 0))
                      println("Modifying distance ", i, j)
                  end =#
            end
        end

        new_distances_symetric = zeros(n, n)
        for i in 1:n
            for j in 1:n
                new_distances_symetric[i, j] = I.distances[j, i] + ret.delta_1_opt[j, i] * (I.th[i] + I.th[j]) + ret.delta_2_opt[j, i] * (I.th[i] * I.th[j])

            end
        end


        cut = @build_constraint(sum(new_distances[i, j] * model[:x][i, j] for i in 1:n, j in 1:n if i != j) <= model[:z])
        cut2 = @build_constraint( sum(new_distances_symetric[i, j] * model[:x][i, j] for i in 1:n, j in 1:n if i != j) <= model[:z])

        MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
        MOI.submit(model, MOI.LazyConstraint(cb_data), cut2)
    end

    set_attribute(
        model,
        MOI.LazyConstraintCallback(),
        cut_callback,
    )
    optimize!(model)
    return model

end

# Solve the model
solved_model = solve_RCVRP_iterative_cuts(instance, model)
#solved_model = solve_RCVRP_callback_cuts(instance, model)

# Extract and print the solution
if termination_status(solved_model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(solved_model))
    x_sol = value.(solved_model[:x])
    solution = [(i, j) for i in 1:n, j in 1:n if value(x_sol[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end


