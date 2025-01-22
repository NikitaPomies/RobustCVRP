
using Random
using JuMP
using Gurobi
using LinearAlgebra


function generate_distance_matrix(n; random_seed=2)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    prices = 100 * rand(rng, n)
    demands = [round(100 * rand(rng), digits=1) for i in 1:n]
    d = [round(sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2),digits = 1) for i in 1:n, j in 1:n]
    return X, Y, d, prices, demands
end

function update_distance_matrix!(d::Matrix{Float64})
    # Generate a random α_i for each row
    alpha = rand(0:80, size(d, 1))
    
    # Update the matrix
    for i in 1:size(d, 1), j in 1:size(d, 2)
        d[i, j] -= alpha[i]
    end

    return d
end

n = 20
X, Y, d, prices, demands = generate_distance_matrix(n-1)
capacity = 15
d= update_distance_matrix!(d)
# Adding a new row (equal to the first row)
new_row = d[1, :]  # Extract the first row
matrix_with_row = vcat(d, new_row')

# Adding a new column (equal to the first column)
new_column = d[:, 1]  # Extract the first column
d_dp = hcat(matrix_with_row, [new_column; new_column[1]])

demands_dp = copy(demands)
push!(demands_dp,0)
demands_dp[1] = 0
d_dp[1,n] = 1000000
d_dp[n,1] = 1000000







mutable struct Label
    dummy_ressource::Vector{Int}
    ressources::Vector{Any}
    cost::Any
    current_node::Int
end



function Base.copy(label::Label)
    return Label(
        copy(label.dummy_ressource),
        copy(label.ressources),
        label.cost,
        label.current_node
    )
end

function extend(label::Label, next_node::Int, distances::Matrix{T}, demands::Vector{T}, capacity::Int) where {T <: Real}
    #Create new label
    next_node_consumption = label.ressources[end] + demands[next_node]
    if next_node_consumption > capacity
        return nothing
    end
    new_label = copy(label)
    new_label.dummy_ressource[next_node] += 1
    new_label.cost += distances[label.current_node, next_node]
    new_label.current_node = next_node
    push!(new_label.ressources, next_node_consumption)
    return new_label
end




function dominates(label1::Label, label2::Label)
    cost_dominance = label1.cost <= label2.cost
    ressource_dominance = label1.ressources[end] <= label2.ressources[end]
    #dummy_ressource_dominance = label1.dummy_ressource[end] <= label2.dummy_ressource
    return all([cost_dominance, ressource_dominance])
end

function add_non_dominated!(label_liste::Vector{Label}, new_label::Label)
    # Check if the new label is dominated by any existing label
    if any(dominates(label, new_label) for label in label_liste)
        return false,label_liste  # The new label is dominated; no changes to label_liste
    end

    # Remove labels that are dominated by the new label
    label_liste = filter(label -> !dominates(new_label, label), label_liste)

    # Add the new label to the list
    push!(label_liste, new_label)

    return true,label_liste
end


function mono_direction_dp(distances,demands,capacity)
    src = 1
    dest = n
    label_lists = [Vector{Label}() for _ in 1:n]
    push!(label_lists[src],Label(zeros(Int,n),Int[0],0,src))
    label_lists[src][1].dummy_ressource[1] = 1
    E = [src]
    #print(label_lists)
    while length(E) > 0
        current_node = pop!(E)
        for label in label_lists[current_node]

            for j in 1:n
                if j!=current_node && label.dummy_ressource[j] == 0   #complete graph

                    #println(label)
                    new_label = extend(label,j,distances,demands,capacity)
                    if new_label === nothing
                        continue
                    else
                        changed,label_lists[j] = add_non_dominated!(label_lists[j],new_label)
                    end
                    if changed
                        push!(E,j)
                    end

                end
            end
        end
        #println(E)
        filter!(e->e!=current_node,E)


    end


return label_lists
end

println(d_dp)
println()
test = mono_direction_dp(d_dp,demands_dp,200)
#println(test[n])
best_route = findmin(x.cost for x in test[n])
println(best_route)
#println(d[1,10])
#println(demands)

capacity  = 200



function build_tsp_model(demandes, capa, n)
    model = Model(Gurobi.Optimizer)
    #set_optimizer_attribute(model, "CPX_PARAM_PREIND", 0)
    #set_optimizer_attribute(model, "CPX_PARAM_ADVIND", 0)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @variable(model, z[1:n], Bin)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 2 * z[i])
    @constraint(model, [i in 1:n], x[i, i] == 0)
    @constraint(model, z[1] == 1) # En foncuton du prix associé au dépot, il faut forcer le pasage par le dépot 
    @constraint(model, sum(demandes .* z) <= capa)
    return model
end



function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

function subtour(x::Matrix{Float64}, z::Vector{Float64})
    return subtour(selected_edges(x, size(x, 1)), size(x, 1), z)
end

subtour(x::AbstractMatrix{VariableRef}, z::Vector{VariableRef}) = subtour(value.(x), value.(z))

function subtour(edges::Vector{Tuple{Int,Int}}, n, z::Vector{Float64})
    selected = Set([i for i in 1:n if z[i] > 0.5])

    shortest_subtour = collect(selected)  # Initialize with all selected nodes
    unvisited = copy(selected)  # Start with selected nodes only

    while !isempty(unvisited)
        this_cycle, neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors =
                [j for (i, j) in edges if i == current && j in unvisited]
        end
        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    return shortest_subtour
end

#lazy_model = build_tsp_model(demands, capacity, n)

function subtour_elimination_callback(cb_data, model)
    status = callback_node_status(cb_data, model)
    if status != MOI.CALLBACK_NODE_STATUS_INTEGER
        return
    end
    
    x_val = callback_value.(cb_data, model[:x])
    y_val = callback_value.(cb_data, model[:z])
    cycle = subtour(x_val, y_val)
    selected_count = sum(y_val)
    
    if !(1 < length(cycle) < selected_count)
        return
    end
    
    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    con = @build_constraint(
        sum(model[:x][i, j] for (i, j) in S) <= length(cycle) - 1
    )
    MOI.submit(model, MOI.LazyConstraint(cb_data), con)
end


function subtour_elimination_callback(cb_data)
    status = callback_node_status(cb_data, lazy_model)
    if status != MOI.CALLBACK_NODE_STATUS_INTEGER
        return
    end

    x_val = callback_value.(cb_data, lazy_model[:x])
    y_val = callback_value.(cb_data, lazy_model[:z])

    cycle = subtour(x_val, y_val)
    selected_count = sum(y_val)

    if !(1 < length(cycle) < selected_count)
        return
    end

    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    con = @build_constraint(
        sum(lazy_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1
    )
    MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
end
#= 
set_attribute(
    lazy_model,
    MOI.LazyConstraintCallback(),
    subtour_elimination_callback,
) =#

function solvepctsp(prices::Vector{T}, model::Model,dist) where {T<:Real}
    
    @objective(model, Min, -sum(prices .* model[:z]) + sum(dist .* model[:x]) / 2)
    optimize!(model)
end

function get_route(x::Matrix{Float64}, z::Vector{Float64})
    n = size(x, 1)
    selected = [i for i in 1:n if ((z[i] > 0.5) && (i!=1))]
    route = [1]  # Start from depot
    current = 1

    
    while length(route) <= length(selected) +1
        # Find next node in route
        next_found = false
        for j in selected
            if !( j in route) && x[current, j] > 0.5
                push!(route, j)
                current = j
                next_found = true
                break
            end
        end
        
        # If no next node found (means we're returning to depot)
        # and we haven't visited all selected nodes, start from depot again
        if !next_found 
            push!(route,1)
            break
        end
    end
    
    return route
end

demands[1] = 0
lazy_model = build_tsp_model(demands, capacity, n-1)
set_silent(lazy_model)


set_attribute(
    lazy_model,
    MOI.LazyConstraintCallback(),
    subtour_elimination_callback,
) 


solvepctsp(zeros(n-1),lazy_model,d)
route = get_route(value.(lazy_model[:x]),value.(lazy_model[:z]))

global S = 0
for x in route
    global S+= demands[x]
end
println("capacity : $(S)")
println("cost $(objective_value(lazy_model))")
println(route)

println(d -transpose(d))
