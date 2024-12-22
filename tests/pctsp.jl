
using JuMP
import CPLEX
import Random
using Plots

function generate_distance_matrix(n; random_seed = 6)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    prices  = 100 * rand(rng,n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return X, Y, d , prices 
end

n = 20
X, Y, d, prices = generate_distance_matrix(n)
println(prices)


function build_tsp_model(d, n)
    model = Model(CPLEX.Optimizer)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @variable(model, z[1:n], Bin)
    @objective(model, Max, sum(prices .* z) - sum(d .* x) / 2)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 2*z[i])
    @constraint(model, [i in 1:n], x[i, i] == 0)
    @constraint(model,z[1]==1) # En foncuton du prix associé au dépot, il faut forcer le pasage par le dépot 
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

lazy_model = build_tsp_model(d, n)




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

set_attribute(
    lazy_model,
    MOI.LazyConstraintCallback(),
    subtour_elimination_callback,
)
optimize!(lazy_model)


function plot_tour(X, Y, x)
    # Create empty plot
    plot = Plots.plot(size=(800,600))
    
    # Plot the edges (routes)
    for (i, j) in selected_edges(x, size(x, 1))
        Plots.plot!([X[i], X[j]], [Y[i], Y[j]], 
                   color=:blue, 
                   linewidth=2, 
                   legend=false)
    end
    
    # Plot the nodes
    Plots.scatter!(X, Y, 
                  color=:red, 
                  markersize=8, 
                  legend=false)
    
    # Add node numbers
    for i in 1:length(X)
        # Annotate with node numbers
        # offset the text slightly above the point for better visibility
        price = prices[i]
        Plots.annotate!(X[i], Y[i]+1, Plots.text("$i price $price", :black, :center, 8))
    end
    
    # Highlight depot (node 1) differently
    Plots.scatter!([X[1]], [Y[1]], 
                  color=:green, 
                  markersize=10, 
                  legend=false)
    
    # Set axis labels
    xlabel!("X Coordinate")
    ylabel!("Y Coordinate")
    title!("CVRP Solution with Node Numbers")
    
    return plot
end

plot_tour(X, Y, value.(lazy_model[:x]))