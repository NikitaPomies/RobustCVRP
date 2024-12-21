
using JuMP
import CPLEX
import Random
using Plots

function generate_distance_matrix(n; random_seed = 1)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    prices  = 100 * rand(rng,n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return X, Y, d , prices 
end

n = 100
X, Y, d, prices = generate_distance_matrix(n)



function build_tsp_model(d, n)
    model = Model(CPLEX.Optimizer)
    @variable(model, x[1:n, 1:n], Bin, Symmetric)
    @objective(model, Min, sum(d .* x) / 2)
    @constraint(model, [i in 1:n], sum(x[i, :]) == 2)
    @constraint(model, [i in 1:n], x[i, i] == 0)
    return model
end

function selected_edges(x::Matrix{Float64}, n)
    return Tuple{Int,Int}[(i, j) for i in 1:n, j in 1:n if x[i, j] > 0.5]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x, size(x, 1)), size(x, 1))
subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))

function subtour(edges::Vector{Tuple{Int,Int}}, n)
    shortest_subtour, unvisited = collect(1:n), Set(collect(1:n))
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
        return  # Only run at integer solutions
    end
    cycle = subtour(callback_value.(cb_data, lazy_model[:x]))
    if !(1 < length(cycle) < n)
        return  # Only add a constraint if there is a cycle
    end
    S = [(i, j) for (i, j) in Iterators.product(cycle, cycle) if i < j]
    con = @build_constraint(
        sum(lazy_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
    )
    MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), con)
    return
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
        Plots.annotate!(X[i], Y[i]+1, Plots.text("$i", :black, :center, 8))
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