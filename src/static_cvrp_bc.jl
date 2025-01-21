using JuMP, CPLEX, LinearAlgebra

file_path = "../data/P-n20-k2.vrp"  # Replace with your instance file path
name, comment, capacity, n, coords, demands, depot, distances = parse_cvrp_instance(file_path)

k = 2

# Create the model
model = Model(CPLEX.Optimizer)

@variable(model, x[i=1:n, j=(i+1):n], Int)  # Binary variable: 1 if arc (i, j) is used

@constraint(model, sum(x[1, :]) == 2 * k)
@constraint(model, visit[i=1:n], sum(x[i, :]) + sum(x[:, i]) == 2)

for i in 2:n
    @constraint(model, 0 <= x[1, i] <= 2)
end

for i in 2:n
    for j in (i+1):n
        @constraint(model, 0 <= x[i, j] <= 1)

    end
end


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
