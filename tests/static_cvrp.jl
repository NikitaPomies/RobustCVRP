
using JuMP, CPLEX, LinearAlgebra

file_path = "tests/data/P-n20-k2.vrp"  # Replace with your instance file path
name, comment, capacity, n, coords, demands, depot, distances = parse_cvrp_instance(file_path)


# Create the model
model = Model(CPLEX.Optimizer)
#set_optimizer_attribute(model, "TimeLimit", 5)


# Decision variables
@variable(model, x[i=1:n, j=1:n], Bin)  # Binary variable: 1 if arc (i, j) is used
@variable(model, u[1:n] >= 0)           # MTZ variables: auxiliary variables for subtour elimination

# Objective: minimize the total distance
#@objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))

for j in 2:n
    @constraint(model, sum(x[i, j] for i in 1:n if i != j) == 1)  # Each node is visited once
end

for i in 2:n
    @constraint(model, sum(x[i, j] for j in 1:n if i != j) == 1)  # Each node is visited once

end

# MTZ subtour elimination constraints
for i in 2:n, j in 2:n
    if i != j
        @constraint(model, u[i] - u[j] <= (capacity) * (1 - x[i, j]) - demands[j])
    end
end

for i in 2:n, j in 2:n
    if i != j
        @constraint(model, x[i, j] + x[j, i] <= 1)
    end
end

for i in 2:n
    @constraint(model, demands[i] <= u[i] <= capacity)

end

@constraint(model, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n)) # 8 vehicles leave the depot


#= 
@objective(model, Min, sum(x[1, j] for j in 2:n ))
optimize!(model)
first_obj_value = objective_value(model)
prinln(objective_value(model))


 =#

# Constraints
@constraint(model, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
@constraint(model, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot

@objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))
optimize!(model)

# Solve the model

# Extract and print the solution
if termination_status(model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model))
    solution = [(i, j) for i in 1:n, j in 1:n if value(x[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end