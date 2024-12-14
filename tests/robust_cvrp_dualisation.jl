
using JuMP, CPLEX, LinearAlgebra

file_path = "tests/data/P-n20-k2.vrp"  # Replace with your instance file path
name, comment, capacity, n, coords, demands, depot, distances = parse_cvrp_instance(file_path)

T = 1

IT = rand(0:10, n)

# Generate a matrix of size n x n with random integers between 0 and 100
matrix = rand(0:10, n, n)

# Print the results
println("Vector IT:")
println(IT)

println("\nMatrix:")
println(matrix)


# Create the model
model = Model(CPLEX.Optimizer)
#set_optimizer_attribute(model, "TimeLimit", 5)


# Decision variables
@variable(model, x[i=1:n, j=1:n], Bin)  # Binary variable: 1 if arc (i, j) is used
@variable(model, u[1:n] >= 0)           # MTZ variables: auxiliary variables for subtour elimination
@variable(model, 0 <= mu_1[i=1:n, j=1:n])
@variable(model, 0 <= mu_2[i=1:n, j=1:n])
@variable(model, lambda_1 >= 0)
@variable(model, lambda_2 >= 0)


for i in 1:n
    for j in 1:n
        if (i != j)
            @constraint(model, mu_1[i, j] + lambda_1 >= (IT[i] + IT[j]) * x[i, j])
            @constraint(model, mu_2[i, j] + lambda_2 >= (IT[i] * IT[j]) * x[i, j])
        end
    end
end


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





# Constraints
@constraint(model, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
@constraint(model, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot

@objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j) + 
lambda_1 * T + lambda_2 * T * T + sum(mu_1[i, j] + 2 * mu_2[i, j] for i in 1:n for j in 1:n if i != j)) 
#@objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))

for i in 1:n
    @constraint(model,x[i,i]==0)
    @constraint(model,mu_1[i,i]==0)
    @constraint(model,mu_2[i,i]==0)
end

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