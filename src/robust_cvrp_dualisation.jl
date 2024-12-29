
using JuMP, CPLEX, LinearAlgebra

include("../data/n_13-euclidean_true")

println("n = $n")
println("t = $t")
println("th = $th")
println("T = $T")
println("d = $d")
println("C = $C")

demands = d
distances = t
capacity = C

function build_RCVRP_dual_model(th,capacity,distances,demands)
    

end



# Create the model_d
model_d = Model(CPLEX.Optimizer)
#set_optimizer_attribute(model_d, "TimeLimit", 5)


# Decision variables
@variable(model_d, x[i=1:n, j=1:n], Bin)  # Binary variable: 1 if arc (i, j) is used
@variable(model_d, u[1:n] >= 0)           # MTZ variables: auxiliary variables for subtour elimination
@variable(model_d, 0 <= mu_1[i=1:n, j=1:n])
@variable(model_d, 0 <= mu_2[i=1:n, j=1:n])
@variable(model_d, lambda_1 >= 0)
@variable(model_d, lambda_2 >= 0)


for i in 1:n
    for j in 1:n
        if (i != j)
            @constraint(model_d, mu_1[i, j] + lambda_1 >= (th[i] + th[j]) * x[i, j])
            @constraint(model_d, mu_2[i, j] + lambda_2 >= (th[i] * th[j]) * x[i, j])
        end
    end
end


# Objective: minimize the total distance
#@objective(model_d, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))

for j in 2:n
    @constraint(model_d, sum(x[i, j] for i in 1:n if i != j) == 1)  # Each node is visited once
end

for i in 2:n
    @constraint(model_d, sum(x[i, j] for j in 1:n if i != j) == 1)  # Each node is visited once

end

# MTZ subtour elimination constraints
for i in 2:n, j in 2:n
    if i != j
        @constraint(model_d, u[i] - u[j] <= (capacity) * (1 - x[i, j]) - demands[j])
    end
end

for i in 2:n, j in 2:n
    if i != j
        @constraint(model_d, x[i, j] + x[j, i] <= 1)
    end
end

for i in 2:n
    @constraint(model_d, demands[i] <= u[i] <= capacity)

end

@constraint(model_d, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n)) # 8 vehicles leave the depot





# Constraints
#@constraint(model_d, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
#@constraint(model_d, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot

@objective(model_d, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j) +
                         lambda_1 * T + lambda_2 * T * T + sum(mu_1[i, j] + 2 * mu_2[i, j] for i in 1:n for j in 1:n if i != j))
#@objective(model_d, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j))

for i in 1:n
    @constraint(model_d, x[i, i] == 0)
    @constraint(model_d, mu_1[i, i] == 0)
    @constraint(model_d, mu_2[i, i] == 0)
end

optimize!(model_d)

# Solve the model_d

# Extract and print the solution
if termination_status(model_d) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model_d))
    solution = [(i, j) for i in 1:n, j in 1:n if value(x[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end