
using JuMP, CPLEX, LinearAlgebra

include("instance.jl")

instance  = read_instance("../data/n_12-euclidean_true")


function build_RCVRP_dual_model(I::Instance)
    n = length(I.demands)
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
                @constraint(model_d, mu_1[i, j] + lambda_1 >= (I.th[i] + I.th[j]) * x[i, j])
                @constraint(model_d, mu_2[i, j] + lambda_2 >= (I.th[i] * I.th[j]) * x[i, j])
            end
        end
    end

    @constraint(model_d, [j in 2:n],sum(x[:, j]) == 1)  # Each node is visited once

    @constraint(model_d, [i in 2:n], sum(x[i, :]) == 1)  # Each node is visited once

    # MTZ subtour elimination constraints
    for i in 2:n, j in 2:n
        if i != j
            @constraint(model_d, u[i] - u[j] <= (I.capacity) * (1 - x[i, j]) - I.demands[j])
        end
    end

    for i in 2:n, j in 2:n
        if i != j
            @constraint(model_d, x[i, j] + x[j, i] <= 1)
        end
    end

   #=  # For any three nodes i, j, k
    for i in 2:n, j in 2:n,k in 2:n
        @constraint(model_d, x[i,j] + x[j,k] + x[k,i] <= 2)
    end =#

    @constraint(model_d,[i in 2:n], I.demands[i] <= u[i] <= I.capacity)


    @constraint(model_d, sum(x[1, j] for j in 2:n) == sum(x[j, 1] for j in 2:n)) # 8 vehicles leave the depot

    # Constraints
    #@constraint(model_d, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
    #@constraint(model_d, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot
    @constraint(model_d,sum(x[1, j] for j in 2:n) >= ceil(sum(I.demands)/I.capacity))


    @objective(model_d, Min, sum(I.distances .*x) + lambda_1 * I.T + lambda_2 * I.T^2 + sum(mu_1 + 2 * mu_2))

    for i in 1:n
        @constraint(model_d, x[i, i] == 0)
        @constraint(model_d, mu_1[i, i] == 0)
        @constraint(model_d, mu_2[i, i] == 0)
    end

    return model_d
end



# Create the model_d
model_d = build_RCVRP_dual_model(instance)
#set_optimizer_attribute(model_d, "TimeLimit", 5)



optimize!(model_d)

# Solve the model_d

# Extract and print the solution
if termination_status(model_d) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model_d))
    x_sol = value.(model_d[:x])
    solution = [(i, j) for i in 1:n, j in 1:n if value(x_sol[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end