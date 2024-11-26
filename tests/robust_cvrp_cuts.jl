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


function solve_subproblem(x_bar)
    sub_model = Model(CPLEX.Optimizer)

    @variable(sub_model, 0 <= delta_1[i=1:n, j=1:n] <= 1)
    @variable(sub_model, 0 <= delta_2[i=1:n, j=1:n] <= 2)

    @constraint(sub_model, sum(delta_1[i, j] for j in 1:n for i in 1:n) <= T)
    @constraint(sub_model, sum(delta_2[i, j] for j in 1:n for i in 1:n) <= T * T)


    for i in 1:n
        @constraint(sub_model, delta_1[i, i] == 0)
        @constraint(sub_model, delta_2[i, i] == 0)
    end

    @objective(sub_model, Max, sum((distances[i, j] + delta_1[i, j] * (IT[i] + IT[j]) + delta_2[i, j] * (IT[i] * IT[j])) * x_bar[i, j] for i in 1:n, j in 1:n if i != j))
    optimize!(sub_model)
    @assert is_solved_and_feasible(sub_model; dual=true)
    return (obj=objective_value(sub_model), delta_1_opt=value.(delta_1),
        delta_2_opt=value.(delta_2))
end




# Create the model
model = Model(CPLEX.Optimizer)
#set_optimizer_attribute(model, "TimeLimit", 5)


# Decision variables
@variable(model, z >= 0)
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





# Constraints
@constraint(model, sum(x[1, j] for j in 2:n) == 2)  # 8 vehicles leave the depot
@constraint(model, sum(x[i, 1] for i in 2:n) == 2)  # 8 vehicles return to the depot

#z constraint 
@constraint(model, obj_cstr,sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n if i != j) <= z)


@objective(model, Min, z)

for i in 1:n
    @constraint(model, x[i, i] == 0)
end
#optimize!(model)

MAXIMUM_ITERATIONS = 100
#set_silent(model)
optimize!(model)
@assert is_solved_and_feasible(model)


println("Iteration  Lower Bound  Upper Bound          Gap")
@assert is_solved_and_feasible(model)

for k in 1:MAXIMUM_ITERATIONS
    println("Iteration $(k)")
    #set_silent(model)


    @assert is_solved_and_feasible(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(x_k)
    #print_iteration(k, lower_bound, upper_bound, gap)
    #=   if gap < ABSOLUTE_OPTIMALITY_GAP
          println("Terminating with the optimal solution")
          break
      end =#
    println(lower_bound - ret.obj)
    for i in 1:n
        for j in 1:n
            if ret.delta_1_opt[i, j] > 0
                println(x_k[i, j])
            end
        end
    end
    #println(ret.delta_1_opt)
    #println(ret.delta_2_opt)
    if lower_bound >= ret.obj
        break
    end


    new_distances = zeros(n, n)
    for i in 1:n
        for j in 1:n
            new_distances[i, j] = distances[i, j] + ret.delta_1_opt[i, j] * (IT[i] + IT[j]) + ret.delta_2_opt[i, j] * (IT[i] * IT[j])
            #new_distances[j,i] = new_distances[i,j]
            #set_normalized_coefficient(obj_cstr,x[i,j],new_distances[i,j])
            if ((x_k[i,j] > 0) && (ret.delta_1_opt[i,j] > 0 || ret.delta_2_opt[i,j] > 0))
                println("Modifying distance ", i, j )
            end
            #println(new_distances[i, j] - distances[i, j],x_k[i,j])
        end
    end
    println("sum before : ", value.(z) - sum(new_distances[i, j]* x_k[i, j] for i in 1:n, j in 1:n if i != j))
    cut = @constraint(model, sum(new_distances[i, j]* x[i, j] for i in 1:n, j in 1:n if i != j) <= z)

    set_silent(model)
    optimize!(model)
    x_k = value.(x)

    println("sum after : ", value.(z) - sum(new_distances[i, j]* x_k[i, j] for i in 1:n, j in 1:n if i != j))
    println("VRP value :", value.(z))


    #cut = @constraint(model, sum((distances[i, j] *10) * x[i, j] for i in 1:n, j in 1:n if i != j) <= z)



    #@constraint(model,z >= 300)
    #@info "Adding the cut $(cut)"
end




# Solve the model

# Extract and print the solution
if termination_status(model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model))
    solution = [(i, j) for i in 1:n, j in 1:n if value(x[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end