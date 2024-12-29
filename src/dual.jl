using JuMP, CPLEX, LinearAlgebra

#Script pour tester la modélisation du problème interne pour cvrp robuste


n = 10

# Define the size of the vector and matrix
n = 10  # You can change this to any integer value
T = 5

# Generate a vector of size n with random integers between 0 and 100
IT = rand(0:10, n)

# Generate a matrix of size n x n with random integers between 0 and 100
matrix = rand(0:10, n, n)

# Print the results
println("Vector IT:")
println(IT)

println("\nMatrix:")
println(matrix)

model_p = Model(CPLEX.Optimizer)

@variable(model_p, 0<=delta_1[i=1:n, j=1:n] <=1)
@variable(model_p,  0<=delta_2[i=1:n, j=1:n] <= 2)

@objective(model_p, Max, sum((delta_1[i,j]*(IT[i] + IT[j]) + delta_2[i,j]*(IT[i]*IT[j]) )* matrix[i, j] for i in 1:n, j in 1:n ))

@constraint(model_p, sum(delta_1[i, j] for j in 1:n for i in 1:n) <=T)  
@constraint(model_p, sum(delta_2[i, j] for j in 1:n for i in 1:n) <=T*T)  



optimize!(model_p)


println("Optimal objective value primal: ", objective_value(model_p))



model_d = Model(CPLEX.Optimizer)

@variable(model_d, 0<=mu_1[i=1:n, j=1:n] )
@variable(model_d,  0<=mu_2[i=1:n, j=1:n] )
@variable(model_d,lambda_1>=0)
@variable(model_d,lambda_2>=0)

@objective(model_d,Min,lambda_1 *T + lambda_2*T*T + sum(mu_1[i,j] + 2*mu_2[i,j] for i in 1:n for j in 1:n))

for i in 1:n
    for j in 1:n

    @constraint(model_d,mu_1[i,j] + lambda_1 >= (IT[i]+IT[j])*matrix[i,j])
    @constraint(model_d,mu_2[i,j] + lambda_2 >= (IT[i]*IT[j])*matrix[i,j])
    end
end

optimize!(model_d)

println("Optimal objective value dual: ", objective_value(model_d))




