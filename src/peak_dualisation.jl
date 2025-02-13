
using JuMP, CPLEX, LinearAlgebra

include("instance.jl")

instance  = read_instance("../data/n_9-euclidean_true")


function build_peak_dual_model(I::Instance)

    n = length(I.demands)
    # Create the model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "TimeLimit", 5)


    # Decision variables
    
    @variable(model, x[i=2:n, j=2:n], Bin)  # Binary variable: 1 if arc (i, j) is used
    @variable(model, x_0[j = 2:n]>=0, Int)
    @constraint(model,[j in 2:n], x_0[j] <=2)
   

    @variable(model, p[i = 2:n], Bin)
    @variable(model, t[i = 2:n]>=0, Int)  

    @variable(model, 0 <= mu_1[i=2:n, j=2:n])
    @variable(model, 0 <= mu_2[i=2:n, j=2:n])

    @variable(model, mu_0_1[j =2:n]>=0)
    @variable(model, mu_0_2[j =2:n]>=0)

    @variable(model, lambda_1 >= 0)
    @variable(model, lambda_2 >= 0)
   

    @variable(model, f[i=2:n, j=2:n]>=0, Int)
    
    @variable(model, u[i=2:n-1]>=0, Int)

    #Contraintes sur X
            
    @constraint(model, [i in 2:n], sum(x[:, i]) + x_0[i] - p[i] == 1) 

    @constraint(model, [i in 2:n], sum(x[i, :]) + p[i] == 1) 

    @constraint(model, [i in 2:n], t[i] + sum(f[i,:]) == sum(f[:,i]) +I.demands[i])

    @constraint(model, [i in 2:n, j in 2:n], I.demands[i]*x[i,j] <= f[i,j])

    @constraint(model, [i in 2:n, j in 2:n], (I.capacity - I.demands[j])*x[i,j] >= f[i,j])

    @constraint(model, [i in 2:n], I.demands[i]*p[i] <= t[i])

    @constraint(model, [i in 2:n], t[i] <= I.capacity*p[i])

    #@constraint(model, sum(p) >= ceil(sum(I.demands)/I.capacity))

    #@constraint(model, sum(x_0) >= 2*ceil(sum(I.demands)/I.capacity))

    @constraint(model, [i in 2:n-1], i <= u[i])

    @constraint(model, [i in 2:n-1], u[i] <= i*p[i] + (n-1)*(1-p[i]))
    
    #@constraint(model, [i in 2:n-1, j in 2:n-1, i!=j],u[i] - u[j] +(n-j-1)*x[i,j] <= n-j-1)
    
    @constraint(model, [i in 2:n], x[i,i]==0)

    #Rounded Peak Count Inequalities
    @constraint(model, [j in 2:n], sum(p[k] for k in j:n) >= ceil(sum(I.demands[k]/I.capacity for k in j:n)))

    #Contraintes qui renforcent la relaxation
    for i in 2:n-1
        for j in 2:n-1
            if i<= j
                @constraint(model, u[i] - u[j] + (n-j-1)*x[i,j] + (n-j-1)x[j,i] <= n-j-1)
            else
                @constraint(model, u[i] - u[j] + (n-j-1)*x[i,j] + (n-i-1)x[j,i] <= n-j-1)
            end
        end
    end
    
    @constraint(model, [i in 2:n-1, j in 2:n-1], u[i] - u[j] + (n-j-1)*x[i,j] + (n-i-1)*p[i] <= n-j-1)

    #Contraintes du DUAL
    for i in 2:n
        for j in 2:n
            if (i != j)
                @constraint(model, mu_1[i, j] + lambda_1 >= (I.th[i] + I.th[j]) * x[i, j])
                @constraint(model, mu_2[i, j] + lambda_2 >= (I.th[i] * I.th[j]) * x[i, j])
            end
        end
    end
    for j in 2:n
        @constraint(model, mu_0_1[j] + lambda_1 >= (I.th[1] + I.th[j]) * x_0[j])
        @constraint(model, mu_0_2[j] + lambda_2 >= (I.th[1] * I.th[j]) * x_0[j])
    end
    
    @objective(model, Min, sum(I.distances[i,j]*x[i,j] for i in 2:n, j in 2:n) + sum(x_0[j]*I.distances[1,j] for j in 2:n) + 
                                lambda_1 * I.T + lambda_2 * I.T^2 + sum(mu_1 + 2 * mu_2)+ sum(mu_0_1 + 2*mu_0_2))

    return model
end

function build_peak_dual_model_test(I::Instance)

    n = length(I.demands)
    # Create the model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "TimeLimit", 5)


    # Decision variables
    
    @variable(model, x[i=1:n, j=1:n]>=0, Int)  # Binary variable: 1 if arc (i, j) is used
    
    @constraint(model,[j in 2:n], x[1,j] <=2)
    @constraint(model,[j in 2:n], x[j,1] <=1)
    @constraint(model,[i = 2:n,j=2:n], x[i,j] <=1)

    @variable(model, p[i = 2:n], Bin)
    @variable(model, t[i = 2:n]>=0, Int)  

    @variable(model, 0 <= mu_1[i=1:n, j=1:n])
    @variable(model, 0 <= mu_2[i=1:n, j=1:n])

    @variable(model, lambda_1 >= 0)
    @variable(model, lambda_2 >= 0)
   

    @variable(model, f[i=2:n, j=2:n]>=0, Int)
    
    @variable(model, u[i=2:n-1]>=0, Int)

    #Contraintes sur X
            
    @constraint(model, [i in 2:n], sum(x[:, i]) - p[i] == 1) 

    @constraint(model, [i in 2:n], sum(x[i, :]) + p[i] == 1) 

    @constraint(model, [i in 2:n], t[i] + sum(f[i,:]) == sum(f[:,i]) +I.demands[i])

    @constraint(model, [i in 2:n, j in 2:n], I.demands[i]*x[i,j] <= f[i,j])

    @constraint(model, [i in 2:n, j in 2:n], (I.capacity - I.demands[j])*x[i,j] >= f[i,j])

    @constraint(model, [i in 2:n], I.demands[i]*p[i] <= t[i])

    @constraint(model, [i in 2:n], t[i] <= I.capacity*p[i])

    #@constraint(model, sum(p) >= ceil(sum(I.demands)/I.capacity))

    #@constraint(model, sum(x_0) >= 2*ceil(sum(I.demands)/I.capacity))

    @constraint(model, [i in 2:n-1], i <= u[i])

    @constraint(model, [i in 2:n-1], u[i] <= i*p[i] + (n-1)*(1-p[i]))
    
    #@constraint(model, [i in 2:n-1, j in 2:n-1, i!=j],u[i] - u[j] +(n-j-1)*x[i,j] <= n-j-1)
    
    @constraint(model, [i in 2:n], x[i,i]==0)

    #Rounded Peak Count Inequalities
    @constraint(model, [j in 2:n], sum(p[k] for k in j:n) >= ceil(sum(I.demands[k]/I.capacity for k in j:n)))

    #Contraintes qui renforcent la relaxation
    for i in 2:n-1
        for j in 2:n-1
            if i<= j
                @constraint(model, u[i] - u[j] + (n-j-1)*x[i,j] + (n-j-1)x[j,i] <= n-j-1)
            else
                @constraint(model, u[i] - u[j] + (n-j-1)*x[i,j] + (n-i-1)x[j,i] <= n-j-1)
            end
        end
    end
    
    @constraint(model, [i in 2:n-1, j in 2:n-1], u[i] - u[j] + (n-j-1)*x[i,j] + (n-i-1)*p[i] <= n-j-1)

    #Contraintes du DUAL
    for i in 1:n
        for j in 1:n
            if (i != j)
                @constraint(model, mu_1[i, j] + lambda_1 >= (I.th[i] + I.th[j]) * x[i, j])
                @constraint(model, mu_2[i, j] + lambda_2 >= (I.th[i] * I.th[j]) * x[i, j])
            end
        end
    end
    for j in 1:n
        @constraint(model, mu_1[j,1] == 0)
        @constraint(model, mu_2[j,1] == 0)
    end

    
    
    @objective(model, Min, sum(I.distances[i,j]*x[i,j] for i in 1:n, j in 1:n) + lambda_1 * I.T + lambda_2 * I.T^2 + sum(mu_1 + 2 * mu_2))

    return model
end




# Create the model
model = build_peak_dual_model_test(instance)
#set_optimizer_attribute(model, "TimeLimit", 5)



optimize!(model)

# Solve the model

# Extract and print the solution
if termination_status(model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model))
    x_sol = value.(model[:x])
    #x0_sol = value.(model[:x_0])
    p_sol = value.(model[:p])
    t_sol = value.(model[:t])

    #solution_x = [(i, j) for i in 2:n, j in 2:n if value(x_sol[i, j]) > 0.5]
    #solution_x0 =  [(1,j) for j in 2:n if value(x0_sol[j]) > 0.5]
    #println(solution_x, solution_x0)
    solution_x = [(i, j) for i in 1:n, j in 1:n if value(x_sol[i, j]) > 0.5]
    println(solution_x)
    println("Les sommets peak sont: ", [i for i in 2:n if value(p_sol[i]) > 0.5])
    println("Les valeurs de t sont : ",[value(t_sol[i]) for i in 2:n])

else
    println("No optimal solution found.")
end