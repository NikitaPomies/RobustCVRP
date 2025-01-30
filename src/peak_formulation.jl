using JuMP, CPLEX, LinearAlgebra


include("instance.jl")

instance = read_instance("../data/n_6-euclidean_true")

function build_BPF_model(I::Instance)

    # Create the model
    model = Model(CPLEX.Optimizer)
    K=n
    #set_optimizer_attribute(model, "TimeLimit", 5)


    # Decision variables

    @variable(model, x[i=2:n, j=2:n], Bin)  # Binary variable: 1 if arc (i, j) is used
    @variable(model, p[i = 2:n], Bin)
    @constraint(model, p[n] == 1)

    @variable(model, t[i = 2:n], Int)
    @variable(model, x_0[j = 2:n], Int)
    @constraint(model,[j in 2:n], x_0[j] <=2)

    @variable(model, f[i=1:n, j=1:n], Int)
    @constraint(model, [j in 1:n], f[1,j] ==0)

    @variable(model, u[i=2:n-1], Int)
            

    @constraint(model, [i in 2:n], sum(x[:, i]) - p[i]== 1)  
    @constraint(model, [i in 2:n], sum(x[i, :]) + p[i]== 1) 

    @constraint(model, [i in 2:n], t[i] + sum(f[i,:]) == sum(f[:,i]) +I.demands[i])

    @constraint(model, [i in 2:n, j in 2:n], I.demands[i]*x[i,j] <= f[i,j])

    @constraint(model, [i in 2:n, j in 2:n], (I.capacity - I.demands[j])*x[i,j] >= f[i,j])

    @constraint(model, [i in 2:n], I.demands[i]*p[i] <= t[i])
    @constraint(model, [i in 2:n], t[i] <= I.capacity*p[i])

    #@constraint(model, sum(x_0) <= 2*K)
    #constraint(model, sum(p) <= K)
    @constraint(model, [i in 2:n-1], i <= u[i])
    @constraint(model, [i in 2:n-1], u[i] <= i*p[i] + (n-1)*(1-p[i]))

    # MTZ subtour elimination constraints
    for i in 2:n-1, j in 2:n-1
        if i != j
            @constraint(model, u[i] - u[j] +(n-j-1)*x[i,j] <= n-j-1)
        end
    end

    @constraint(model, [i in 2:n], x[i,i]==0)
    
    @objective(model, Min, 
    sum(I.distances[i,j]*x[i,j] for i in 2:n, j in 2:n) + sum(x_0[j]*I.distances[1,j] for j in 2:n))

    return model
end

model = build_BPF_model(instance)

optimize!(model)
println("Nombre de variables x: ", length(all_variables(model)))


if termination_status(model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model))
    x_sol = value.(model[:x])
    
    solution = [(i, j) for i in 2:n, j in 2:n if value(x_sol[i, j]) > 0.5]
    println(solution)

else
    println("No optimal solution found.")
end

