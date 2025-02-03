using JuMP, CPLEX, LinearAlgebra, GLPK


include("instance.jl")

instance = read_instance("../data/n_6-euclidean_true")
"""
function build_BPF_model(I::Instance)

    # Create the model
    model = Model(CPLEX.Optimizer)
    K=n
    #set_optimizer_attribute(model, "TimeLimit", 5)


    # Decision variables

    @variable(model, x[i=2:n, j=2:n], Bin)  # Binary variable: 1 if arc (i, j) is used
    @variable(model, p[i = 2:n], Bin)
    #@constraint(model, p[n] == 1)

    @variable(model, t[i = 2:n], Int)
    @constraint(model, [i in 2:n], t[i] >=0)

    @variable(model, x_0[j = 2:n], Int)
    @constraint(model,[j in 2:n], x_0[j] <=2)
    @constraint(model, [j in 2:n], x_0[j] >=0)

    @variable(model, f[i=2:n, j=2:n], Int)
    #@constraint(model, [j in 1:n], f[1,j] ==0)
    @constraint(model, [i=2:n, j=2:n], f[i,j] >=0)

    @variable(model, u[i=2:n-1], Int)
            
    @constraint(model, [i in 2:n], sum(x[:, i]) - p[i]== 1)  
    @constraint(model, [i in 2:n], sum(x[i, :]) + p[i]== 1) 

    @constraint(model, [i in 2:n], t[i] + sum(f[i,:]) == sum(f[:,i]) +I.demands[i])

    @constraint(model, [i in 2:n, j in 2:n], I.demands[i]*x[i,j] <= f[i,j])

    @constraint(model, [i in 2:n, j in 2:n], (I.capacity - I.demands[j])*x[i,j] >= f[i,j])

    @constraint(model, [i in 2:n], I.demands[i]*p[i] <= t[i])

    @constraint(model, [i in 2:n], t[i] <= I.capacity*p[i])

    #@constraint(model, sum(x_0) <= 2*K)
    #@constraint(model, sum(p) <= K)
    #@constraint(model, sum(x_0) >= 2*ceil(sum(I.demands)/I.capacity))
    #@constraint(model, sum(p) >= ceil(sum(I.demands)/I.capacity))

    #@constraint(model,[j in 2:n], sum(p[k] for k in j:n) >=ceil(sum(I.demands[k] for k in j:n)/I.capacity))

    @constraint(model, [i in 2:n-1], i <= u[i])
    @constraint(model, [i in 2:n-1], u[i] <= i*p[i] + (n-1)*(1-p[i]))

    
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
"""

function build_BPF_model(I::Instance)
    model = Model(CPLEX.Optimizer)
    Q = I.capacity
    demands = I.demands
    distances = I.distances
    # Ensemble des clients (excluant le dépôt)
    customers = 2:n
    
    # Variables de décision
    @variable(model, x[i in 1:n, j in 1:n], Bin)    # 1 si arc (i, j) utilisé
    @variable(model, p[i in customers], Bin)        # 1 si i est un client peak
    @variable(model, t[i in customers] >= 0, Int)   # Charge transportée pour chaque peak
    @variable(model, f[i in 1:n, j in 1:n] >= 0, Int)  # Flux transporté sur l’arc (i, j)
    @variable(model, u[i in customers], Int)        # Indice max visité sur le chemin

    # Contraintes d'équilibre de flux (chaque client est servi une seule fois)
    @constraint(model, [i in customers], sum(x[:, i]) + p[i] == 1)  
    @constraint(model, [i in customers], sum(x[i, :]) - p[i] == 1)

    # Contraintes de charge transportée
    @constraint(model, [i in customers], demands[i] * p[i] <= t[i])
    @constraint(model, [i in customers], t[i] <= Q * p[i])

    # Conservation du flux
    @constraint(model, [i in customers], t[i] + sum(f[i, :]) == sum(f[:, i]) + demands[i])

    # Limites de capacité sur les flux
    @constraint(model, [i in customers, j in customers], demands[i] * x[i, j] <= f[i, j])
    @constraint(model, [i in customers, j in customers], (Q - demands[j]) * x[i, j] >= f[i, j])

    # Contrainte de sous-tour MTZ
    @constraint(model, [i in customers], i <= u[i])
    @constraint(model, [i in customers], u[i] <= i * p[i] + (n - 1) * (1 - p[i]))

    @constraint(model, [i in customers, j in customers; i != j], 
                u[i] - u[j] + (n - j - 1) * x[i, j] <= n - j - 1)

    # Contrainte sur le nombre de véhicules utilisés
    #@constraint(model, sum(p) == K)
    #@constraint(model, sum(x[1, j] for j in customers) == 2 * K)

    # Objectif : minimiser la distance parcourue
    @objective(model, Min, sum(distances[i, j] * x[i, j] for i in 1:n, j in 1:n))

    return model
end

model = build_BPF_model(instance)

set_optimizer_attribute(model, "CPXPARAM_Conflict_Algorithm", 1)  # Active l'analyse de conflit
set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_LowerCutoff", -1e20)
set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_UpperCutoff", 1e20)
set_optimizer_attribute(model, "CPXPARAM_Conflict_Display", 2)  # Niveau de détail de l'analyse
set_optimizer_attribute(model, "CPXPARAM_MIP_Display", 4)
set_optimizer_attribute(model, "CPXPARAM_Conflict_Display", 2)

optimize!(model)
println("Nombre de variables x: ", length(all_variables(model)))



if termination_status(model) == MOI.OPTIMAL
    println("Optimal objective value: ", objective_value(model))
    x_sol = value.(model[:x])
    x0_sol = value.(model[:x_0])
    p_sol = value.(model[:p])
    t_sol = value.(model[:t])

    solution_x = [(i, j) for i in 2:n, j in 2:n if value(x_sol[i, j]) > 0.5]
    solution_x0 =  [(1,j) for j in 2:n if value(x0_sol[j]) > 0.5]
    println(solution_x, solution_x0)
    println([i for i in 2:n if value(p_sol[i]) < 0.5])
    println([value(t_sol[i]) for i in 2:n])

else
    println("⚠️ Le modèle est infaisable. Analyse des conflits en cours...")
    compute_conflict!(model)  # Active l'analyse d'infaisabilité
end

