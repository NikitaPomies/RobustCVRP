using JuMP, CPLEX, MathOptInterface
const MOI = MathOptInterface

include("instance.jl")

instance = read_instance("../data/n_6-euclidean_true")

function build_BPF_model(I::Instance)
    # n : nombre total de nœuds (1 = dépôt, 2..n = clients)
    # Q : capacité des véhicules
    # demands : vecteur de demandes (demands[1] = 0 pour le dépôt)
    # distances : matrice n×n des coûts (distance) entre nœuds
    # K : nombre de véhicules (routes)
    Q = I.capacity
    demands = I.demands
    distances = I.distances
    
    model = Model(CPLEX.Optimizer)
    
    # Définition des ensembles
    depot = 1
    customers = 2:n           # ensemble des clients
    nonpeak_customers = customers[1:end-1]  # pour les contraintes MTZ (tous sauf le client de plus haut indice)
    
    # --- Variables ---
    # x0[j] : variable indiquant le nombre d'arcs utilisés entre le dépôt et le client j (j ∈ customers).
    # Dans le papier, x0[j] ∈ {0,1,2}. On peut déclarer une variable entière avec borne supérieure 2.
    @variable(model, x0[j in customers], Bin)  # On part du principe que x0 est binaire (souvent, il suffit de savoir si l'arc est utilisé).
    # Si nécessaire, on peut autoriser la valeur 2 en modifiant la borne supérieure :
    for j in customers
        set_upper_bound(x0[j], 2)
    end

    # x[i,j] : pour i,j ∈ customers, i ≠ j, indique si l'arc (i, j) est utilisé.
    @variable(model, x[i in customers, j in customers], Bin)
    for i in customers, j in customers
        if i == j
            fix(x[i,j], 0; force = true)
        end
    end

    # p[i] : 1 si le client i est le peak de sa route.
    @variable(model, p[i in customers], Bin)
    # Selon le papier, le client de plus haut indice est forcement un peak :
    fix(p[n], 1; force = true)
    
    # t[i] : charge (total demand) sur la route se terminant en le client i (si i est peak).
    @variable(model, t[i in customers] >= 0)

    # f[i,j] : flux (quantité livrée) sur l'arc (i,j) entre clients.
    @variable(model, f[i in customers, j in customers] >= 0)
    
    # u[i] : variable pour éliminer les sous-tours (pour i dans customers sauf le plus grand indice)
    @variable(model, u[i in nonpeak_customers], Int)

    # --- Contraintes ---
    # (10) Contrainte de départ : pour chaque client i, la somme des arcs partant de i (vers d'autres clients)
    # plus l'indicateur p[i] vaut 1.
    @constraint(model, [i in customers], sum(x[i, j] for j in customers if j != i) + p[i] == 1)

    # (11) Contrainte d'arrivée : pour chaque client i, la somme des arcs entrant dans i
    # moins l'indicateur p[i] vaut 1.
    @constraint(model, [i in customers], sum(x[j, i] for j in customers if j != i) - p[i] == 1)
    
    # (15) Contrainte sur le dépôt : la somme des arcs du dépôt vers les clients vaut 2K.
    #@constraint(model, sum(x0[j] for j in customers) == 2*K)
    
    # (16) Contrainte sur le nombre de peaks : il doit y avoir exactement K peaks.
    #@constraint(model, sum(p[i] for i in customers) == K)
    
    # (12) Conservation du flux pour chaque client i.
    # La charge t[i] plus le flux sortant de i doit être égale au flux entrant en i plus la demande du client i.
    @constraint(model, [i in customers], t[i] + sum(f[i, j] for j in customers) == sum(f[j, i] for j in customers) + demands[i])
    
    # (13) Pour chaque arc entre clients, si l'arc (i,j) est utilisé, alors
    # f[i,j] est compris entre demands[i] et Q - demands[j].

    @constraint(model, [i in customers, j in customers; i != j], f[i,j] >= demands[i] * x[i,j])
    @constraint(model, [i in customers, j in customers; i != j], f[i,j] <= (Q - demands[j]) * x[i,j])
    
    # (14) Pour chaque client, si ce client est peak, alors t[i] est entre demands[i] et Q.
    @constraint(model, [i in customers], t[i] >= demands[i] * p[i])
    @constraint(model, [i in customers], t[i] <= Q * p[i])
    
    # Contraintes MTZ (pour éliminer les sous-tours) sur les variables u.
    @constraint(model, [i in nonpeak_customers], u[i] >= i)  # borne inférieure
    @constraint(model, [i in nonpeak_customers], u[i] <= i * p[i] + (n - 1)*(1 - p[i]))
    @constraint(model, [i in nonpeak_customers, j in nonpeak_customers; i != j],
        u[i] - u[j] + (n - j - 1)* x[i,j] <= n - j - 1)
    
    # --- Objectif ---
    # Coût = coût des arcs du dépôt vers les clients + coût des arcs entre clients.
    @objective(model, Min, 
        sum(distances[depot, j] * x0[j] for j in customers) + 
        sum(distances[i,j] * x[i,j] for i in customers, j in customers if i != j))
    
    return model
end


model = build_BPF_model(instance)

println("Optimisation en cours ...")
optimize!(model)

if termination_status(model) == MOI.OPTIMAL
    println("Valeur optimale de l'objectif: ", objective_value(model))
    println("\nArcs départ du dépôt :")
    for j in 2:n
        if value(model[:x0][j]) > 0.5
            println("Dépôt -> $j : ", value(model[:x0][j]))
        end
    end
    println("\nArcs entre clients :")
    for i in 2:n, j in 2:n
        if i != j && value(model[:x][i,j]) > 0.5
            println("$i -> $j")
        end
    end
    println("\nClients identifiés comme peaks : ", [i for i in 2:n if value(model[:p][i]) > 0.5])
else
    println("⚠️ Le modèle est infaisable.")
end
