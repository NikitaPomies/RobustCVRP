
using Random


function generate_distance_matrix(n; random_seed=2)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    prices = 100 * rand(rng, n)
    demands = 100 * rand(rng, n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1:n, j in 1:n]
    return X, Y, d, prices, demands
end

function update_distance_matrix!(d::Matrix{Float64})
    # Generate a random α_i for each row
    alpha = rand(0:80, size(d, 1))
    
    # Update the matrix
    for i in 1:size(d, 1), j in 1:size(d, 2)
        d[i, j] -= alpha[i]
    end

    return d
end

n = 10
X, Y, d, prices, demands = generate_distance_matrix(n)
capacity = 15

d = update_distance_matrix!(d)





mutable struct Label
    dummy_ressource::Vector{Int}
    ressources::Vector{Any}
    cost::Any
    current_node::Int
end



function Base.copy(label::Label)
    return Label(
        copy(label.dummy_ressource),
        copy(label.ressources),
        label.cost,
        label.current_node
    )
end

function extend(label::Label, next_node::Int, distances::Matrix{T}, demands::Vector{T}, capacity::Int) where {T <: Real}
    #Create new label
    next_node_consumption = label.ressources[end] + demands[next_node]
    if next_node_consumption > capacity
        return nothing
    end
    new_label = copy(label)
    new_label.dummy_ressource[next_node] += 1
    new_label.cost += distances[label.current_node, next_node]
    new_label.current_node = next_node
    push!(new_label.ressources, next_node_consumption)
    return new_label
end




function dominates(label1::Label, label2::Label)
    cost_dominance = label1.cost <= label2.cost
    ressource_dominance = label1.ressources[end] <= label2.ressources[end]
    #dummy_ressource_dominance = label1.dummy_ressource[end] <= label2.dummy_ressource
    return all([cost_dominance, ressource_dominance])
end

function add_non_dominated!(label_liste::Vector{Label}, new_label::Label)
    # Check if the new label is dominated by any existing label
    if any(dominates(label, new_label) for label in label_liste)
        return false,label_liste  # The new label is dominated; no changes to label_liste
    end

    # Remove labels that are dominated by the new label
    label_liste = filter(label -> !dominates(new_label, label), label_liste)

    # Add the new label to the list
    push!(label_liste, new_label)

    return true,label_liste
end


function mono_direction_dp(distances,demands,capacity)
    src = 1
    dest = n
    label_lists = [Vector{Label}() for _ in 1:n]
    push!(label_lists[src],Label(zeros(Int,n),Int[0],0,src))
    label_lists[src][1].dummy_ressource[1] = 1
    E = [src]
    print(label_lists)
    while length(E) > 0
        current_node = pop!(E)
        for label in label_lists[current_node]

            for j in 1:n
                if j!=current_node && label.dummy_ressource[j] == 0   #complete graph

                    #println(label)
                    new_label = extend(label,j,distances,demands,capacity)
                    if new_label === nothing
                        continue
                    else
                        changed,label_lists[j] = add_non_dominated!(label_lists[j],new_label)
                    end
                    if changed
                        push!(E,j)
                    end

                end
            end
        end
        #println(E)
        filter!(e->e!=current_node,E)


    end


return label_lists
end

d[1,10] = 33333
demands[10] = 0
test = mono_direction_dp(d,demands,200)
println(test[10])
#println(d[1,10])
#println(demands)