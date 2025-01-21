

mutable struct Label{T <: Real}
    dummy_ressource::Vector{Int}
    ressources::Vector{T}
    cost::T
    current_node::Int
end

function Label{T}(dummy_ressource::Vector{Any}, ressources::Vector{Any}, cost::Any, current_node::Int) where T <: Real
    return Label{T}(
        convert(Vector{Int}, dummy_ressource),
        convert(Vector{T}, ressources),
        convert(T, cost),
        current_node
    )
end

function Base.copy(label::Label{T}) where T <: Real
    return Label{T}(
        copy(label.dummy_ressource),
        copy(label.ressources),
        label.cost,
        label.current_node
    )
end

function extend(label::Label,next_node::Int,distances::Matrix{Int},demands::Vector{Int},capacity::C)
    #Create new label
    next_node_consumption = label.ressources[end] + demands[next_node]
    if next_node_consumption > capacity
        return nothing
    end
    new_label = copy(label)
    new_label.dummy_ressource[next_node] +=1
    new_label.cost += distances[label.current_node,next_node]
    new_label.current_node = next_node
    push!(new_label.ressources,next_node_consumption)
    return new_label
end


function dominates(label1::Label, label2:: Label)
    cost_dominance = label1.cost <= label2.cost
    ressource_dominance = all(label1.ressources .<= label2.ressources)
    dummy_ressource_dominance = all(label1.dummy_ressource .<= label2.dummy_ressource)
    return all([cost_dominance,ressource_dominance,dummy_ressource_dominance])
end

