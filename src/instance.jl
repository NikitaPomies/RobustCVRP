Base.@kwdef struct Instance
    distances::Matrix{Int64}
    demands::Vector{Int64}
    th::Vector{Int64}
    T::Int64
    capacity::Int64
end


function read_instance(file_path::String)::Instance
    
    include(file_path)
    
    return Instance(
        distances = t,
        demands = d,
        th = th,
        T = T,
        capacity = C
    )
end
