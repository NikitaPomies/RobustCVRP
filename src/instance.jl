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
        distances=t,
        demands=d,
        th=th,
        T=T,
        capacity=C
    )
end

Base.@kwdef struct Solution
    routes::Vector{Vector{Int64}}
    cost::Int64
end

function read_solution(file_path::String)::Solution
    routes = Vector{Vector{Int64}}()
    cost = 0

    open(file_path, "r") do file
        for line in eachline(file)
            if startswith(line, "Route")
                # Extract the route as a vector of integers
                route = parse.(Int, split(strip(split(line, ":")[2])))
                push!(routes, route)
            elseif startswith(line, "Cost")
                # Extract the cost
                cost = parse(Int, split(line)[2])
            end
        end
    end

    return Solution(routes=routes, cost=cost)
end

function compute_robust_cost(I::Instance, solution::Solution)::Int64

    n = length(I.demands)

    selected_e = []
    for route in solution.routes
        push!(selected_e, (1, route[1] + 1))
        for i in 1:length(route)-1
            push!(selected_e, (route[i] + 1, route[i+1] + 1))
        end
        push!(selected_e, (route[end] + 1, 1))
    end

    sort!(selected_e, by=x -> (I.th[x[1]] + I.th[x[2]]), rev=true)




    S1 = 0
    robust_cost_1 = 0
    for (i, j) in selected_e
        if S1 + 1 <= I.T
            # Can take full item
            robust_cost_1 += I.th[i] + I.th[j]
            S1 += 1.0
        else
            # Take fractional part
            remaining = I.T - S1
            if remaining > 0
                robust_cost_1 += remaining * (I.th[i] + I.th[j])
                S1 += remaining
            end
            break  # We've reached I.T, no need to continue
        end
    end

    sort!(selected_e, by=x -> (I.th[x[1]] * I.th[x[2]]), rev=true)


    S2 = 0
    robust_cost_2 = 0
    for (i, j) in selected_e
        if S2 + 2 <= I.T * I.T
            # Can take full item
            robust_cost_2 += 2 * (I.th[i] * I.th[j])
            S2 += 2.0
        else
            # Take fractional part
            remaining = I.T * I.T - S2
            if remaining > 0
                robust_cost_2 += remaining * (I.th[i] * I.th[j])
            end
            break  # We've reached T, no need to continue
        end
    end
    obj = sum(I.distances[i, j] for (i, j) in selected_e) + robust_cost_1 + robust_cost_2
    return obj

end
