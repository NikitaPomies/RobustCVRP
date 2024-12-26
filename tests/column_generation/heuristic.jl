using Random

function nearest_neighbor_routes(
    distances::Matrix{T1},
    demands::Vector{T2},
    capacity::T3
) where {T1<:Real,T2<:Real,T3<:Real}
    n_customers = length(demands)  # Excluding depot
    unvisited = Set(2:n_customers)  # Customer indices (depot is 1)
    routes = Vector{Vector{Int}}()

    while !isempty(unvisited)
        current_route = [1]  # Start at depot
        current_load = 0.0

        while !isempty(unvisited)
            current = current_route[end]
            # Find nearest unvisited customer that doesn't exceed capacity
            nearest = nothing
            min_dist = Inf

            for next_customer in unvisited
                if current_load + demands[next_customer] <= capacity
                    dist = distances[current, next_customer]  # +1 for Julia indexing
                    if dist < min_dist
                        min_dist = dist
                        nearest = next_customer
                    end
                end
            end

            # If no feasible customer found, close the route
            if isnothing(nearest)
                break
            end

            # Add customer to route
            push!(current_route, nearest)
            current_load += demands[nearest]
            delete!(unvisited, nearest)
        end

        # Return to depot
        push!(current_route, 1)
        #=         println( current_route)
                indices = 2:length(current_route)-1  # Indices of the middle elements
                middle_elements = current_route[indices]  # Extract the middle elements
                shuffle!(middle_elements)  # Shuffle the middle elements
                current_route[indices] .= middle_elements 
                println(current_route) =#
        push!(routes, current_route)
    end

    # Generate random routes
    #routes = Vector{Vector{Int}}()
    num_random_routes = 500
    for _ in 1:num_random_routes
        unvisited = Set(2:n_customers)  # Reset unvisited customers

        while !isempty(unvisited)
            current_route = [1]  # Start at depot
            current_load = 0.0

            while !isempty(unvisited)
                # Random route generation
                feasible_customers = [c for c in unvisited if current_load + demands[c] <= capacity]

                if isempty(feasible_customers)
                    break
                end

                # Select a random customer from the feasible list
                random_customer = rand(feasible_customers)

                # Add customer to route
                push!(current_route, random_customer)
                current_load += demands[random_customer]
                delete!(unvisited, random_customer)
            end

            # Return to depot
            push!(current_route, 1)
            push!(routes, current_route)
        end
    end



    return routes
end


function compute_route_cost(route::Vector{T1}, distances::Matrix{Float64}) where {T1<:Real}
    S = 0
    for i in 1:length(route)-1
        u, v = route[i], route[i+1]
        S += distances[u, v]
    end
    return S
end

function calculate_tour_length(tour, distances)
    n = length(tour)
    return sum(distances[tour[i], tour[i+1]] for i in 1:n-1)
end

function two_opt_swap!(tour, i, j)
    reverse!(@view tour[i:j])
    return tour
end

function lkh_two_opt(route, distances)
    current_tour = copy(route)
    improved = true
    max_iterations = 100
    iterations = 0

    while improved && iterations < max_iterations
        improved = false
        current_length = calculate_tour_length(current_tour, distances)

        for i in 2:length(current_tour)-1
            for j in i+1:length(current_tour)-1
                new_tour = copy(current_tour)
                two_opt_swap!(new_tour, i, j)
                new_length = calculate_tour_length(new_tour, distances)

                if new_length < current_length
                    current_tour = new_tour
                    current_length = new_length
                    improved = true
                    break
                end
            end
            improved && break
        end
        iterations += 1
    end

    return current_tour
end


function reverse_segment!(tour, i, j)
    reverse!(@view tour[i:j])
    return tour
end



function three_opt_swap!(tour, i, j, k, case)
    n = length(tour)
    if case == 1
        # Reverse segment i->j
        reverse_segment!(tour, i, j-1)
    elseif case == 2
        # Reverse segment j->k
        reverse_segment!(tour, j, k-1)
    elseif case == 3
        # Reverse both segments
        reverse_segment!(tour, i, j-1)
        reverse_segment!(tour, j, k-1)
    elseif case == 4
        # Rearrange segments without reversing
        temp_tour = copy(tour)
        # Copy first segment
        idx = 1
        for x in 1:i-1
            tour[idx] = temp_tour[x]
            idx += 1
        end
        # Copy second segment to third position
        for x in j:k-1
            tour[idx] = temp_tour[x]
            idx += 1
        end
        # Copy third segment to second position
        for x in i:j-1
            tour[idx] = temp_tour[x]
            idx += 1
        end
        # Copy rest of the tour
        for x in k:n
            tour[idx] = temp_tour[x]
            idx += 1
        end
    end
end

function lkh_2opt_3opt(route, distances)
    current_tour = copy(route)
    improved = true
    max_iterations = 100
    iterations = 0
    
    while improved && iterations < max_iterations
        improved = false
        current_length = calculate_tour_length(current_tour, distances)
        
        # Try 2-opt moves first
        for i in 2:length(current_tour)-1
            for j in i+1:length(current_tour)-1
                new_tour = copy(current_tour)
                two_opt_swap!(new_tour, i, j)
                new_length = calculate_tour_length(new_tour, distances)
                
                if new_length < current_length
                    current_tour = new_tour
                    current_length = new_length
                    improved = true
                    @goto next_iteration
                end
            end
        end
        
        # Try 3-opt moves if no 2-opt improvement was found
        for i in 2:length(current_tour)-2
            for j in i+1:length(current_tour)-1
                for k in j+1:length(current_tour)
                    # Try all possible 3-opt configurations
                    for case in 1:4
                        new_tour = copy(current_tour)
                        three_opt_swap!(new_tour, i, j, k, case)
                        new_length = calculate_tour_length(new_tour, distances)
                        
                        if new_length < current_length
                            current_tour = new_tour
                            current_length = new_length
                            improved = true
                            @goto next_iteration
                        end
                    end
                end
            end
        end
        
        @label next_iteration
        iterations += 1
    end
    
    return current_tour
end



function init_column_pool(
    distances::Matrix{T1},
    demands::Vector{T2},
    capacity::T3,
    n::Int
) where {T1<:Real,T2<:Real,T3<:Real}
    n_customers = length(demands)
    unvisited = [i for i in 2:n_customers]
    routes = Vector{Vector{Int}}()


    for i in 1:n
        shuffle!(unvisited)

        num_visited_clients = 0
        customer_idx = 1

        while num_visited_clients < n_customers - 1
            current_route = [1]  # Start at depot
            current_load = 0.0

            while (customer_idx <= length(unvisited)) && (current_load + demands[unvisited[customer_idx]] <= capacity)
                next_customer = unvisited[customer_idx]
                push!(current_route, next_customer)
                current_load += demands[next_customer]
                customer_idx += 1
                num_visited_clients += 1
            end
            push!(current_route, 1)

            # optimize the route here 
            current_route = lkh_2opt_3opt(current_route, distances)
            push!(routes, current_route)

        end
    end



    return routes
end


