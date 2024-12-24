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
        u,v = route[i], route[i + 1]
        S += distances[u, v]
    end
    return S
end