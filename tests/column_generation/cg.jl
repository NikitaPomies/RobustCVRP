using JuMP, CPLEX, LinearAlgebra

file_path = "../data/P-n20-k2.vrp"  # Replace with your instance file path

name, comment, capacity, n, coords, demands, depot, distances = parse_cvrp_instance(file_path)

demands = [demands[j] for j in 1:n]


##Construct initial routes

#routes = nearest_neighbor_routes(distances, demands, capacity)

routes = init_column_pool(distances, demands, capacity, 200)
#println(routes)

#= test = 0
for route in routes
    cost = compute_route_cost(route, distances)
    global test += cost

end
println(test)

test = 0

for route in routes
    println("Changing new route")
    println(route)
    cost = compute_route_cost(route,distances)
    println(cost)
    new_route = lkh_two_opt(route,distances)
    println(new_route)
    new_cost = compute_route_cost(new_route,distances)
    println(new_cost)
    

end =#
 #Define restricted master problem

rmasterpb = Model(CPLEX.Optimizer)

R = length(routes)
@variable(rmasterpb, x[1:R]>=0)


@constraint(rmasterpb, c[i = 2:n], sum(x[r] for r in 1:R if i in routes[r]) >= 1)


@constraint(rmasterpb, con, sum(x) <= 5) # nombre max de véhicules

@objective(rmasterpb, Min, sum(compute_route_cost(routes[r], distances) * x[r] for r in 1:R))

#Define the sub_model

submodel = build_tsp_model(demands, capacity, n)
set_attribute(
    submodel,
    MOI.LazyConstraintCallback(),
    cb_data -> subtour_elimination_callback(cb_data, submodel)
)

MAXIMUM_ITERATIONS = 100
#set_silent(model)
optimize!(rmasterpb)

@assert is_solved_and_feasible(rmasterpb)
#println("Routes : ",routes)
println(value.(rmasterpb[:x]))


for k in 1:200
    println("Iteration $(k)")
    #println(routes)
    #set_silent(rmasterpb)

    @assert is_solved_and_feasible(rmasterpb;dual=true)
    lower_bound = objective_value(rmasterpb)
    #println(value.(rmasterpb[:x]))
    #println(con)
    #println(c)
    #println(objective_function(rmasterpb))
    println("obj value : ", lower_bound)

    local prices = [0.0 for i in 1:n]
    for i in 2:n
        prices[i] = dual(c[i])
    end
    #println(prices)
    set_silent(submodel)
    solvepctsp(prices, submodel,distances)
    println("valeur pctsp ",objective_value(submodel))
    
    #println(prices)
    route = get_route(value.(submodel[:x]), value.(submodel[:z]))
    if route in routes
        println("duplicate_routes")
        #break
        continue
    end
    route_cost = compute_route_cost(route,distances)
    println(route, route_cost)
    push!(routes, route)
    push!(x, @variable(rmasterpb, lower_bound = 0))
    set_objective_coefficient(rmasterpb, x[end], route_cost)
    for i in 2:n
        if i in route
            # Get the constraint for customer i
            #customer_con = c[i] # i-1 because we start from customer 2
            # Add variable to constraint with coefficient 1
            set_normalized_coefficient(c[i], x[end], 1.0)
        else
            set_normalized_coefficient(c[i],x[end],0)
        end
    end
    set_normalized_coefficient(con, x[end], 1.0)
    optimize!(rmasterpb)
    #println(value.(rmasterpb[:x]))



end

sol = value.(rmasterpb[:x])
for i in 1:length(routes)-1
    if sol[i] > 0
        println(routes[i])
    end
end


for i in 1:length(routes) - 1
    set_integer(rmasterpb[:x][i])
end
optimize!(rmasterpb)



sol = value.(rmasterpb[:x])
for i in 1:length(routes)-1
    if sol[i] > 0
        println(routes[i])
    end
end

println(objective_value(rmasterpb))



 
