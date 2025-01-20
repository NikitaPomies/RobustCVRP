using JuMP, CPLEX,Gurobi, LinearAlgebra

include("../instance.jl")
include("heuristic.jl")
include("subproblem.jl")

#name, comment, capacity, n, coords, demands, depot, distances = parse_cvrp_instance(file_path)


instance  = read_instance("../../data/n_40-euclidean_true")


println(instance.distances)
println(instance.demands)

##Construct initial routes

#routes = nearest_neighbor_routes(distances, demands, capacity)

routes = init_column_pool(instance.distances, instance.demands, instance.capacity, 100)

function check_consecutive_pair(arr::Vector{Int}, pair::Tuple{Int,Int})
    i, j = pair
    n = length(arr)
    
    # Check each consecutive pair in the array
    for idx in 1:(n-1)
        if arr[idx] == i && arr[idx + 1] == j
            return true
        end
    end
    
    return false
end
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

    new_route = lkh_2opt_3opt(route,distances)
    println(new_route)
    new_cost = compute_route_cost(new_route,distances)
    println(new_cost)

    

end  =#
 #Define restricted master problem

rmasterpb = Model(Gurobi.Optimizer)
#set_optimizer_attribute(rmasterpb, "CPX_PARAM_LPMETHOD", 4)
#set_optimizer_attribute(rmasterpb, "CPX_PARAM_SOLUTIONTYPE", 2)
set_silent(rmasterpb)

R = length(routes)
@variable(rmasterpb, x[1:R]>=0)
@variable(rmasterpb, 0 <= mu_1[i=1:n, j=1:n])
@variable(rmasterpb, 0 <= mu_2[i=1:n, j=1:n])
@variable(rmasterpb, lambda_1 >= 0)
@variable(rmasterpb, lambda_2 >= 0)


@constraint(rmasterpb, c[i = 2:n], sum(x[r] for r in 1:R if i in routes[r]) >= 1)
@constraint(rmasterpb,z1[i=1:n,j=1:n],mu_1[i,j] + lambda_1 >= (th[i] + th[j]) * sum(x[r] for r in 1:R if check_consecutive_pair(routes[r],(i,j))))
println(z1)



@constraint(rmasterpb,z2[i=1:n,j=1:n],mu_2[i,j] + lambda_2 >= (th[i] * th[j]) * sum(x[r] for r in 1:R if check_consecutive_pair(routes[r],(i,j))))

@constraint(rmasterpb, con, sum(x) <=15) # nombre max de véhicules

@objective(rmasterpb, Min, sum(compute_route_cost(routes[r], instance.distances) * x[r] for r in 1:R) + lambda_1 * instance.T + lambda_2 * instance.T^2 + sum(mu_1 + 2 * mu_2))

#Define the sub_model

submodel = build_tsp_model(instance.demands, instance.capacity, n)
set_silent(submodel)
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
    new_distances = Float64.(deepcopy(instance.distances))
    for i in 1:n
        for j in 1:n
            if i!=j
                dual1 = dual(z1[i,j])
                dual2 = dual(z2[i,j])
                #println(dual1)
                new_distances[i,j]  +=  dual1*(th[i]+th[j])  + dual2*(th[i]*th[j]) 
                #println(new_distances[i,j] - instance.distances[i,j])
            end
        end
    end
    println("startin to solve subpb")
    solvepctsp(prices, submodel,new_distances)
    println("valeur pctsp ",objective_value(submodel))
    
    #println(prices)
    route = get_route(value.(submodel[:x]), value.(submodel[:z]))
    if route in routes
        println("duplicate_routes")
        break
        continue
    end
    route_cost = compute_route_cost(route,instance.distances)
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

    for i in 1:n
        for j in 1:n
            if i!=j
                if check_consecutive_pair(route,(i,j))
                    set_normalized_coefficient(z1[i,j],x[end],-(th[i]+th[j]))
                    set_normalized_coefficient(z2[i,j],x[end],(-th[i]*th[j]))
                else
                    set_normalized_coefficient(z1[i,j],x[end],0)
                    set_normalized_coefficient(z2[i,j],x[end],0)

                end
                
            end 
        end
    end
    optimize!(rmasterpb)
    #println(value.(rmasterpb[:x]))



end

sol = value.(rmasterpb[:x])
for i in 1:length(routes)
    if sol[i] > 0
        println(routes[i])
    end
end


println("value of relaxation : ", objective_value(rmasterpb))
for i in 1:length(routes)
    set_integer(rmasterpb[:x][i])
end
optimize!(rmasterpb)


println("heuristic_solution")
global S = 0
sol = value.(rmasterpb[:x])
for i in 1:length(routes)
    if sol[i] > 0
        println(routes[i],compute_route_cost(routes[i],instance.distances))
        global S+= compute_route_cost(routes[i],instance.distances)
    end
end

println(objective_value(rmasterpb))
println(S)



 
 