using LinearAlgebra

# Function to parse the CVRP instance file
function parse_cvrp_instance(file_path)
    # Initialize placeholders
    name = ""
    comment = ""
    capacity = 0
    dimension = 0
    coords = Dict()
    demands = Dict()
    depot = 0

    # Read the file
    open(file_path, "r") do file
        section = ""
        for line in eachline(file)
            # Remove leading/trailing whitespace
            line = strip(line)
            
            # Skip empty lines
            if isempty(line)
                continue
            end

            # Parse the sections
            if startswith(line, "NAME")
                name = split(line, ":")[2] |> strip
            elseif startswith(line, "COMMENT")
                comment = split(line, ":")[2] |> strip
            elseif startswith(line, "CAPACITY")
                capacity = parse(Int, split(line, ":")[2] |> strip)
            elseif startswith(line, "DIMENSION")
                dimension = parse(Int, split(line, ":")[2] |> strip)
            elseif startswith(line, "NODE_COORD_SECTION")
                section = "NODE_COORD_SECTION"
            elseif startswith(line, "DEMAND_SECTION")
                section = "DEMAND_SECTION"
            elseif startswith(line, "DEPOT_SECTION")
                section = "DEPOT_SECTION"
            elseif startswith(line, "EOF")
                break
            elseif section == "NODE_COORD_SECTION"
                parts = split(line)
                node = parse(Int, parts[1])
                coords[node] = (parse(Float64, parts[2]), parse(Float64, parts[3]))
            elseif section == "DEMAND_SECTION"
                parts = split(line)
                node = parse(Int, parts[1])
                demands[node] = parse(Int, parts[2])
            elseif section == "DEPOT_SECTION"
                if line == "-1"
                    break
                else
                    depot = parse(Int, line)
                end
            end
        end
    end

    # Compute the distance matrix
    distance_matrix = zeros(dimension, dimension)
    for i in 1:dimension
        for j in 1:dimension
            if i != j
                xi, yi = coords[i]
                xj, yj = coords[j]
                distance_matrix[i, j] = round(Int,sqrt((xi - xj)^2 + (yi - yj)^2))
            end
        end
    end

    return name, comment, capacity, dimension, coords, demands, depot, distance_matrix
end

