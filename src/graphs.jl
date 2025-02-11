using Graphs, GraphPlot, Compose

# Création d'un graphe orienté avec 5 sommets
g = DiGraph(5)  # Graphe dirigé (orienté)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)
add_edge!(g, 1, 4)
add_edge!(g, 1, 5)
add_edge!(g, 3, 4)
add_edge!(g, 2, 5)


labels = [string(i) for i in 1:nv(g)]
gplot(g, nodelabel=labels, layout=spring_layout)

