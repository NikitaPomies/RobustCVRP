# Cutting Planes

Our robust formulation is equivalent to :
$$
\begin{align*}
 {min}_{\, x, Z} ~& \; {Z } \\
{s.t} & \sum_{(i,j)\in A} t'_{ij} \, x_{ij} \leq Z & \forall \, t' \in \mathcal{U} \\
& x \in \mathcal{X} 
\end{align*}
$$

where $\mathcal{X}$ is the set of admissible CVRP solutions.

Our cutting planes algorithm works by starting with a restricted set  $\, \mathcal{U'}$ of scenarios in $\mathcal{U}$. In the B&B scheme, new scenarios are added each time a integer solution is found by solving the following subproblem : 

$$
\begin{align*}
{max}_{\, \delta^1,\delta^2} ~& \;\sum_{(i,j)\in A} t_{ij} \, x^*_{ij} + (\delta_{ij}^1(\hat{t}_i + \hat{t}_j) + \delta_{ij}^2 \hat{t_i}\hat{t_j}) \, x^*_{ij} \\
 {s.t} & \sum_{(i,j)\in A} \delta_{ij}^1 \leq T \\
& \sum_{(i,j)\in A} \delta_{ij}^2 \leq T^2 \\
& \delta_{ij}^1 \in [0,1], \delta_{ij}^2 \in [0,2] & \forall (i,j) \in A 
\end{align*}
$$

This subproblem is solved by sorting the selected edges of the current node solution and then added via  LazyCallbacks.