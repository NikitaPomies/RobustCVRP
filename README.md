# Robust Capacitated Vehicle Routing Problem

An implementation of various methods to solve a RCVRP with travel times uncertainty. 

We tested four different methods : 
- [Cutting planes](https://github.com/NikitaPomies/RobustCVRP/tree/main/src/cutting_planes)
- [Dualisation](https://github.com/NikitaPomies/RobustCVRP/tree/main/src/dualisation)
- A Hybrid Genetic Search meta heuristic
- [Column generation](https://github.com/NikitaPomies/RobustCVRP/tree/main/src/column_generation).

For symmetric instances, we also tried a formulation from litterature that break symmetry, called [Peak Formulation](https://github.com/NikitaPomies/RobustCVRP/blob/main/src/peak_formulation).

For development convenience, the metaheuristic code is in another repo : https://github.com/NikitaPomies/Robust-HGS-CVRP.

A technical report (in French) presenting the results is available here.

##  Description of the problem 
The uncertainty set it described as follow : 


```math
\mathcal{U} = \{ t'_{ij} = t_{ij} +  \delta_{ij}^1 (\hat{t}_i + \hat{t}_j) + \delta_{ij}^2 \,  \hat{t}_i \, \hat{t}_j  \: \text{ s.t} 
```

```math
\sum_{(i,j)\in A} \delta_{ij}^1 \leq T  , \sum_{(i,j)\in A} \delta_{ij}^2 \leq T^{\,2}, \delta_{ij}^1 \in [0,1], \;  \delta_{ij}^2 \in [0,2] \;\forall (i,j) \in A \}
```
<br>


T and $\hat{t}$ are parameters of the problem 

```math
\begin{align*}
 {min}_{\,x}~& \; {max}_{\, \delta^1,\delta^2} \sum_{(i,j)\in A} t_{ij} \, x_{ij} + (\delta_{ij}^1(\hat{t}_i + \hat{t}_j) + \delta_{ij}^2 \hat{t_i}\hat{t_j}) \, x_{ij} & \\
 {s.t} &\sum_{i \in [n]} x_{ij} = 1& \forall j \in [2,n]\\
  &\sum_{i \in [n]} x_{ji} = 1& \forall j \in [2,n]\\
 &w_i - w_j \leq C \, (1 - x_{ij}) - d_j &\forall (i,j) \in A\\
& d_i\leq w_i \leq C  & \forall i \in [n]\\ 
& \sum_{(i,j)\in A} \delta_{ij}^1 \leq T \\
& \sum_{(i,j)\in A} \delta_{ij}^2 \leq T^2 \\
& x_{ij} \in \{0,1\}& \forall (i,j) \in A\\
& \delta_{ij}^1 \in [0,1], \delta_{ij}^2 \in [0,2] & \forall (i,j) \in A 
\end{align*}
```
