# Dualization

After dualizing the following maximization problem,


```math

\begin{align}
 \eta(x) = {max}_{\, \delta^1,\delta^2}~& \;\sum_{(i,j)\in A} \left[  \delta_{ij}^1 (\hat{t}_i + \hat{t}_j) + \delta_{ij}^2 \,  \hat{t}_i \, \hat{t}_j  \right ] \, x_{ij}  & \\
 {s.t. } &\sum_{(i,j)\in A} \delta_{ij}^1 \leq T  & \\
  &\sum_{(i,j)\in A} \delta_{ij}^2 \leq T^{\,2} & \\
\
& \delta_{ij}^1 \leq 1& \forall  (i,j)\in A\\
&\delta_{ij}^2 \leq 2 & \forall  (i,j)\in A\\
& \delta_{ij}^1, \delta_{ij}^2  \geq 0 & \forall (i,j)\in A
\end{align}
```

our robust problem can be formulated as follow : 



```math
\begin{align*}
 {min}_{\, x, \lambda, \mu}~& \;\sum_{(i,j)\in A} t_{ij} \,  x_{ij} + \lambda_1 T  + \lambda_2 T^{2} + \sum_{(i,j)\in A}   \mu_{ij}^1 +  2\sum_{(i,j)\in A}  \mu_{ij}^2& \\
 {s.t} &\sum_{i \in [n]} x_{ij} = 1& \forall j \in [2,n]\\
  &\sum_{i \in [n]} x_{ji} = 1& \forall j \in [2,n]\\
 &w_i - w_j \leq C \, (1 - x_{ij}) - d_j &\forall (i,j) \in A\\
%  & x_{ij} = x_{ji} &\forall (i,j) \in A \\
& d_i\leq w_i \leq C  & \forall i \in [n]\\ 
  &\mu_{ij}^1 + \lambda_1 \geq (\hat{t}_i + \hat{t}_j) \, x_{ij}  & \forall (i,j)\in A\\
 &\mu_{ij}^2 + \lambda_2 \geq \hat{t}_i \, \hat{t}_j \, x_{ij}  & \forall (i,j)\in A\\
& \mu_{ij}^1, \mu_{ij}^2  \geq 0 & \forall (i,j)\in A\\
& \lambda_1, \lambda_2  \geq 0 & \\
& x_{ij} \in \{0,1\}& \forall (i,j) \in A
\end{align*}
```

This formulation can now be solved using any milp solver.
