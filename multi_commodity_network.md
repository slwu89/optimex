# The network multi-commodity flow problem
Sean L. Wu
2023-12-19

The network multi-commondity flow problem is an extension of the [The
multi-commodity flow problem](@ref), where instead of having a bipartite
graph of supply and demand nodes, the graph can contains a set of nodes,
$i \in \mathcal{N}$ , which each have a (potentially zero) supply
capacity, $u^s_{i,p}$, and (potentially zero) a demand, $d_{i,p}$ for
each commodity $p \in P$. The nodes are connected by a set of edges
$(i, j) \in \mathcal{E}$, which have a shipment cost $c^x_{i,j,p}$ and a
total flow capacity of $u^x_{i,j}$.

Our take is to choose an optimal supply for each node $s_{i,p}$, as well
as the optimal transshipment $x_{i,j,p}$ that minimizes the total cost.

The mathematical formulation is:

``` math
\begin{aligned}
\min \;\; & \sum_{(i,j)\in\mathcal{E}, p \in P} c^x_{i,j,p} x_{i,j,p} + \sum_{i\in\mathcal{N}, p \in P} c^s_{i,p} s_{i,p} \\
s.t. \;\; & s_{i,p} + \sum_{(j, i) \in \mathcal{E}} x_{j,i,p} - \sum_{(i,j) \in \mathcal{E}} x_{i,j,p} = d_{i,p} & \forall i \in \mathcal{N}, p \in P \\
          & x_{i,j,p} \ge 0           & \forall (i, j) \in \mathcal{E}, p \in P \\
          & \sum_{p \in P} x_{i,j,p} \le u^x_{i,j}           & \forall (i, j) \in \mathcal{E} \\
          & 0 \le s_{i,p} \le u^s_{i,p} & \forall i \in \mathcal{N}, p \in P.
\end{aligned}
```
