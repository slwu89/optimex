using Catlab, DataFrames
using JuMP, HiGHS
const VarType = Union{JuMP.VariableRef,Float64}

# the data types
@present SchProjGraph <: SchLabeledGraph begin
  Duration::AttrType
  duration::Attr(V,Duration)
  es::Attr(V,Duration) # earliest possible starting times
  ef::Attr(V,Duration) # earliest possible finishing times
  ls::Attr(V,Duration) # latest possible starting times
  lf::Attr(V,Duration) # latest possible finishing times
  float::Attr(V,Duration) # "float" time of task
end

to_graphviz(SchProjGraph)

@abstract_acset_type AbstractProjGraph <: AbstractGraph

@acset_type ProjGraph(SchProjGraph, index=[:src,:tgt]) <: AbstractProjGraph

"""
    Make a `ProjGraph` acset from a `DataFrame` with columns for
`Activity`, `Predecessor`, and `Duration`.
"""
function make_ProjGraph(input::DataFrame)
    g = @acset ProjGraph{Symbol,Int} begin
        V = size(input,1)
        label = input.Activity
        duration = input.Duration
    end

    for r in eachrow(input)
        if length(r.Predecessor) > 0
            prednodes = vcat(incident(g, r.Predecessor, :label)...)
            node = only(incident(g, r.Activity, :label))
            add_edges!(g, prednodes, fill(node, length(r.Predecessor)))
        else
            continue
        end
    end

    return g
end

# da codez

"""
    Perform a forward pass of the critical path sweep. This identifies
the earliest starting times `es` and finishing times `ef` for each node.
This assumes that the node with part id `1` is the start, and that the
node with highest part id is the end, and the graph is a DAG. Returns
the topological sort of the nodes.
"""
function forward_pass!(g::T) where {T<:AbstractProjGraph}
    g[1, :es] = 0
    g[1, :ef] = 0

    toposort = topological_sort(g)

    for v in toposort[2:end]
        pred = collect(inneighbors(g, v))
        g[v, :es] = maximum(g[pred, :ef])
        g[v, :ef] = g[v, :es] + g[v, :duration]
    end

    return toposort
end

"""
    Perform a backward pass of the critical path sweep. This identifies
the latest starting times `ls` and finishing times `lf` for each node.
It also computes the `float`, the maximum acceptable delay for each node.
This assumes that the graph is a DAG.
"""
function backward_pass!(g::T, toposort) where {T<:AbstractProjGraph}
    reverse!(toposort)

    g[toposort[1], :lf] = g[toposort[1], :es]
    g[toposort[1], :ls] = g[toposort[1], :es]
    
    for v in toposort[2:end]
        succ = collect(outneighbors(g, v))
        g[v, :lf] = minimum(g[succ, :ls])
        g[v, :ls] = g[v, :lf] - g[v, :duration]
    end

    # float: delay acceptable for each task w/out delaying completion of project
    g[:,:float] = g[:,:lf] - g[:,:ef]
end

"""
    Find a critical path, assumes that the forward and backward passes have been
completed.
"""
function find_critical_path(g::T) where {T<:AbstractProjGraph}
    # set of critical activities 
    ca = incident(g, 0, :float)

    # stack of nodes for dfs search
    S = [1]
    seen = zeros(Bool, nv(g))
    # edges on the critical path
    ce = Int[]

    while !isempty(S)
        v = pop!(S)
        if !seen[v]
            seen[v] = true
            # only iterate over neighbors who are critical actvities
            for w in intersect(outneighbors(g, v), ca)
                if g[v, :lf] == g[w, :es]
                    push!(S, w)
                    push!(ce, only(edges(g, v, w)))
                end
            end
        end
    end

    return findall(seen), ce
end

# ex: table 7.1
proj_df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:E,:F,:G,:H,:I,:J,:end],
    Predecessor = [
        [], [:start], [:A], [:A,:B], [:B], [:B,:C], [:C,:D,:E],
        [:D], [:F,:G], [:F,:G], [:I], [:H,:J]
    ],
    Duration = [0,5,3,7,4,6,4,2,9,6,2,0]
)

projnet = make_ProjGraph(proj_df)
to_graphviz(projnet, node_labels=:label)

toposort = forward_pass!(projnet)
backward_pass!(projnet, toposort)
cV, cE = find_critical_path(projnet)

cg = Subobject(projnet, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# ex: fig 7.4
proj_df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:end],
    Predecessor = [
        [], [:start], [:start], [:A,:B], [:A,:B], [:C,:D]
    ],
    Duration = [0,5,4,7,8,0]
)

projnet = make_ProjGraph(proj_df)
to_graphviz(projnet, node_labels=:label)

toposort = forward_pass!(projnet)
backward_pass!(projnet, toposort)
cV, cE = find_critical_path(projnet)

cg = Subobject(projnet, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# "reduce" time for D to 7
proj_df[proj_df.Activity .== :D, :Duration] .= 7
projnet = make_ProjGraph(proj_df)

toposort = forward_pass!(projnet)
backward_pass!(projnet, toposort)
cV, cE = find_critical_path(projnet)

cg = Subobject(projnet, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# CPM with acceleration as a LP problem
@present SchAccelProjGraph <: SchLabeledGraph begin
    VarType::AttrType
    TimeType::AttrType
    CostType::AttrType
    x::Attr(V,VarType) # decision var: duration
    y::Attr(V,VarType) # decision var: start time
    tmax::Attr(V,TimeType) # normal duration
    tmin::Attr(V,TimeType) # minimum duration
    Δc::Attr(V,CostType) # cost to reduce duration by one unit
end

to_graphviz(SchAccelProjGraph)

@abstract_acset_type AbstractAccelProjGraph <: AbstractGraph

@acset_type AccelProjGraph(SchAccelProjGraph, index=[:src,:tgt]) <: AbstractAccelProjGraph

"""
    Make a `AccelProjGraph` acset from a `DataFrame` with columns for
`Activity`, `Predecessor`, `Max`, `Min`, and `Cost`.
"""
function make_AccelProjGraph(input::DataFrame)
    g = @acset AccelProjGraph{Symbol,VarType,Int,Float64} begin
        V = size(input,1)
        label = input.Activity
        tmax = input.Max
        tmin = input.Min
        Δc = input.Cost
    end

    for r in eachrow(input)
        if length(r.Predecessor) > 0
            prednodes = vcat(incident(g, r.Predecessor, :label)...)
            node = only(incident(g, r.Activity, :label))
            add_edges!(g, prednodes, fill(node, length(r.Predecessor)))
        else
            continue
        end
    end

    return g
end

# ex: Fig III.12
proj_df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:E,:end],
    Predecessor = [
        [], [:start], [:start], [:A], [:A], [:B,:C], [:D,:E]
    ],
    Max = [0,3,5,4,4,5,0],
    Min = [0,3,2,1,1,2,0],
    Cost = [0,0,200,200,100,600,0]
)

projnet = make_AccelProjGraph(proj_df)
to_graphviz(projnet, node_labels=:label)

T = 11 # generalize somehow

# make the LP model

jumpmod = JuMP.Model(HiGHS.Optimizer)

# the decision vars
@variable(
    jumpmod, 
    projnet[v,:tmin] ≤ x_j[v ∈ vertices(projnet)] ≤ projnet[v,:tmax] 
)

projnet[:,:x] = jumpmod[:x_j]

@variable(
    jumpmod, 
    0 ≤ y_j[v ∈ vertices(projnet)]
)

projnet[:,:y] = jumpmod[:y_j]

# final time constraint
nt = only(incident(projnet, :end, :label))
@constraint(
    jumpmod,
    final_time,
    projnet[nt,:y] ≤ T
)

# precedence constraints
for j in vertices(projnet)
    pred = collect(inneighbors(projnet,j))
    if length(j) > 0
        @constraint(
            jumpmod,
            [i ∈ pred],
            projnet[j,:y] ≥ projnet[i,:y] + projnet[i,:x]
        )
    end
end

# minimize costs
@objective(
    jumpmod,
    Max,
    sum(projnet[v,:Δc] * projnet[v,:x] for v in vertices(projnet))
)

optimize!(jumpmod)

projnet[:,:x] = value.(projnet[:,:x])
projnet[:,:y] = value.(projnet[:,:y])