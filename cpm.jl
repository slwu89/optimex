using Catlab, DataFrames
# import Catlab.Graphs.BasicGraphs

# the data types
@present SchActivityOnNodeGraph <: SchLabeledGraph begin
  Duration::AttrType
  duration::Attr(V,Duration)
  es::Attr(V,Duration) # earliest possible starting times
  ef::Attr(V,Duration) # earliest possible finishing times
  ls::Attr(V,Duration) # latest possible starting times
  lf::Attr(V,Duration) # latest possible finishing times
  float::Attr(V,Duration) # "float" time of task
end

to_graphviz(SchActivityOnNodeGraph)

@abstract_acset_type AbstractActivityOnNodeGraph <: AbstractGraph

@acset_type ActivityOnNodeGraph(SchActivityOnNodeGraph, index=[:src,:tgt]) <: AbstractActivityOnNodeGraph

"""
    Make a `ActivityOnNodeGraph` acset from a `DataFrame` with columns for
`Activity`, `Predecessor`, and `Duration`.
"""
function make_aon_graph(input::DataFrame)
    g = @acset ActivityOnNodeGraph{Symbol,Int} begin
        V = size(proj_df,1)
        label = proj_df.Activity
        duration = proj_df.Duration
    end

    for r in eachrow(proj_df)
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
function forward_pass!(g::T) where {T<:AbstractActivityOnNodeGraph}
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
function backward_pass!(g::T, toposort) where {T<:AbstractActivityOnNodeGraph}
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
function find_critical_path(g::T) where {T<:AbstractActivityOnNodeGraph}
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

projnet = make_aon_graph(proj_df)
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

projnet = make_aon_graph(proj_df)
to_graphviz(projnet, node_labels=:label)

toposort = forward_pass!(projnet)
backward_pass!(projnet, toposort)
cV, cE = find_critical_path(projnet)

cg = Subobject(projnet, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# "reduce" time for D to 7
proj_df[proj_df.Activity .== :D, :Duration] .= 7
projnet = make_aon_graph(proj_df)

toposort = forward_pass!(projnet)
backward_pass!(projnet, toposort)
cV, cE = find_critical_path(projnet)

cg = Subobject(projnet, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)