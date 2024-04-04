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

# table 7.1
proj_df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:E,:F,:G,:H,:I,:J,:end],
    Predecessor = [
        [], [:start], [:A], [:A,:B], [:B], [:B,:C], [:C,:D,:E],
        [:D], [:F,:G], [:F,:G], [:I], [:H,:J]
    ],
    Duration = [0,5,3,7,4,6,4,2,9,6,2,0]
)

projnet = @acset ActivityOnNodeGraph{Symbol,Int} begin
    V = size(proj_df,1)
    label = proj_df.Activity
    duration = proj_df.Duration
end

for r in eachrow(proj_df)
    if length(r.Predecessor) > 0
        prednodes = vcat(incident(projnet, r.Predecessor, :label)...)
        node = only(incident(projnet, r.Activity, :label))
        add_edges!(projnet, prednodes, fill(node, length(r.Predecessor)))
    else
        continue
    end
end

to_graphviz(projnet, node_labels=:label)

# forward pass

# starting node ES = EF = 0
projnet[1, :es] = 0
projnet[1, :ef] = 0

sorted = topological_sort(projnet)

for v in sorted[2:end]
    pred = collect(inneighbors(projnet, v))
    projnet[v, :es] = maximum(projnet[pred, :ef])
    projnet[v, :ef] = projnet[v, :es] + projnet[v, :duration]
end

# backward pass
reverse!(sorted)

projnet[sorted[1], :lf] = projnet[sorted[1], :es]
projnet[sorted[1], :ls] = projnet[sorted[1], :es]

for v in sorted[2:end]
    succ = collect(outneighbors(projnet, v))
    projnet[v, :lf] = minimum(projnet[succ, :ls])
    projnet[v, :ls] = projnet[v, :lf] - projnet[v, :duration]
end

# float: delay acceptable for each task w/out delaying completion of project
projnet[:,:float] = projnet[:,:lf] - projnet[:,:ef]

# find critical path(s)
s = 1
g = deepcopy(projnet)

# dfs code
parents = zeros(Int, nv(g))
seen = zeros(Bool, nv(g))
S = [s]
seen[s] = true
parents[s] = s
while !isempty(S)
    v = S[end]
    u = 0
    for n in outneighbors(g, v)
        if !seen[n]
            u = n
            break
        end
    end
    # if all outneighbors of v have been seen already, remove it from S
    if u == 0
        pop!(S)
    else
        # otherwise we have now seen outneighbor n, it's parent is v
        seen[u] = true
        push!(S, u)
        parents[u] = v
    end
end
return parents

# our attempt

# critical activities
ca = incident(g, 0, :float)

collect(inneighbors(g, nv(g)))

g[collect(inneighbors(g, nv(g))), :lf]

# outneighbors who are critical activities
intersect(outneighbors(g, 1), ca)

# set of nodes (algo terminates when this is empty)
s = 1
seen = zeros(Bool, nv(g))
S = [s]
# edges on the critical path
critical_edges = Int[]

while !isempty(S)
    v = pop!(S)
    if !seen[v]
        seen[v] = true
        for w in intersect(outneighbors(g, v), ca)
            if g[v, :lf] == g[w, :es]
                push!(S, w)
                push!(critical_edges, only(edges(g, v, w)))
            end
        end
    end
end

cg = Subobject(g, V=findall(seen), E=critical_edges)

to_graphviz(cg,node_labels=:label)
