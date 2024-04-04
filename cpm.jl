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
end

to_graphviz(SchActivityOnNodeGraph)

@abstract_acset_type AbstractActivityOnNodeGraph <: AbstractGraph

@acset_type ActivityOnNodeGraph(SchActivityOnNodeGraph, index=[:src,:tgt]) <: AbstractActivityOnNodeGraph

# table 7.1
proj_df = DataFrame(
    Activity = [:A,:B,:C,:D,:E,:F,:G,:H,:I,:J],
    Predecessor = [
        [], [:A], [:A,:B], [:B], [:B,:C], [:C,:D,:E],
        [:D], [:F,:G], [:F,:G], [:I]
    ],
    Duration = [5,3,7,4,6,4,2,9,6,2]
)

projnet = @acset ActivityOnNodeGraph{Symbol,Int} begin
    V = size(proj_df,1) + 2
    label = [:start; proj_df.Activity; :end]
    duration = [0; proj_df.Duration; 0]
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

add_edge!(projnet, only(incident(projnet, :start, :label)), only(incident(projnet, first(proj_df).Activity, :label)))
add_edge!(projnet, only(incident(projnet, last(proj_df).Activity, :label)), only(incident(projnet, :end, :label)))

to_graphviz(projnet, node_labels=:label)

# starting node ES = EF = 0
node = only(incident(projnet, :start, :label))
projnet[node, :es] = 0
projnet[node, :ef] = 0

remaining = setdiff(parts(projnet, :V), node)

collect(outneighbors(projnet, node))