using Catlab, DataFrames
using JuMP, HiGHS

include("scheduling.jl")

# look at the schema for projects
to_graphviz(SchProjGraph)

# ex: table 7.1 from Eiselt, H. A., & Sandblom, C. L. (2022). Operations research: A model-based approach.
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

# ex: fig 7.4 from Eiselt, H. A., & Sandblom, C. L. (2022). Operations research: A model-based approach.
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
to_graphviz(SchAccelProjGraph)

# ex: Fig III.12 from Eiselt, H. A., & Sandblom, C. L. (2013). Decision analysis, location models, and scheduling problems.
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

proj_jump = optimize_AccelProjGraph!(projnet, 10)

# migrate the AccelProjGraph to ProjGraph to do CPM
projnet_cpm = ProjGraph{Symbol,Int}()
copy_parts!(projnet_cpm, projnet)

projnet_cpm[:,:duration] = projnet[:,:x]

toposort = forward_pass!(projnet_cpm)
backward_pass!(projnet_cpm, toposort)
cV, cE = find_critical_path(projnet_cpm)

cg = Subobject(projnet_cpm, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)
