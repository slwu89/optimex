using Catlab, DataFrames
using JuMP, HiGHS

include("scheduling.jl")

# --------------------------------------------------------------------------------
# basic CPM examples

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

projnet_cpm = ProjGraphCPM{Symbol,Int}()
copy_parts!(projnet_cpm, projnet)

toposort = forward_pass!(projnet_cpm)
backward_pass!(projnet_cpm, toposort)
cV, cE = find_critical_path(projnet_cpm)

cg = Subobject(projnet_cpm, V=cV, E=cE)
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

projnet_cpm = ProjGraphCPM{Symbol,Int}()
copy_parts!(projnet_cpm, projnet)

toposort = forward_pass!(projnet_cpm)
backward_pass!(projnet_cpm, toposort)
cV, cE = find_critical_path(projnet_cpm)

cg = Subobject(projnet_cpm, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# "reduce" time for D to 7
proj_df[proj_df.Activity .== :D, :Duration] .= 7
projnet = make_ProjGraph(proj_df)

projnet_cpm = ProjGraphCPM{Symbol,Int}()
copy_parts!(projnet_cpm, projnet)

toposort = forward_pass!(projnet_cpm)
backward_pass!(projnet_cpm, toposort)
cV, cE = find_critical_path(projnet_cpm)

cg = Subobject(projnet_cpm, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# --------------------------------------------------------------------------------
# CPM with acceleration as a LP problem

# view the schema
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
projnet_cpm = ProjGraphCPM{Symbol,Int}()
copy_parts!(projnet_cpm, projnet)

projnet_cpm[:,:duration] = projnet[:,:x]

toposort = forward_pass!(projnet_cpm)
backward_pass!(projnet_cpm, toposort)
cV, cE = find_critical_path(projnet_cpm)

cg = Subobject(projnet_cpm, V=cV, E=cE)
to_graphviz(cg, node_labels=:label)

# --------------------------------------------------------------------------------
# basic CPM/scheduling as LP

# ex: sec 4.3 from Ulusoy, G., Hazır, Ö., Ulusoy, G., & Hazır, Ö. (2021). Introduction to Project Modeling and Planning 
proj_df = DataFrame(
    Activity = [:start,:A,:B,:C,:D,:E,:F,:G,:end],
    Predecessor = [
        [], [:start], [:start], [:start], [:A], [:C], [:C], [:D,:B,:E], [:F,:G]
    ],
    Duration = [0,5,6,4,5,3,6,7,0]
)

projnet = make_ProjGraph(proj_df)
to_graphviz(projnet, node_labels=:label)

projnet_lp = ProjGraphLP{Symbol,Int,VarType}()
copy_parts!(projnet_lp, projnet)

jumpmod = JuMP.Model(HiGHS.Optimizer)

# dv: when does each task start?
@variable(
    jumpmod,
    0 ≤ t[v in vertices(projnet_lp)]
)

projnet_lp[:,:t] = jumpmod[:t]

# wlog, start project begins at 0
@constraint(
    jumpmod,
    projnet_lp[1,:t] == 0
)

@constraint(
    jumpmod,
    [e ∈ edges(projnet_lp)],
    projnet_lp[e, (:tgt, :t)] - projnet_lp[e, (:src, :t)] ≥ projnet_lp[e, (:src, :duration)]
)

# minimize overall time
@objective(
    jumpmod,
    Min,
    projnet_lp[nv(projnet_lp),:t]
)

optimize!(jumpmod)
    
projnet_lp[:,:t] = value.(projnet_lp[:,:t])


x

# # --------------------------------------------------------------------------------
# # resource constrained scheduling problems

# proj_df = DataFrame(
#     Activity = [:start,:A,:B,:C,:D,:E,:F,:end],
#     Predecessor = [
#         [], [:start], [:start], [:start], [:A,:B], [:C], [:E], [:D,:F]
#     ],
#     Duration = [0,4,1,2,1,2,2,0]
# )

# projnet = make_ProjGraph(proj_df)
# to_graphviz(projnet, node_labels=:label)

# # optimization requires we do the forward and backward passes
# toposort = forward_pass!(projnet)
# backward_pass!(projnet, toposort)

# # figure out the optimization problem
# projnet_optim = ProjGraphRes{Symbol,Int,Int,VectorVarType}()
# copy_parts!(projnet_optim, projnet)
# projnet_optim[:,:resource] = [0,1,2,1,2,1,3,0]

# jumpmod = JuMP.Model(HiGHS.Optimizer)

# # the time it would take if everything was scheduled sequentially
# up_bound = sum(projnet_optim[:,:duration])

# # the decision vars
# @variable(
#     jumpmod, 
#     g[v,:tmin] ≤ x_j[v ∈ vertices(g)] ≤ g[v,:tmax] 
# )

# # this is probably good to read: https://hal.science/hal-00361395/document
# # the local solver ex https://www.hexaly.com/docs/last/exampletour/resource-constrained-project-scheduling-problem-rcpsp.html