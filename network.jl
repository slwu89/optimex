# from http://home.ubalt.edu/ntsbarsh/opre640a/partIII.htm

using JuMP, HiGHS
using ACSets

# the transportation problem
TransSch = BasicSchema(
    [:Warehouse,:Customer,:Edge], 
    [
        (:src,:Edge,:Customer),
        (:tgt,:Edge,:Customer)
    ], 
    [:NumAttr], 
    [
        (:supply,:Warehouse,:NumAttr),
        (:demand,:Customer,:NumAttr),
        (:cost,:Edge,:NumAttr)
    ]
)

@acset_type TransData(TransSch, index = [:src,:tgt])

transdat = @acset TransData{Int} begin
    Warehouse = 3
    Customer = 4
    Edge = 12
    src = vcat([fill(i,4) for i in 1:3]...)
    tgt = repeat(1:4, 3)
    supply = [1200,1000,800]
    demand = [1100,400,750,750]
    cost = [35,30,40,32,37,40,42,25,40,15,20,28]
end

transjump = JuMP.Model(HiGHS.Optimizer)

# decision variable is amount of stuff shipped on each edge
@variable(
    transjump, 
    shipped[parts(transdat, :Edge)] >= 0
)

# demand is fulfilled
@constraint(
    transjump,
    demand[c=parts(transdat, :Customer)],
    sum(shipped[incident(transdat, c, :tgt)]) >= transdat[c,:demand]
)

# supply is respected
@constraint(
    transjump,
    supply[w=parts(transdat, :Warehouse)],
    sum(shipped[incident(transdat, w, :src)]) <= transdat[w,:supply]
)

# minimize the shipping costs
@objective(
    transjump,
    Min,
    sum(shipped[e] * transdat[e,:cost] for e in parts(transdat, :Edge))
)

optimize!(transjump)

# reminder for me to finish work on https://github.com/AlgebraicJulia/ACSets.jl/issues/17 =)
incident_mult(acs, parts, f) = begin
    return intersect([incident(acs, parts[i], f[i]) for i in eachindex(f)]...)
end

# same as online
shipped_vals = [value(shipped[only(incident_mult(transdat, (w,c), (:src,:tgt)))]) for w=parts(transdat, :Warehouse), c=parts(transdat, :Customer)]