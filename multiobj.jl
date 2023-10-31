# ------------------------------------------------------------
# multiobj and other topics

using JuMP, HiGHS
using ACSets
using Catlab

# ------------------------------------------------------------
# goal prog: diamond ex
# from https://link.springer.com/book/10.1007/978-3-642-31054-6
# sec 3.3

DiamondSch = BasicSchema(
    [:Store,:Partition], 
    [
        (:partition,:Store,:Partition)
    ], 
    [:NumAttr,:NameAttr], 
    [
        (:partname,:Partition,:NameAttr),
        (:partalloc,:Partition,:NumAttr),
        (:theft,:Store,:NumAttr)
    ]
)

@acset_type DiamondData(DiamondSch, index = [:partition])

diamonddat = @acset DiamondData{Float64,Symbol} begin
    Store=5
    Partition=2
    partition=[1,1,1,2,2]
    partname=[:Mall,:NotMall]
    partalloc=[0.8,0.2]
    theft=[0.001,0.001,0.09,0.02,0.03]
end

# make the jump model
jumpmod = JuMP.Model(HiGHS.Optimizer)

# d.v. is the quantity of diamonds allocated to each store
@variable(
    jumpmod, 
    diamonds[parts(diamonddat, :Store)] >= 0
)

# absolute constraints
@constraint(
    jumpmod,
    sum(diamonds) >= 1000
)

@constraint(
    jumpmod,
    sum(diamonds) <= 1200
)

@constraint(
    jumpmod,
    diamonds[5] >= 300
)

# desired goal c: the stores in the malls should receive at least 80% of the diamonds, if possible
@variables(
    jumpmod, begin
    d1_n >= 0
    d1_p >= 0
end)

jumpmod[:goal_con_c] = @constraint(
    jumpmod,
    sum(diamonddat[2,:partalloc] * diamonds[incident(diamonddat, 1, :partition)]) -
        sum(diamonddat[1,:partalloc] * diamonds[incident(diamonddat, 2, :partition)]) + 
        d1_n - d1_p == 0
)

# desired goal d: the allocations to the stores in the mall should be equal to each other, if possible
@variables(
    jumpmod, begin
    d2_n >= 0
    d2_p >= 0
    d3_n >= 0
    d3_p >= 0
    d4_n >= 0
    d4_p >= 0
end)

jumpmod[:goal_con_d_x1] = @constraint(
    jumpmod,
    2/3*diamonds[1] - 1/3*diamonds[2] - 1/3*diamonds[3] + d2_n - d2_p == 0
)

jumpmod[:goal_con_d_x2] = @constraint(
    jumpmod,
    -1/3*diamonds[1] + 2/3*diamonds[2] - 1/3*diamonds[3] + d3_n - d3_p == 0
)

jumpmod[:goal_con_d_x3] = @constraint(
    jumpmod,
    -1/3*diamonds[1] - 1/3*diamonds[2] + 2/3*diamonds[3] + d4_n - d4_p == 0
)

# desired goal e: we would like to minimize the expected loss based on the given probabilities
@variables(
    jumpmod, begin
    d5_n >= 0
    d5_p >= 0
end)

jumpmod[:goal_con_e] = @constraint(
    jumpmod,
    sum(subpart(diamonddat, :theft) .* diamonds) + d5_n - d5_p == 0
)

# objective
@objective(
    jumpmod,
    Min,
    50*d5_p + 2(d2_n + d2_p) + 2(d3_n + d3_p) + 2(d4_n + d4_p) + d1_n
)

optimize!(jumpmod)

@assert termination_status(jumpmod) == OPTIMAL

value.(diamonds)
sum(value.(diamonds))
value(d5_p)
