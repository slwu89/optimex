# ------------------------------------------------------------
# from http://home.ubalt.edu/ntsbarsh/opre640a/partIII.htm

using JuMP, HiGHS
using ACSets

# ------------------------------------------------------------
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

objective_value(transjump)


# ------------------------------------------------------------
# the assignment problem

# please note that the structure of this problem is such that we could have reused the schema
# from the transportation problem, but we make another schema so that the names match nicely

AssignSch = BasicSchema(
    [:Applicant,:Job,:Assignment], 
    [
        (:applicant,:Assignment,:Applicant),
        (:job,:Assignment,:Job)
    ], 
    [:NumAttr], 
    [
        (:cost,:Assignment,:NumAttr)
    ]
)

@acset_type AssignData(AssignSch, index = [:applicant,:job])

assigndat = @acset AssignData{Int} begin
    Applicant = 5
    Job = 5
    Assignment = 25
    applicant = vcat([fill(i,5) for i in 1:5]...)
    job = repeat(1:5, 5)
    cost = [
        10,4,6,10,12,
        11,7,7,9,14,
        13,8,12,14,15,
        14,16,13,17,17,
        19,11,17,20,19
    ]
end

# the jump model
assignjump = JuMP.Model(HiGHS.Optimizer)

# decision variable is assignment of applicants to jobs
@variable(
    assignjump, 
    assignment[parts(assigndat, :Assignment)],
    Bin
)

# each applicant can only be assigned one job
@constraint(
    assignjump,
    one_job[a=parts(assigndat, :Applicant)],
    sum(assignment[incident(assigndat, a, :applicant)]) == 1
)

# each job can only be assigned one applicant
@constraint(
    assignjump,
    one_applicant[j=parts(assigndat, :Job)],
    sum(assignment[incident(assigndat, j, :job)]) == 1
)

# minimize cost
@objective(
    assignjump,
    Min,
    sum(assignment[a] * assigndat[a,:cost] for a in parts(assigndat, :Assignment))
)

optimize!(assignjump)

# check it
assignment_vals = [value(assignment[only(incident_mult(assigndat, (a,j), (:applicant,:job)))]) for a=parts(assigndat, :Applicant), j=parts(assigndat, :Job)]
objective_value(assignjump)