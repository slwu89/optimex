# --------------------------------------------------------------------------------
# load libraries
using ACSets
using JuMP, HiGHS


# --------------------------------------------------------------------------------
# https://jump.dev/JuMP.jl/stable/tutorials/linear/knapsack/

KnapsackSch = BasicSchema(
    [:Object], 
    [], 
    [:NumAttr], 
    [
        (:profit,:Object,:NumAttr),
        (:weight,:Object,:NumAttr)
    ]
)

@acset_type KnapsackData(KnapsackSch)

knapsackdat = @acset KnapsackData{Float64} begin
    Object = 5
    profit = [5.0, 3.0, 2.0, 7.0, 4.0]
    weight = [2.0, 8.0, 4.0, 2.0, 5.0]
end

capacity = 10.0

# jump model
jumpmod = JuMP.Model(HiGHS.Optimizer)

# decision variable is which things chosen
@variable(
    jumpmod, 
    x[parts(knapsackdat, :Object)],
    Bin
)

@constraint(
    jumpmod,
    sum(knapsackdat[i,:weight] * x[i] for i in parts(knapsackdat, :Object)) <= capacity
)

@objective(
    jumpmod, 
    Max, 
    sum(knapsackdat[i,:profit] * x[i] for i in parts(knapsackdat, :Object))
)

optimize!(jumpmod)

value.(x)
objective_value(jumpmod)


# --------------------------------------------------------------------------------
# https://jump.dev/JuMP.jl/stable/tutorials/linear/diet/

DietSch = BasicSchema(
    [:Object], 
    [], 
    [
        :NameAttr,
        :NumAttr
    ], 
    [
        (:name,:Object,:NameAttr),
        (:cost,:Object,:NumAttr),
        (:calories,:Object,:NumAttr),
        (:protein,:Object,:NumAttr),
        (:fat,:Object,:NumAttr),
        (:sodium,:Object,:NumAttr)
    ]
)

@acset_type DietData(DietSch)

knapsackdat = @acset DietData{String,Float64} begin
    Object = 9
    name = ["hamburger","chicken","hot dog","fries","macaroni","pizza","salad","milk","ice cream"]
    cost = [2.49,2.89,1.50,1.89,2.09,1.99,2.49,0.89,1.59]
    calories = [410,420,560,380,320,320,320,100,330]
    protein = [24,32,20,4,12,15,31,8,8]
    fat = [26,10,32,19,10,12,12,2.5,10]
    sodium = [730,1190,1800,270,930,820,1230,125,180]
end













