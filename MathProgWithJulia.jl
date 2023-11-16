using JuMP, HiGHS
using LinearAlgebra

# LP programming

# ------------------------------------------------------------
# incredible chairs 1
IC = Model(HiGHS.Optimizer)

@variable(IC, xA≥0)
@variable(IC, xB≥0)
@objective(IC, Max, 4*xA + 6*xB)
IC[:line1] = @constraint(IC, 2*xA≤14)
IC[:line2] = @constraint(IC, 3*xB≤15)
IC[:line3] = @constraint(IC, 4*xA+3*xB≤36)

print(IC)
optimize!(IC)
value(xA)
value(xB)

# ------------------------------------------------------------
# incredible chairs 2
chairs = string.(collect(Iterators.take('a':'z', 10)))
C = length(chairs)
prodlines = collect(1:5)
P = length(prodlines)

Profit = [6,5,9,5,6,3,4,7,4,3]
Capacity = [47,19,36,13,46]
ResourceUsage = [
    6  4 2  3 1 10 2 9  3 5;
    5  6 1  1 7  2 9 1  8 6;
    8 10 7  2 9  6 9 6  5 6;
    8  4 8 10 5  4 1 5  3 5;
    1  4 7  2 4  1 2 3 10 1]
@assert size(ResourceUsage) == (P,C)

IC = Model(HiGHS.Optimizer)

@variable(IC, x[c=1:C]≥0)

# prettier with linear algebra
IC[:prod_capacity] = @constraint(IC, ResourceUsage * x ≤ Capacity)
@objective(IC, Max, Profit ⋅ x)

optimize!(IC)
value.(x)
objective_value(IC)

# ------------------------------------------------------------
# class jobs
Wish = [
    1 3 2 5 5;
    5 2 1 1 2;
    1 5 1 1 1;
    4 5 4 4 4;
    3 5 3 5 3
]
J = 5
C = 5

@assert size(Wish) == (J,C)

model = Model(HiGHS.Optimizer)

@variable(model, Assignment[1:J,1:C], Bin)

model[:each_child_has_one_job] = @constraint(model, sum(Assignment,dims=1) .== 1)
model[:each_job_has_one_child] = @constraint(model, sum(Assignment,dims=2) .== 1)

@objective(model, Max, sum(Wish .* Assignment))

optimize!(model)
objective_value(model)
value.(Assignment)

# ------------------------------------------------------------
# chair transport
P = 2
D = 4
Distances = [
    137 92  48  173;
    54  109 111 85
]
Costs = Distances .* 0.0375

Capacity = [7500, 8500]
Demand = [3250, 3500, 3500, 3000]

model = Model(HiGHS.Optimizer)

# how each plant makes for each depot
@variable(model, x[1:P,1:D] ≥ 0)

# cannot produce more than capacity for each plant
model[:capacity] = @constraint(model, sum(x,dims=2) .≤ Capacity)
# must fulfill demand
model[:demand] = @constraint(model, sum(x,dims=1) .== transpose(Demand))

@objective(model, Min, sum(Costs .* x))

optimize!(model)
objective_value(model)
value.(x)

# ------------------------------------------------------------
# jewellery production
Profit = [50, 35, 85, 60, 55]
Machine = [
    7 0 0 9 0;
    5 7 11 0 5;
    0 3 8 15 3
]
Assembly = [12, 3, 11, 9, 6]

M = 3
N = 5
E = 2

model = Model(HiGHS.Optimizer)

@variable(model, x[1:N] ≥ 0 )
model[:machine_work] = @constraint(model, Machine * x .≤ 7.5*60)
model[:human_work] = @constraint(model, Assembly ⋅ x .≤ 7.5*60*2)
@objective(model, Max, x ⋅ Profit)

optimize!(model)
objective_value(model)

# ------------------------------------------------------------
# micro brewery

Demand = [15,30,25,55,75,115,190,210,105,65,20,20]

model = Model(HiGHS.Optimizer)

@variable(model, 0 ≤ storage[1:12] ≤ 200)
@variable(model, 0 ≤ production[1:12] ≤ 120)

model[:balance_constraint] = @constraint(model, [m=1:12], storage[m] == (m>1 ? storage[m-1] : 0) + production[m] - Demand[m])

# storage cost is $1/liter, so overall cost is just the sum
@objective(model, Min, sum(storage))

optimize!(model)
objective_value(model)
value.(storage)
value.(production)

# ------------------------------------------------------------
# chair distribution 2

R = 6
P = 2
D = 4
UnitTransportCost = 0.0375

DistancesProdDepot = [
    137 92  48  173;
    54  109 111 85
]
CostsProdDepot = DistancesProdDepot .* UnitTransportCost

DistancesProdRetail = [
    307 260 215 196 148 268;
    234 173 194 264 204 218
]
CostsProdRetail = DistancesProdRetail .* UnitTransportCost

DistancesDepotRetail = [
    109 58 65 187 128 88;
    214 163 54 89 26 114;
    223 173 97 71 29 162;
    81 51 133 239 170 155;    
]
CostsDepotRetail = DistancesDepotRetail .* UnitTransportCost

CapacityProd = [7500, 8500]
CapacityDepot = [3250, 3500, 3500, 3000]

DemandRetail = [1500, 2500, 2000, 3000, 2000, 3000]

@assert size(DistancesProdDepot) == (P,D)
@assert size(DistancesProdRetail) == (P,R)
@assert size(DistancesDepotRetail) == (D,R)

model = Model(HiGHS.Optimizer)

# how much each plant makes
@variable(model, 0 ≤ chairs[i=1:P] ≤ CapacityProd[i])
@variable(model, Prod2Depot[1:P,1:D] ≥ 0)
@variable(model, Prod2Retail[1:P,1:R] ≥ 0)
@variable(model, Depot2Retail[1:D,1:R] ≥ 0)

# amount shipped out should equal amount produced
model[:prod_consistency] = @constraint(model, sum(Prod2Depot,dims=2) + sum(Prod2Retail,dims=2) .== chairs)

# amount shipped to depots should not exceed capacity of the depot
model[:depot_capacity] = @constraint(model, sum(Prod2Depot,dims=1) .≤ transpose(CapacityDepot))

# depots can only ship out what was shipped to them
model[:depot_consistency] = @constraint(model, sum(Prod2Depot,dims=1) .== transpose(sum(Depot2Retail,dims=2)))

# products at depots must come from production sites
model[:prod_to_depot] = @constraint(model, sum(Depot2Retail) == sum(Prod2Depot))

# amount at retail needs to equal or exceed their demand
model[:retail_demand] = @constraint(model, sum(Prod2Retail,dims=1) + sum(Depot2Retail,dims=1) .≥ transpose(DemandRetail))

# objective: minimize costs
@objective(model, Min, sum(CostsProdDepot .* Prod2Depot) + sum(CostsProdRetail .* Prod2Retail) + sum(CostsDepotRetail .* Depot2Retail))

# optimize
optimize!(model)
objective_value(model)
value.(chairs)
value.(Prod2Depot)
value.(Prod2Retail)
value.(Depot2Retail)

# MIP programming

# ------------------------------------------------------------
# IC with integer

model = Model(HiGHS.Optimizer)

@variable(model, xA≥0, Int)
@variable(model, xB≥0, Int)

model[:plant1] = @constraint(model, 2*xA ≤ 14)
model[:plant2] = @constraint(model, 3*xB ≤ 15)
model[:plant3] = @constraint(model, 4*xA + 3*xB ≤ 36)

@objective(model, Max, 4*xA + 6*xB)

optimize!(model)
objective_value(model)
value(xA)
value(xB)

# ------------------------------------------------------------
# 6.2 Chair Logistics

R = 6
P = 2
D = 4
UnitTransportCap = 40
UnitTransportCost = 1.5

DistancesProdDepot = [
    137 92  48  173;
    54  109 111 85
]
CostsProdDepot = DistancesProdDepot .* UnitTransportCost

DistancesProdRetail = [
    307 260 215 196 148 268;
    234 173 194 264 204 218
]
CostsProdRetail = DistancesProdRetail .* UnitTransportCost

DistancesDepotRetail = [
    109 58 65 187 128 88;
    214 163 54 89 26 114;
    223 173 97 71 29 162;
    81 51 133 239 170 155;    
]
CostsDepotRetail = DistancesDepotRetail .* UnitTransportCost

CapacityProd = [7500, 8500]
CapacityDepot = [3250, 3500, 3500, 3000]

DemandRetail = [1500, 2500, 2000, 3000, 2000, 3000]

@assert size(DistancesProdDepot) == (P,D)
@assert size(DistancesProdRetail) == (P,R)
@assert size(DistancesDepotRetail) == (D,R)

model = Model(HiGHS.Optimizer)

# how much each plant makes
@variable(model, 0 ≤ chairs[i=1:P] ≤ CapacityProd[i])

# movement of chairs
@variable(model, Prod2Depot[1:P,1:D] ≥ 0)
@variable(model, Prod2Retail[1:P,1:R] ≥ 0)
@variable(model, Depot2Retail[1:D,1:R] ≥ 0)

# movement of trucks
@variable(model, TruckProd2Depot[1:P,1:D] ≥ 0, Int)
@variable(model, TruckProd2Retail[1:P,1:R] ≥ 0, Int)
@variable(model, TruckDepot2Retail[1:D,1:R] ≥ 0, Int)

# amount shipped out should equal amount produced
model[:prod_consistency] = @constraint(model, sum(Prod2Depot,dims=2) + sum(Prod2Retail,dims=2) .== chairs)

# amount shipped to depots should not exceed capacity of the depot
model[:depot_capacity] = @constraint(model, sum(Prod2Depot,dims=1) .≤ transpose(CapacityDepot))

# depots can only ship out what was shipped to them
model[:depot_consistency] = @constraint(model, sum(Prod2Depot,dims=1) .== transpose(sum(Depot2Retail,dims=2)))

# products at depots must come from production sites
model[:prod_to_depot] = @constraint(model, sum(Depot2Retail) == sum(Prod2Depot))

# amount at retail needs to equal or exceed their demand
model[:retail_demand] = @constraint(model, sum(Prod2Retail,dims=1) + sum(Depot2Retail,dims=1) .≥ transpose(DemandRetail))

# shipping of chairs must respect the constraints of having integer number of trucks available
model[:truck_prod_to_depot] = @constraint(model, Prod2Depot .≤ TruckProd2Depot .* UnitTransportCap)
model[:truck_prod_to_depot] = @constraint(model, Prod2Retail .≤ TruckProd2Retail .* UnitTransportCap)
model[:truck_prod_to_depot] = @constraint(model, Depot2Retail .≤ TruckDepot2Retail .* UnitTransportCap)

# objective: minimize transportation costs
@objective(model, Min, sum(CostsProdDepot .* TruckProd2Depot) + sum(CostsProdRetail .* TruckProd2Retail) + sum(CostsDepotRetail .* TruckDepot2Retail))

# optimize
optimize!(model)
objective_value(model)
value.(chairs)
value.(Prod2Depot)
value.(Prod2Retail)
value.(Depot2Retail)

# ------------------------------------------------------------
# 6.3 Class jobs 2

Time = [
    1 2 1 4 4;
    6 2 4 2 2;
    3 3 2 4 4;
    1 1 4 4 2;
    7 2 2 3 1
]

Wish = [
    1 3 2 5 5;
    5 2 1 1 2;
    1 5 1 1 1;
    4 5 4 4 4;
    3 5 3 5 3
]
J = 5
C = 5

@assert size(Wish) == (J,C)
@assert size(Time) == (J,C)

model = Model(HiGHS.Optimizer)

@variable(model, Assignment[1:J,1:C], Bin)

model[:each_child_has_one_job] = @constraint(model, sum(Assignment,dims=1) .== 1)
model[:each_job_has_one_child] = @constraint(model, sum(Assignment,dims=2) .== 1)
model[:job_time] = @constraint(model, Assignment .* Time .≤ 3)

@objective(model, Max, sum(Wish .* Assignment))

optimize!(model)
objective_value(model)
value.(Assignment)

# ------------------------------------------------------------
# 6.4 Startup fund

InvestDetails = [
    1.7 1.4 1.3 2.1 1.9 1.8 1.5 2.2 1.6;
    17  25  19  25  28  23  29  31  18
]

AvailableFunds = 100

model = Model(HiGHS.Optimizer)

@variable(model, Invest[1:9], Bin)

model[:capital_constraints] = @constraint(model, sum(Invest .* InvestDetails[2,:]) ≤ AvailableFunds)

@objective(model, Max, sum(Invest .* InvestDetails[1,:]))

optimize!(model)
objective_value(model)
value.(Invest)

# modify problem
# 1. cannot simultaneously invest in 1 and 4
# 2. cannot simultaneously invest in 6 and 8
# 3. 6 OR 9 can only be invested in if 2 OR 3 is invested in
model = Model(HiGHS.Optimizer)

@variable(model, Invest[1:9], Bin)
model[:capital_constraints] = @constraint(model, sum(Invest .* InvestDetails[2,:]) ≤ AvailableFunds)
@objective(model, Max, sum(Invest .* InvestDetails[1,:]))

model[:OneOrFour] = @constraint(model, Invest[1] + Invest[4] ≤ 1)
model[:SixOrEight] = @constraint(model, Invest[6] + Invest[8] <= 1)
model[:SixOrNine] = @constraint(model, Invest[6] + Invest[9] <= 2*(Invest[2] + Invest[3]))

optimize!(model)
objective_value(model)
value.(Invest)

# --------------------------------------------------------------------------------
# chair distribution with storage


# distribution with storage in a smaller example
Times = 1:12
R = 4
P = 2
D = 3
UnitTransportCost = 0.0375

DistancesProdDepot = [
    137 92  48 ;
    54  109 111
]
CostsProdDepot = DistancesProdDepot .* UnitTransportCost

DistancesProdRetail = [
    307 260 215 196;
    234 173 194 264
]
CostsProdRetail = DistancesProdRetail .* UnitTransportCost

DistancesDepotRetail = [
    109 58 65 187;
    214 163 54 89;
    223 173 97 71
]
CostsDepotRetail = DistancesDepotRetail .* UnitTransportCost

CapacityProd = [3750, 5000]
CapacityDepot = [3250, 3500, 3500]
CostsDepot = [0.75, 0.5, 1.5]

DemandRetail = [1500, 2500, 2000, 3000]

DemandRetailTime = vcat([transpose(DemandRetail - DemandRetail*mult) for mult in abs.(collect(Times) .- 7) ./ 6]...)
DemandRetailTime = round.(DemandRetailTime)

@assert size(DistancesProdDepot) == (P,D)
@assert size(DistancesProdRetail) == (P,R)
@assert size(DistancesDepotRetail) == (D,R)
@assert size(DemandRetailTime) == (length(Times),R)

model = Model(HiGHS.Optimizer)

# how much each plant makes
@variable(model, 0 ≤ chairs[i=1:P,t=Times] ≤ CapacityProd[i])

# shipping
@variable(model, Prod2Depot[1:P,1:D,Times] ≥ 0)
@variable(model, Prod2Retail[1:P,1:R,Times] ≥ 0)
@variable(model, Depot2Retail[1:D,1:R,Times] ≥ 0)

# storage
@variable(model, DepotStorage[1:D,Times] ≥ 0)

# each time period, amount shipped out should equal amount produced
model[:prod_consistency] = @constraint(
    model, 
    [t=Times], 
    sum(Prod2Depot[:,d,t] for d in 1:D) + sum(Prod2Retail[:,r,t] for r in 1:R) .== chairs[:,t]
)

# depot storage has capacity that cannot be exceeded
model[:depot_capacity] = @constraint(model, [t=Times], DepotStorage[:,t] ≤ CapacityDepot)

# balance equation for storage
model[:depot_balance] = @constraint(
    model, 
    [t=Times], 
    DepotStorage[:,t] == (t>1 ? DepotStorage[:,t-1] : zeros(Int,D))
        + sum(Prod2Depot[p,:,t] for p in 1:P)
        - sum(Depot2Retail[:,r,t] for r in 1:R)
)

# amount at retail needs to equal or exceed their demand
model[:retail_demand] = @constraint(
    model, 
    [t=Times], 
    sum(Prod2Retail[p,:,t] for p in 1:P) 
        + sum(Depot2Retail[d,:,t] for d in 1:D) ≥ DemandRetailTime[t,:]
)

# objective: minimize costs
@objective(
    model, 
    Min,
    sum(CostsProdDepot .* sum(Prod2Depot[:,:,t] for t in Times))
        + sum(CostsProdRetail .* sum(Prod2Retail[:,:,t] for t in Times))
        + sum(CostsDepotRetail .* sum(Depot2Retail[:,:,t] for t in Times))
        + sum(CostsDepot .* sum(DepotStorage[:,t] for t in Times))
)

optimize!(model)
objective_value(model)

value.(DepotStorage)
value.(Prod2Retail)
value.(Depot2Retail)
value.(Prod2Depot)

# total amount recieved by customers is equal to the sum of these two
sum(value.(Depot2Retail)[d,:,:] for d in 1:D)
sum(value.(Prod2Retail)[p,:,:] for p in 1:P)
# it needs to equal
transpose(DemandRetailTime)

# total amount shipped out by producers is equal to amount created
sum(value.(Prod2Depot)[:,d,:] for d in 1:D)
sum(value.(Prod2Retail)[:,r,:] for r in 1:R)
# it needs to equal
value.(chairs)


# --------------------------------------------------------------------------------
# 6.6 scrap removal

item_weights = [35,10,45,53,37,22,26,38,63,17,44,54,62,42,39,51,24,52,46,29]
@assert length(item_weights) == 20

Items = 20
Bags = 10

model = Model(HiGHS.Optimizer)

@variable(model, bag_assignment[1:Bags,1:Items], Bin)
@variable(model, bag_used[1:Bags], Bin)

model[:each_item_goes_in_one_bag] = @constraint(model, sum(bag_assignment,dims=1) .== 1)
model[:bag_weight_limit] = @constraint(model, [b=1:Bags], sum(bag_assignment[b,:] .* item_weights) ≤ bag_used[b] * 100)

# objective is to max money saved, it looks weird because for some reason the following doesn't work (counting bags used)
# sum(bag_assignment,dims=2) .> 0
@objective(model, Min, sum(bag_used)*50)

optimize!(model)
objective_value(model)
value.(bag_assignment)


# --------------------------------------------------------------------------------
# 6.8 micro brewery 2

cafe_orders = [
    35 20 15 45 25  65  40  50  35 85 50 55;
     15 10 20 15 15  55  90  80  25 45  5 30;
     5 20 20 35 35  80   60  30  35 20 20 40
]
init_inventory = [25,65,75]
max_prod = 120
max_storage = 300
cost = 0.1

Time = 12
Beers = 3

model = Model(HiGHS.Optimizer)

@variable(model, 0 ≤ storage[1:Beers,1:Time])
@variable(model, 0 ≤ production[1:Beers,1:Time])
@variable(model, production_choice[1:Beers,1:Time], Bin)

# storage balance
model[:storage_balance] = @constraint(
    model, 
    [b=1:Beers,t=1:Time], 
    storage[b,t] == (t>1 ? storage[b,t-1] : init_inventory[b]) +
        production[b,t] - cafe_orders[b,t]
)

# producing the beer
model[:production_limit] = @constraint(model, [b=1:Beers,t=1:Time], production[b,t] ≤ production_choice[b,t]*max_prod)

# storage constraint
model[:storage_constraint] = @constraint(model, [t=1:Time], sum(storage[b,t] for b in 1:Beers) ≤ max_storage)

# production choice
model[:production_choice] = @constraint(model, [t=1:Time], sum(production_choice[b,t] for b in 1:Beers) ≤ 1)

@objective(model, Min, sum(storage)*cost)

optimize!(model)
objective_value(model)

# ------------------------------------------------------------
# 6.9 Chair Logistics 2

R = 6
P = 2
D = 6
UnitTransportCap = 40
UnitTransportCost = 1.5

DistancesProdDepot = [
    137 92  48  173 88  109;
    54  109 111 85  128 105
]
CostsProdDepot = DistancesProdDepot .* UnitTransportCost

DistancesProdRetail = [
    307 260 215 196 148 268;
    234 173 194 264 204 218
]
CostsProdRetail = DistancesProdRetail .* UnitTransportCost

DistancesDepotRetail = [
    109 58 65 187 128 88;
    214 163 54 89 26 114;
    223 173 97 71 29 162;
    81 51 133 239 170 155;
    223 172 62 80 13 124;
    185 135 31 113 54 96
]
CostsDepotRetail = DistancesDepotRetail .* UnitTransportCost

CapacityProd = [7500, 8500]
CapacityDepot = [3250, 3500, 3500, 3000, 3750, 3250]

DemandRetail = [1500, 2500, 2000, 3000, 2000, 3000]

# related to opening/closing of depots
DepotOpenCost = [0, 0, 0, 0, 5000, 5000]
DepotCloseProfit = [5000, 5000, 5000, 5000, 0, 0]

@assert size(DistancesProdDepot) == (P,D)
@assert size(DistancesProdRetail) == (P,R)
@assert size(DistancesDepotRetail) == (D,R)

model = Model(HiGHS.Optimizer)

# how much each plant makes
@variable(model, 0 ≤ chairs[i=1:P] ≤ CapacityProd[i])

# movement of chairs
@variable(model, Prod2Depot[1:P,1:D] ≥ 0)
@variable(model, Prod2Retail[1:P,1:R] ≥ 0)
@variable(model, Depot2Retail[1:D,1:R] ≥ 0)

# movement of trucks
@variable(model, TruckProd2Depot[1:P,1:D] ≥ 0, Int)
@variable(model, TruckProd2Retail[1:P,1:R] ≥ 0, Int)
@variable(model, TruckDepot2Retail[1:D,1:R] ≥ 0, Int)

# do we open a depot or not
@variable(model, OpenDepot[1:D], Bin)

# constraints

# amount shipped out should equal amount produced
model[:prod_consistency] = @constraint(model, sum(Prod2Depot,dims=2) + sum(Prod2Retail,dims=2) .== chairs)

# amount shipped to depots should not exceed capacity of the depot (if opened)
model[:depot_capacity] = @constraint(model, sum(Prod2Depot,dims=1) .≤ transpose(CapacityDepot .* OpenDepot))

# depots can only ship out what was shipped to them
model[:depot_consistency] = @constraint(model, sum(Prod2Depot,dims=1) .== transpose(sum(Depot2Retail,dims=2)))

# products at depots must come from production sites
model[:prod_to_depot] = @constraint(model, sum(Depot2Retail) == sum(Prod2Depot))

# amount at retail needs to equal or exceed their demand
model[:retail_demand] = @constraint(model, sum(Prod2Retail,dims=1) + sum(Depot2Retail,dims=1) .≥ transpose(DemandRetail))

# shipping of chairs must respect the constraints of having integer number of trucks available
model[:truck_prod_to_depot] = @constraint(model, Prod2Depot .≤ TruckProd2Depot .* UnitTransportCap)
model[:truck_prod_to_depot] = @constraint(model, Prod2Retail .≤ TruckProd2Retail .* UnitTransportCap)
model[:truck_prod_to_depot] = @constraint(model, Depot2Retail .≤ TruckDepot2Retail .* UnitTransportCap)

# objective: minimize transportation costs
@objective(
    model, 
    Min, 
    sum(CostsProdDepot .* TruckProd2Depot) + sum(CostsProdRetail .* TruckProd2Retail) + sum(CostsDepotRetail .* TruckDepot2Retail) +
    sum(DepotOpenCost .* OpenDepot) - sum(DepotCloseProfit .* (1 .- OpenDepot))
)

# optimize
optimize!(model)
objective_value(model)
value.(chairs)
value.(Prod2Depot)
value.(Prod2Retail)
value.(Depot2Retail)
value.(OpenDepot)