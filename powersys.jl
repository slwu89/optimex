# https://jump.dev/JuMP.jl/stable/tutorials/applications/power_systems/

using JuMP
import DataFrames
import HiGHS
import Plots
import StatsPlots

# economic dispatch

# decision variables
# w: power generated from wind
# g: power generated from thermal

# constraints
# g_min < g < g_max: thermal power generation between min and max
# 0 < w < wf: wind power generation between 0 and forecast
# \sum g + w = df: sum of power generation is equal to (forecasted) demand

# obj
# Min \sum cg * g + cw * w: minimize total costs of production


function ThermalGenerator(
    min::Float64,
    max::Float64,
    fixed_cost::Float64,
    variable_cost::Float64,
)
    return (
        min = min,
        max = max,
        fixed_cost = fixed_cost,
        variable_cost = variable_cost,
    )
end

generators = [
    ThermalGenerator(0.0, 1000.0, 1000.0, 50.0),
    ThermalGenerator(300.0, 1000.0, 0.0, 100.0),
]

WindGenerator(variable_cost::Float64) = (variable_cost = variable_cost,)

wind_generator = WindGenerator(50.0)

function Scenario(demand::Float64, wind::Float64)
    return (demand = demand, wind = wind)
end

scenario = Scenario(1500.0, 200.0)

function solve_economic_dispatch(generators::Vector, wind, scenario)
    # Define the economic dispatch (ED) model
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    # Define decision variables
    # power output of generators
    N = length(generators)
    @variable(model, generators[i].min <= g[i = 1:N] <= generators[i].max)
    # wind power injection
    @variable(model, 0 <= w <= scenario.wind)
    # Define the objective function
    @objective(
        model,
        Min,
        sum(generators[i].variable_cost * g[i] for i in 1:N) +
        wind.variable_cost * w,
    )
    # Define the power balance constraint
    @constraint(model, sum(g[i] for i in 1:N) + w == scenario.demand)
    # Solve statement
    optimize!(model)
    # return the optimal value of the objective function and its minimizers
    return (
        g = value.(g),
        w = value(w),
        wind_spill = scenario.wind - value(w),
        total_cost = objective_value(model),
    )
end

solution = solve_economic_dispatch(generators, wind_generator, scenario);

println("Dispatch of Generators: ", solution.g, " MW")
println("Dispatch of Wind: ", solution.w, " MW")
println("Wind spillage: ", solution.wind_spill, " MW")
println("Total cost: \$", solution.total_cost)


# acsets for organization
using ACSets

PowerSch = BasicSchema(
    [:Generator,:Wind,:Thermal], 
    [
        # inclusions of wind and thermal into generator objs
        (:wind_i,:Wind,:Generator),
        (:thermal_i,:Thermal,:Generator)
    ], 
    [:NumAttr], 
    [
        # attrs for generator
        (:variable_cost,:Generator,:NumAttr),
        # attrs for wind
        (:wind,:Wind,:NumAttr),
        # attrs for thermal
        (:min,:Thermal,:NumAttr),
        (:max,:Thermal,:NumAttr),
        (:fixed_cost,:Thermal,:NumAttr)
    ]
)

@acset_type PowerData(PowerSch, index=[:wind_i,:thermal_i])

pwr = @acset PowerData{Float64} begin
    Generator=3    
    Thermal=2
    Wind=1
    thermal_i=[1,2]
    wind_i=[3]    
    variable_cost=[50,100,50]
    wind=[200]
    min=[0,300]
    max=[1000,1000]
    fixed_cost=[1000,0]
end

demand::Float64 = 1500.0


function solve_economic_dispatch(pwr::PowerData, demand)
    # Define the economic dispatch (ED) model
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # decision variable: pwr from thermal
    @variable(
        model,
        pwr[i,:min] <= g[i=parts(pwr,:Thermal)] <= pwr[i,:max]
    )

    # decision variable: pwr from wind
    @variable(
        model,
        0 <= w[i=parts(pwr,:Wind)] <= pwr[i,:wind]
    )

    # power balance constraint
    @constraint(
        model,
        sum(g[i] for i in parts(pwr,:Thermal)) + sum(w[i] for i in parts(pwr,:Wind)) == demand
    )

    # objective
    @objective(
        model,
        Min,
        sum(g[i] * pwr[i,[:thermal_i,:variable_cost]] for i in parts(pwr,:Thermal)) +
            sum(w[i] * pwr[i,[:wind_i,:variable_cost]] for i in parts(pwr,:Wind))
    )

    optimize!(model)
    # return the optimal value of the objective function and its minimizers
    return (
        g = collect(value.(g)),
        w = collect(value.(w)),
        wind_spill = collect(pwr[:,:wind] - value.(w)),
        total_cost = objective_value(model),
    )
end

sol_acset = solve_economic_dispatch(pwr, demand)

# inefficient use of wind generators

demand_scale_df = DataFrames.DataFrame(;
    demand = Float64[],
    dispatch_G1 = Float64[],
    dispatch_G2 = Float64[],
    dispatch_wind = Float64[],
    spillage_wind = Float64[],
    total_cost = Float64[],
)

function scale_demand(scenario, scale)
    return Scenario(scale * scenario.demand, scenario.wind)
end

for demand_scale in 0.2:0.1:1.4
    new_scenario = scale_demand(scenario, demand_scale)
    sol = solve_economic_dispatch(generators, wind_generator, new_scenario)
    push!(
        demand_scale_df,
        (
            new_scenario.demand,
            sol.g[1],
            sol.g[2],
            sol.w,
            sol.wind_spill,
            sol.total_cost,
        ),
    )
end

demand_scale_df

dispatch_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_G1, :dispatch_G2],
        labels = ["G1" "G2"],
        title = "Thermal Dispatch",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand",
        ylabel = "Dispatch [MW]",
    )
)

wind_plot = StatsPlots.@df(
    demand_scale_df,
    Plots.plot(
        :demand,
        [:dispatch_wind, :spillage_wind],
        labels = ["Dispatch" "Spillage"],
        title = "Wind",
        legend = :bottomright,
        linewidth = 3,
        xlabel = "Demand [MW]",
        ylabel = "Energy [MW]",
    )
)

Plots.plot(dispatch_plot, wind_plot)