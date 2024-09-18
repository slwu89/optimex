using JuMP
import DataFrames
import HiGHS
import SQLite
import SQLite.DBInterface
import Test

filename = joinpath(@__DIR__, "commodity_nz.db");
db = SQLite.DB(filename)
function get_table(db, table)
    query = DBInterface.execute(db, "SELECT * FROM $table")
    return DataFrames.DataFrame(query)
end

df_shipping = get_table(db, "shipping")
df_products = get_table(db, "products")
df_supply = get_table(db, "supply")
df_demand = get_table(db, "demand")

model = Model(HiGHS.Optimizer)
set_silent(model)

df_shipping.x_flow = @variable(model, x[1:size(df_shipping, 1)] >= 0)
df_shipping

df_supply.x_supply = @variable(model, s[1:size(df_supply, 1)] >= 0)
set_upper_bound.(df_supply.x_supply, df_supply.capacity)
df_supply

df_cost = DataFrames.leftjoin(df_shipping, df_products; on = [:product])
df_cost.flow_cost = df_cost.cost_per_km .* df_cost.distance_km
df_cost

@objective(
    model,
    Min,
    df_cost.flow_cost' * df_shipping.x_flow +
    df_supply.cost' * df_supply.x_supply
);

capacity = 30
for df in DataFrames.groupby(df_shipping, [:origin, :destination])
    @constraint(model, sum(df.x_flow) <= capacity)
end

df_flow_out = DataFrames.DataFrame(
    (node = i.origin, product = i.product, x_flow_out = sum(df.x_flow)) for
    (i, df) in pairs(DataFrames.groupby(df_shipping, [:origin, :product]))
)

df_flow_in = DataFrames.DataFrame(
    (node = i.destination, product = i.product, x_flow_in = sum(df.x_flow))
    for (i, df) in
    pairs(DataFrames.groupby(df_shipping, [:destination, :product]))
)

df = DataFrames.outerjoin(df_flow_in, df_flow_out; on = [:node, :product])

df = DataFrames.leftjoin(
    df,
    DataFrames.select(df_supply, [:origin, :product, :x_supply]);
    on = [:node => :origin, :product],
)

df = DataFrames.leftjoin(
    df,
    DataFrames.select(df_demand, [:destination, :product, :demand]);
    on = [:node => :destination, :product],
)

@constraint(
    model,
    [r in eachrow(df)],
    coalesce(r.x_supply, 0.0) + coalesce(r.x_flow_in, 0.0) -
    coalesce(r.x_flow_out, 0.0) == coalesce(r.demand, 0.0),
);

optimize!(model)
Test.@test is_solved_and_feasible(model)
solution_summary(model)

df_shipping.x_flow = value.(df_shipping.x_flow)
df_supply.x_supply = value.(df_supply.x_supply);

DataFrames.select(
    filter!(row -> row.x_flow > 0.0, df_shipping),
    [:origin, :destination, :product, :x_flow],
)