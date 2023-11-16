# the infamous GAMS blog post: https://www.gams.com/blog/2023/07/performance-in-optimization-models-a-comparative-analysis-of-gams-pyomo-gurobipy-and-jump/
# the JuMP dev respose blog post: https://jump.dev/2023/07/20/gams-blog/
# the original thread in discourse: https://discourse.julialang.org/t/performance-julia-jump-vs-python-pyomo/92044/9
# the JuMP dev announcement in discourse: https://discourse.julialang.org/t/jump-developers-response-to-benchmarks-by-gams/101920
# the GH repo where it can be reproduced: https://github.com/justine18/performance_experiment

using DataFrames
using Distributions

SampleBinomialVec = function(A,B,C,p=0.05)
    vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    while sum(vec) == 0
        vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    end
    return vec
end

n=100 # something large
m=20 # 20

# Sets IJKLM 
I = ["i$x" for x in 1:n]
J = ["j$x" for x in 1:m]
K = ["k$x" for x in 1:m]
L = ["l$x" for x in 1:m]
M = ["m$x" for x in 1:m]

# make IJK
IJK = DataFrame(Iterators.product(I,J,K))
rename!(IJK, [:i,:j,:k])
IJK.value = SampleBinomialVec(I,J,K)

# make JKL
JKL = DataFrame(Iterators.product(J,K,L))
rename!(JKL, [:j,:k,:l])
JKL.value = SampleBinomialVec(J,K,L)

# make KLM
KLM = DataFrame(Iterators.product(K,L,M))
rename!(KLM, [:k,:l,:m])
KLM.value = SampleBinomialVec(K,L,M)

# make the products just general sparse relations
IJK_sparse = [(x.i, x.j, x.k) for x in eachrow(IJK) if x.value == true]
JKL_sparse = [(x.j, x.k, x.l) for x in eachrow(JKL) if x.value == true]
KLM_sparse = [(x.k, x.l, x.m) for x in eachrow(KLM) if x.value == true]

# make the optim model
using JuMP, HiGHS

# -----------------------------------------------
# the "intuitive" formulation --- the one to beat
x_list = [
    (i, j, k, l, m)
    for (i, j, k) in IJK_sparse
    for (jj, kk, l) in JKL_sparse if jj == j && kk == k
    for (kkk, ll, m) in KLM_sparse if kkk == k && ll == l
]
model = Model(HiGHS.Optimizer)
@variable(model, x[x_list] >= 0)
@constraint(
    model,
    [i in I], 
    sum(x[k] for k in x_list if k[1] == i) >= 0
)
optimize!(model)

# ----------------------
# the DataFrames version
IJK_sparse_df = filter(x -> x.value == 1, IJK)
select!(IJK_sparse_df, Not(:value))

JKL_sparse_df = filter(x -> x.value == 1, JKL)
select!(JKL_sparse_df, Not(:value))

KLM_sparse_df = filter(x -> x.value == 1, KLM)
select!(KLM_sparse_df, Not(:value))

ijklm = DataFrames.innerjoin(
    DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
    KLM_sparse_df;
    on = [:k, :l],
)

# build model
# model = Model(Gurobi.Optimizer)
# set_silent(model)
# ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
# for df in DataFrames.groupby(ijklm, :i)
#     @constraint(model, sum(df.x) >= 0)
# end
# optimize!(model)