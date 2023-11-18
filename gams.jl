# the infamous GAMS blog post: https://www.gams.com/blog/2023/07/performance-in-optimization-models-a-comparative-analysis-of-gams-pyomo-gurobipy-and-jump/
# the JuMP dev respose blog post: https://jump.dev/2023/07/20/gams-blog/
# the original thread in discourse: https://discourse.julialang.org/t/performance-julia-jump-vs-python-pyomo/92044/9
# the JuMP dev announcement in discourse: https://discourse.julialang.org/t/jump-developers-response-to-benchmarks-by-gams/101920
# the GH repo where it can be reproduced: https://github.com/justine18/performance_experiment

using DataFrames
using Distributions
using BenchmarkTools

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
@benchmark let 
    x_list = [
        (i, j, k, l, m)
        for (i, j, k) in IJK_sparse
        for (jj, kk, l) in JKL_sparse if jj == j && kk == k
        for (kkk, ll, m) in KLM_sparse if kkk == k && ll == l
    ]
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[x_list] >= 0)
    @constraint(
        model,
        [i in I], 
        sum(x[k] for k in x_list if k[1] == i) >= 0
    )
    optimize!(model)
end

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
model = JuMP.Model()
set_silent(model)
ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
for df in DataFrames.groupby(ijklm, :i)
    @constraint(model, sum(df.x) >= 0)
end
# optimize!(model)

# ------------------
# the acsets version

using ACSets

IJKLMSch = BasicSchema(
    [:I,:J,:K,:L,:M,:IJK,:JKL,:KLM], 
    [
        (:IJK_I,:IJK,:I),
        (:IJK_J,:IJK,:J),
        (:IJK_K,:IJK,:K),
        (:JKL_J,:JKL,:J),
        (:JKL_K,:JKL,:K),
        (:JKL_L,:JKL,:L),
        (:KLM_K,:KLM,:K),
        (:KLM_L,:KLM,:L),
        (:KLM_M,:KLM,:M)
    ], 
    [:IntAttr], 
    [
        (:value_ijk,:IJK,:IntAttr),
        (:value_jkl,:JKL,:IntAttr),
        (:value_klm,:KLM,:IntAttr)
    ]
)

@acset_type IJKLMData(IJKLMSch, index=[:IJK_I,:IJK_J,:IJK_K,:JKL_J,:JKL_K,:JKL_L,:KLM_K,:KLM_L,:KLM_M])

n=8 # something large
m=10 # 20

# the basic sets
I = collect(1:n)
J = collect(1:m)
K = collect(1:m)
L = collect(1:m)
M = collect(1:m)

# make the data
ijklm_dat = IJKLMData{Int}()

add_parts!(ijklm_dat, :I, length(I))
add_parts!(ijklm_dat, :J, length(J))
add_parts!(ijklm_dat, :K, length(K))
add_parts!(ijklm_dat, :L, length(L))
add_parts!(ijklm_dat, :M, length(M))

# add the relations...first as product...then sparsify

# IJK
ijk_prod = Iterators.product(I,J,K)

add_parts!(
    ijklm_dat, 
    :IJK,
    length(ijk_prod),
    IJK_I=vec([e[1] for e in ijk_prod]),
    IJK_J=vec([e[2] for e in ijk_prod]),
    IJK_K=vec([e[3] for e in ijk_prod]),
    value_ijk=SampleBinomialVec(I,J,K)
)

rem_parts!(
    ijklm_dat, 
    :IJK, 
    findall(ijklm_dat[:,:value_ijk] .== 0)
)

# JKL
jkl_prod = Iterators.product(J,K,L)

add_parts!(
    ijklm_dat, 
    :JKL,
    length(jkl_prod),
    JKL_J=vec([e[1] for e in jkl_prod]),
    JKL_K=vec([e[2] for e in jkl_prod]),
    JKL_L=vec([e[3] for e in jkl_prod]),
    value_jkl=SampleBinomialVec(J,K,L)
)

rem_parts!(
    ijklm_dat, 
    :JKL, 
    findall(ijklm_dat[:,:value_jkl] .== 0)
)

# KLM
klm_prod = Iterators.product(K,L,M)

add_parts!(
    ijklm_dat, 
    :KLM,
    length(klm_prod),
    KLM_K=vec([e[1] for e in klm_prod]),
    KLM_L=vec([e[2] for e in klm_prod]),
    KLM_M=vec([e[3] for e in klm_prod]),
    value_klm=SampleBinomialVec(K,L,M)
)

rem_parts!(
    ijklm_dat, 
    :KLM, 
    findall(ijklm_dat[:,:value_klm] .== 0)
)

connected_paths_query = @relation (i=i,j=j,k=k,l=l,m=m) begin
    IJK(IJK_I=i, IJK_J=j, IJK_K=k)
    JKL(JKL_J=j, JKL_K=k, JKL_L=l)
    KLM(KLM_K=k, KLM_L=l, KLM_M=m)
end

ijklm = query(ijklm_dat, connected_paths_query)



# we want to check this against the DataFrames version

IJK_sparse_df = DataFrame(tables(ijklm_dat).IJK)
select!(IJK_sparse_df, Not(:value_ijk))
rename!(IJK_sparse_df, [:i,:j,:k])

JKL_sparse_df = DataFrame(tables(ijklm_dat).JKL)
select!(JKL_sparse_df, Not(:value_jkl))
rename!(JKL_sparse_df, [:j,:k,:l])

KLM_sparse_df = DataFrame(tables(ijklm_dat).KLM)
select!(KLM_sparse_df, Not(:value_klm))
rename!(KLM_sparse_df, [:k,:l,:m])

ijklm_df = DataFrames.innerjoin(
    DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
    KLM_sparse_df;
    on = [:k, :l],
)

sort(ijklm) == sort(ijklm_df)




# # remove redundant things in I,J,K,L,M

# # find parts of `ob` that are not hit by any incoming homs and remove them
# function not_in_image(acs, ob)
#     in_homs = homs(acset_schema(acs), to=ob, just_names=true)    
#     hit_by = [length.(incident(acs, parts(acs, ob), f)) for f in in_homs]
#     to_rem = Int[]
#     for i in parts(acs, ob)
#         if all([h[i] for h in hit_by] .== 0)
#             push!(to_rem, i)
#         end
#     end
#     return to_rem
# end

# # for I, find things with no in-homs
# not_in_image(ijklm_dat, :I)

# # for J, find things with no in-homs
# not_in_image(ijklm_dat, :J)

# # for K, find things with no in-homs
# not_in_image(ijklm_dat, :K)

# # L
# not_in_image(ijklm_dat, :L)

# # M
# not_in_image(ijklm_dat, :M)