# the infamous GAMS blog post: https://www.gams.com/blog/2023/07/performance-in-optimization-models-a-comparative-analysis-of-gams-pyomo-gurobipy-and-jump/
# the JuMP dev respose blog post: https://jump.dev/2023/07/20/gams-blog/
# the original thread in discourse: https://discourse.julialang.org/t/performance-julia-jump-vs-python-pyomo/92044/9
# the JuMP dev announcement in discourse: https://discourse.julialang.org/t/jump-developers-response-to-benchmarks-by-gams/101920
# the GH repo where it can be reproduced: https://github.com/justine18/performance_experiment

using DataFrames
using Distributions
using JuMP, HiGHS
using Catlab
using BenchmarkTools, MarkdownTables

# --------------------------------------------------------------------------------
# the IJKLM model

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

IJK_sparse_df = filter(x -> x.value == 1, IJK)
select!(IJK_sparse_df, Not(:value))

JKL_sparse_df = filter(x -> x.value == 1, JKL)
select!(JKL_sparse_df, Not(:value))

KLM_sparse_df = filter(x -> x.value == 1, KLM)
select!(KLM_sparse_df, Not(:value))

ijklm_df = DataFrames.innerjoin(
    DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
    KLM_sparse_df;
    on = [:k, :l],
)

# acsets
@present IJKLMSch(FreeSchema) begin
    (I,J,K,L,M,IJK,JKL,KLM)::Ob
    IJK_I::Hom(IJK,I)
    IJK_J::Hom(IJK,J)
    IJK_K::Hom(IJK,K)
    JKL_J::Hom(JKL,J)
    JKL_K::Hom(JKL,K)
    JKL_L::Hom(JKL,L)
    KLM_K::Hom(KLM,K)
    KLM_L::Hom(KLM,L)
    KLM_M::Hom(KLM,M)
    IntAttr::AttrType
    value_ijk::Attr(IJK,IntAttr)
    value_jkl::Attr(JKL,IntAttr)
    value_klm::Attr(KLM,IntAttr)
end

Catlab.to_graphviz(IJKLMSch, graph_attrs=Dict(:dpi=>"72",:ratio=>"expand",:size=>"8"))

@acset_type IJKLMData(IJKLMSch, index=[:IJK_I,:IJK_J,:IJK_K,:JKL_J,:JKL_K,:JKL_L,:KLM_K,:KLM_L,:KLM_M])

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
    value_ijk=IJK.value
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
    value_jkl=JKL.value
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
    value_klm=KLM.value
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

Catlab.to_graphviz(connected_paths_query, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:dpi=>"72",:size=>"3.5",:ratio=>"expand"))

query(ijklm_dat, connected_paths_query)


# --------------------------------------------------------------------------------
# the supply chain model
# translated from https://github.com/justine18/performance_experiment/blob/master/supply_chain/data_generation.py

m=3
n=4

# "sets"
J = ["j$x" for x in 1:m+1]
K = ["k$x" for x in 1:m+1]
L = ["l$x" for x in 1:m+1]
M = ["m$x" for x in 1:m+1]

I = ["f$x" for x in 1:n+1]

share = Int(ceil(length(J) * 0.05))

# IJ
# draw a set of units j about to process product i
IJ = Set([(i,j) for i in I for j in sample(J, share)])
# make sure that every unit j is able to process at least one product i
used_j = Set([j for (i,j) in IJ])
for j in J
    if j ∉ used_j
        push!(
            IJ, 
            (only(sample(I,1)), j)
        )
    end
end


# JK
JK = Set([(j,k) for j in J for k in sample(K, share)])
# make sure that every unit j is able to process at least one product i
used_k = Set([k for (j,k) in JK])
for k in K
    if k ∉ used_k
        push!(
            JK,
            (only(sample(J,1)), k)
        )
    end
end

# IK
df_IJ = DataFrame(IJ, [:i, :j])
df_JK = DataFrame(JK, [:j, :k])
df_IJK = innerjoin(df_IJ, df_JK, on=[:j])
IJK = Set([Tuple(r) for r in eachrow(df_IJK)])
# reduce IJK by around 50%
reduced_IJK = sample(collect(IJK), Int(ceil(length(IJK) * 0.5))) |> Set
IK = Set([(i,k) for (i,j,k) in reduced_IJK])

# KL & LM
KL = Set{Tuple{String,String}}()
LM = Set{Tuple{String,String}}()
for k in K
    for m in M
        ll = sample(L, share)
        for l in ll
            push!(KL, (k,l))
            push!(LM, (l,m))
        end
    end
end
# does every l has a k
used_l = Set([l for (k,l) in KL])
for l in L
    if l ∉ used_l
        push!(KL, (only(sample(K,1)),l))
    end
end
# does every l has an m
used_l = Set([l for (l,m) in LM])
for l in L
    if l ∉ used_l
        push!(LM, (l, only(sample(M,1))))
    end
end

# IL, IM
df_KL = DataFrame(KL, [:k,:l])
df_LM = DataFrame(LM, [:l,:m])
df_IJKLM = innerjoin(
    innerjoin(
        df_IJK, df_KL, 
        on=[:k]
    ),
    df_LM, on=[:l]
)
IJKLM = [Tuple(r) for r in eachrow(df_IJKLM)]
IL = Set([(i, l) for (i, j, k, l, m) in IJKLM])
IM = Set([(i, m) for (i, j, k, l, m) in IJKLM])


# IKL, ILM
IKL = Set([(i, k, l) for (i, j, k, l, m) in IJKLM])
ILM = Set([(i, l, m) for (i, j, k, l, m) in IJKLM])

# Demand
D = Dict([(i,m) => rand(0:100) for (i,m) in IM])





# # --------------------------------------------------------------------------------
# # what about with data migration?
# using DataMigrations

# M = @migration IJKLMSch IJKLMSch begin
#     # each Ob gets a diagram
#     # "set" objects
#     I => I
#     J => J
#     K => K
#     L => L
#     M => M
#     # "relation" objects
#     IJK => @join begin
#         i::I, j::J, k::K, ijk::IJK
#         IJK_I(ijk) == i
#         IJK_J(ijk) == j
#         IJK_K(ijk) == k
#     end
#     JKL => @join begin
#         j::J, k::K, l::K, jkl::JKL
#         JKL_J(jkl) == j
#         JKL_K(jkl) == k
#         JKL_L(jkl) == l
#     end
#     KLM => @join begin
#         k::K, l::L, m::M, klm::KLM
#         KLM_K(klm) == k
#         KLM_L(klm) == l
#         KLM_M(klm) == m
#     end
#     # each Hom gets a morphism of diagrams
#     IJK_I => i⋅id
#     IJK_J => j⋅id
#     IJK_K => k⋅id
#     JKL_J => j⋅id
#     JKL_K => k⋅id
#     JKL_L => l⋅id
#     KLM_K => k⋅id
#     KLM_L => l⋅id
#     KLM_M => m⋅id
#   end