# the infamous GAMS blog post: https://www.gams.com/blog/2023/07/performance-in-optimization-models-a-comparative-analysis-of-gams-pyomo-gurobipy-and-jump/
# the JuMP dev respose blog post: https://jump.dev/2023/07/20/gams-blog/
# the original thread in discourse: https://discourse.julialang.org/t/performance-julia-jump-vs-python-pyomo/92044/9
# the JuMP dev announcement in discourse: https://discourse.julialang.org/t/jump-developers-response-to-benchmarks-by-gams/101920
# the GH repo where it can be reproduced: https://github.com/justine18/performance_experiment

using DataFrames
using Distributions
using JuMP, HiGHS
using Catlab, DataMigrations
using BenchmarkTools, MarkdownTables


# --------------------------------------------------------------------------------
# synethetic data for the IJKLM model
SampleBinomialVec = function(A,B,C,p=0.05)
    vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    while sum(vec) == 0
        vec = rand(Binomial(1, p), length(A) * length(B) * length(C))
    end
    return vec
end

SampleIJKLM = function(n=100,m=20)
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

    return I,J,K,L,M,IJK,JKL,KLM,IJK_sparse,JKL_sparse,KLM_sparse
end

I,J,K,L,M,IJK,JKL,KLM,IJK_sparse,JKL_sparse,KLM_sparse = SampleIJKLM()

# --------------------------------------------------------------------------------
# the "intuitive" formulation
IntuitiveIJKLM = function(IJK_sparse,JKL_sparse,KLM_sparse)
    x_list = [
        (i, j, k, l, m)
        for (i, j, k) in IJK_sparse
        for (jj, kk, l) in JKL_sparse if jj == j && kk == k
        for (kkk, ll, m) in KLM_sparse if kkk == k && ll == l
    ]
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[x_list] >= 0)
    @constraint(
        model,
        [i in I], 
        sum(x[k] for k in x_list if k[1] == i) >= 0
    )
    optimize!(model)
end

# @benchmark IntuitiveIJKLM(IJK_sparse,JKL_sparse,KLM_sparse)


# --------------------------------------------------------------------------------
# the "dataframes" formulation
sparsify_df = function(IJK, JKL, KLM)
    IJK_sparse_df = filter(x -> x.value == 1, IJK)
    select!(IJK_sparse_df, Not(:value))

    JKL_sparse_df = filter(x -> x.value == 1, JKL)
    select!(JKL_sparse_df, Not(:value))

    KLM_sparse_df = filter(x -> x.value == 1, KLM)
    select!(KLM_sparse_df, Not(:value))

    return IJK_sparse_df, JKL_sparse_df, KLM_sparse_df
end

DataBaseIJKLM = function(IJK_sparse_df, JKL_sparse_df, KLM_sparse_df)
    ijklm = DataFrames.innerjoin(
        DataFrames.innerjoin(IJK_sparse_df, JKL_sparse_df; on = [:j, :k]),
        KLM_sparse_df;
        on = [:k, :l],
    )
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
    for df in DataFrames.groupby(ijklm, :i)
        @constraint(model, sum(df.x) >= 0)
    end
    optimize!(model)
end

# benchmark it
IJK_sparse_df,JKL_sparse_df,KLM_sparse_df = sparsify_df(IJK,JKL,KLM)
# @benchmark DataBaseIJKLM(IJK_sparse_df,JKL_sparse_df,KLM_sparse_df) setup=(DataBaseIJKLM(IJK_sparse_df,JKL_sparse_df,KLM_sparse_df))

# --------------------------------------------------------------------------------
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

make_acset = function(I,J,K,L,M,IJK,JKL,KLM)
    ijklm_dat = IJKLMData{Int}()

    add_parts!(ijklm_dat, :I, length(I))
    add_parts!(ijklm_dat, :J, length(J))
    add_parts!(ijklm_dat, :K, length(K))
    add_parts!(ijklm_dat, :L, length(L))
    add_parts!(ijklm_dat, :M, length(M))

    # IJK
    ijk_prod = Iterators.product(1:length(I),1:length(J),1:length(K))

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
    jkl_prod = Iterators.product(1:length(J),1:length(K),1:length(L))

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
    klm_prod = Iterators.product(1:length(K),1:length(L),1:length(M))

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

    return ijklm_dat
end

ijklm_acs = make_acset(I,J,K,L,M,IJK,JKL,KLM)

# --------------------------------------------------------------------------------
# conjunctive query on acset method

ijklm_query = @relation (i=i,j=j,k=k,l=l,m=m) begin
    IJK(IJK_I=i, IJK_J=j, IJK_K=k)
    JKL(JKL_J=j, JKL_K=k, JKL_L=l)
    KLM(KLM_K=k, KLM_L=l, KLM_M=m)
end

Catlab.to_graphviz(ijklm_query, box_labels=:name, junction_labels=:variable, graph_attrs=Dict(:dpi=>"72",:size=>"3.5",:ratio=>"expand"))

ijklm_query_df = query(ijklm_acs, ijklm_query)

# benchmark it
QueryIJKLM = function(acs, uwd_query)
    ijklm = query(acs, uwd_query)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
    for df in DataFrames.groupby(ijklm, :i)
        @constraint(model, sum(df.x) >= 0)
    end
    optimize!(model)
end

# @benchmark QueryIJKLM(ijklm_acs,ijklm_query) setup=(QueryIJKLM(ijklm_acs,ijklm_query) )

# --------------------------------------------------------------------------------
# data migration

# schema for D
@present IJKLMRelSch(FreeSchema) begin
    (IJKLM,I,J,K,L,M)::Ob
    i::Hom(IJKLM,I)
    j::Hom(IJKLM,J)
    k::Hom(IJKLM,K)
    l::Hom(IJKLM,L)
    m::Hom(IJKLM,M)
end

Catlab.to_graphviz(IJKLMRelSch, graph_attrs=Dict(:dpi=>"72",:ratio=>"expand",:size=>"3.5"))

@acset_type IJKLMRelType(IJKLMRelSch, index=[:i,:j,:k,:l,:m])

# F: D -> Diag^{op}(C)
ijklm_migration = @migration IJKLMRelSch IJKLMSch begin
    IJKLM => @join begin
        ijk::IJK
        jkl::JKL
        klm::KLM
        i::I
        j::J
        k::K
        l::L
        m::M
        IJK_I(ijk) == i
        IJK_J(ijk) == j
        JKL_J(jkl) == j
        IJK_K(ijk) == k
        JKL_K(jkl) == k
        KLM_K(klm) == k
        JKL_L(jkl) == l
        KLM_L(klm) == l
        KLM_M(klm) == m
    end
    I => I
    J => J
    K => K
    L => L
    M => M
    i => i
    j => j
    k => k
    l => l
    m => m
end

# F = functor(ijklm_migration)

# # look at the most important diagram
# to_graphviz(F.ob_map[:IJKLM],node_labels=true)

# ijklm_migrate_acset = migrate(IJKLMRelType, ijklm_dat, ijklm_migration)
# nparts(ijklm_migrate_acset, :IJKLM) == size(ijklm_query_df,1)

# pretty_tables(ijklm_migrate_acset, tables=[:IJKLM], max_num_of_rows=5)


# lens = [length(incident(ijklm_migrate_acset, i, :i)) for i in parts(ijklm_migrate_acset,:I)]


# # benchmark it
# ijklm_migrate_acset = migrate(IJKLMRelType, ijklm_dat, ijklm_migration)
# model = JuMP.Model(HiGHS.Optimizer)
# set_silent(model)
# @variable(model, x[parts(ijklm_migrate_acset,:IJKLM)] >= 0)
# for i in parts(ijklm_migrate_acset,:I)
#     @constraint(model, sum(x[incident(ijklm_migrate_acset,i,:i)]) >= 0)
# end
# optimize!(model)


MigrateIJKLM = function(acs,migration)
    acs_migrate = migrate(IJKLMRelType, acs, migration)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[parts(acs_migrate,:IJKLM)] >= 0)
    for i in parts(acs_migrate,:I)
        @constraint(model, sum(x[incident(acs_migrate,i,:i)]) >= 0)
    end
    optimize!(model)
end

# @benchmark MigrateIJKLM(ijklm_acs, ijklm_migration) setup=(MigrateIJKLM(ijklm_acs, ijklm_migration))

# --------------------------------------------------------------------------------
# comparison of results across range of model sizes

nrange = 10 .^ (1:5)
mrange = Int.(nrange ./ 5)

# we run everything once before starting benchmarking to avoid including compile times, just care about runtime performance
BenchmarkIJKLM = function(in_query,in_migration,n=100,m=20)
    I,J,K,L,M,IJK,JKL,KLM,IJK_sparse,JKL_sparse,KLM_sparse = SampleIJKLM(n,m)

    b_intuit = @benchmark IntuitiveIJKLM(IJK_sparse,JKL_sparse,KLM_sparse) setup=(IntuitiveIJKLM(IJK_sparse,JKL_sparse,KLM_sparse))

    IJK_sparse_df, JKL_sparse_df, KLM_sparse_df = sparsify_df(IJK,JKL,KLM)
    b_dataframes = @benchmark DataBaseIJKLM(IJK_sparse_df,JKL_sparse_df,KLM_sparse_df) setup=(DataBaseIJKLM(IJK_sparse_df,JKL_sparse_df,KLM_sparse_df))

    ijklm_acs = make_acset(I,J,K,L,M,IJK,JKL,KLM)
    b_query = @benchmark QueryIJKLM(ijklm_acs,$(in_query)) setup=(QueryIJKLM(ijklm_acs,$(in_query)))

    b_migrate = @benchmark MigrateIJKLM(ijklm_acs,$(in_migration)) setup=(MigrateIJKLM(ijklm_acs,$(in_migration)))

    return (i=b_intuit,d=b_dataframes,q=b_query,m=b_migrate)
end

# benchmark it
@time ijklm_benchmark_results = [BenchmarkIJKLM(ijklm_query,ijklm_migration,nrange[i],mrange[i]) for i in 1:length(nrange)]

sdfdsfsd=5

# # --------------------------------------------------------------------------------
# # "rebuild" the subacset from the query result
# df_recoded = deepcopy(ijklm_dat_res)

# I = sort(unique(ijklm_dat_res.i))
# J = sort(unique(ijklm_dat_res.j))
# K = sort(unique(ijklm_dat_res.k))
# L = sort(unique(ijklm_dat_res.l))
# M = sort(unique(ijklm_dat_res.m))

# # rename the integers to be 1:size(Set)
# I_recode = Dict([i => ix for (ix,i) in enumerate(I)])
# J_recode = Dict([i => ix for (ix,i) in enumerate(J)])
# K_recode = Dict([i => ix for (ix,i) in enumerate(K)])
# L_recode = Dict([i => ix for (ix,i) in enumerate(L)])
# M_recode = Dict([i => ix for (ix,i) in enumerate(M)])

# df_recoded.i = map(x->I_recode[x],df_recoded.i)
# df_recoded.j = map(x->J_recode[x],df_recoded.j)
# df_recoded.k = map(x->K_recode[x],df_recoded.k)
# df_recoded.l = map(x->L_recode[x],df_recoded.l)
# df_recoded.m = map(x->M_recode[x],df_recoded.m)

# ijklm_dat_rebuild = IJKLMData{Int}()

# add_parts!(ijklm_dat_rebuild, :I, length(I))
# add_parts!(ijklm_dat_rebuild, :J, length(J))
# add_parts!(ijklm_dat_rebuild, :K, length(K))
# add_parts!(ijklm_dat_rebuild, :L, length(L))
# add_parts!(ijklm_dat_rebuild, :M, length(M))

# for r in eachrow(unique(df_recoded[:,[:i,:j,:k]]))
#     add_part!(ijklm_dat_rebuild, :IJK, IJK_I = r.i, IJK_J = r.j, IJK_K = r.k)
# end

# for r in eachrow(unique(df_recoded[:,[:j,:k,:l]]))
#     add_part!(ijklm_dat_rebuild, :JKL, JKL_J = r.j, JKL_K = r.k, JKL_L = r.l)
# end

# for r in eachrow(unique(df_recoded[:,[:k,:l,:m]]))
#     add_part!(ijklm_dat_rebuild, :KLM, KLM_K = r.k, KLM_L = r.l, KLM_M = r.m)
# end

# ijklm_dat_rebuild_res = query(ijklm_dat_rebuild, connected_paths_query)

# # should have same number of results
# size(ijklm_dat_res,1) == size(ijklm_dat_rebuild_res,1)