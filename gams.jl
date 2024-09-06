using DataFrames
using Distributions
using JuMP, HiGHS
using Catlab, DataMigrations
using BenchmarkTools, MarkdownTables

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
filter!(:value => v -> v != 0, IJK)
select!(IJK, Not(:value))

# make JKL
JKL = DataFrame(Iterators.product(J,K,L))
rename!(JKL, [:j,:k,:l])
JKL.value = SampleBinomialVec(J,K,L)
filter!(:value => v -> v != 0, JKL)
select!(JKL, Not(:value))

# make KLM
KLM = DataFrame(Iterators.product(K,L,M))
rename!(KLM, [:k,:l,:m])
KLM.value = SampleBinomialVec(K,L,M)
filter!(:value => v -> v != 0, KLM)
select!(KLM, Not(:value))

# slow one
@benchmark let 
    x_list = [
        (i, j, k, l, m)
        for (i, j, k) in eachrow(IJK)
        for (jj, kk, l) in eachrow(JKL) if jj == j && kk == k
        for (kkk, ll, m) in eachrow(KLM) if kkk == k && ll == l
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

# DataFrames
@benchmark let
    ijklm = DataFrames.innerjoin(
        DataFrames.innerjoin(IJK, JKL; on = [:j, :k]),
        KLM;
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
end

@acset_type IJKLMData(IJKLMSch, index=[:IJK_I,:IJK_J,:IJK_K,:JKL_J,:JKL_K,:JKL_L,:KLM_K,:KLM_L,:KLM_M])


ijklm_dat = @acset IJKLMData begin
    I = n
    J = m
    K = m
    L = m
    M = m

    IJK = nrow(IJK)
    IJK_I = [parse(Int, i[2:end]) for i in IJK.i]
    IJK_J = [parse(Int, j[2:end]) for j in IJK.j]
    IJK_K = [parse(Int, k[2:end]) for k in IJK.k]

    JKL = nrow(JKL)
    JKL_J = [parse(Int, j[2:end]) for j in JKL.j]
    JKL_K = [parse(Int, k[2:end]) for k in JKL.k]
    JKL_L = [parse(Int, l[2:end]) for l in JKL.l]

    KLM = nrow(KLM)
    KLM_K = [parse(Int, k[2:end]) for k in KLM.k]
    KLM_L = [parse(Int, l[2:end]) for l in KLM.l]
    KLM_M = [parse(Int, m[2:end]) for m in KLM.m]
end

connected_paths_query = @relation (i=i,j=j,k=k,l=l,m=m) begin
    IJK(IJK_I=i, IJK_J=j, IJK_K=k)
    JKL(JKL_J=j, JKL_K=k, JKL_L=l)
    KLM(KLM_K=k, KLM_L=l, KLM_M=m)
end

@benchmark let
    ijklm = query(ijklm_dat, connected_paths_query)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    ijklm[!, :x] = @variable(model, x[1:size(ijklm, 1)] >= 0)
    for df in DataFrames.groupby(ijklm, :i)
        @constraint(model, sum(df.x) >= 0)
    end
    optimize!(model)
end



@present IJKLMRelSch(FreeSchema) begin
    (IJKLM,I,J,K,L,M)::Ob
    i::Hom(IJKLM,I)
    j::Hom(IJKLM,J)
    k::Hom(IJKLM,K)
    l::Hom(IJKLM,L)
    m::Hom(IJKLM,M)
end

@acset_type IJKLMRelType(IJKLMRelSch)

M = @migration IJKLMRelSch IJKLMSch begin
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
end;


@benchmark let
    ijklm = migrate(IJKLMRelType, ijklm_dat, M)
    model = JuMP.Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[parts(ijklm,:IJKLM)] >= 0)
    for i in parts(ijklm,:I)
        @constraint(model, sum(x[incident(ijklm,i,:i)]) >= 0)
    end
    optimize!(model)
end