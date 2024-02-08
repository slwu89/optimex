using Catlab
using Distributions
using Plots

@present ProductionSch(FreeSchema) begin
    (M,F,A,Mconv,Fin,Fout,Fact,Ain,Aout)::Ob
    # projections from M_conv
    m_from_conv::Hom(Mconv,M)
    m_to_conv::Hom(Mconv,M)
    # projections from F_in
    f_Fin_::Hom(Fin,F)
    m_Fin::Hom(Fin,M)
    # projections from F_out
    f_Fout::Hom(Fout,F)
    m_Fout::Hom(Fout,M)
    # projections from F_act
    f_Fact::Hom(Fact,F)
    a_Fact::Hom(Fact,A)
    # projections from A_in
    f_Ain::Hom(Ain,F)
    m_Ain::Hom(Ain,M)
    a_Ain::Hom(Ain,A)
    # projections from A_out
    f_Aout::Hom(Aout,F)
    m_Aout::Hom(Aout,M)
    a_Aout::Hom(Aout,A)

    # data attributes
    RealType::AttrType
    NameType::AttrType
    mat_name::Attr(M,NameType)
    buy_min::Attr(M,RealType)
    buy_opt::Attr(M,RealType)
    buy_max::Attr(M,RealType)
    buy_cost::Attr(M,RealType)
    sell_min::Attr(M,RealType)
    sell_opt::Attr(M,RealType)
    sell_max::Attr(M,RealType)
    sell_cost::Attr(M,RealType)

    fac_name::Attr(F,NameType)
    cap_min::Attr(F,RealType)
    cap_max::Attr(F,RealType)

    in_min::Attr(Fin,RealType)
    in_opt::Attr(Fin,RealType)
    in_max::Attr(Fin,RealType)

    out_min::Attr(Fout,RealType)
    out_opt::Attr(Fout,RealType)
    out_max::Attr(Fout,RealType)

    conv_yield::Attr(Mconv,RealType)
    conv_cost::Attr(Mconv,RealType)
    conv_opt::Attr(Mconv,RealType)

    act_min::Attr(Fact,RealType)
    act_opt::Attr(Fact,RealType)
    act_max::Attr(Fact,RealType)
    act_cost::Attr(Fact,RealType)
    act_cap_rate::Attr(Fact,RealType)

    act_in_rate::Attr(Ain,RealType)

    act_out_rate::Attr(Aout,RealType)

    act_name::Attr(A,NameType)
end

Catlab.to_graphviz(ProductionSch,graph_attrs=Dict(:dpi=>"60",:size=>"8",:ratio=>"expand"))

# we need a synthetic data generator
@acset_type ProductionData(ProductionSch)



get_gamma(μ,σ) = begin
    θ = (σ^2)/μ
    k = μ/θ
    return k, θ
end

k,θ = get_gamma(100,10)

histogram(rand(Gamma(k,θ),1000))

production_data = ProductionData{Float64,String}()

# set M
add_parts!(
    production_data, :M, 9,
    mat_name=["raw" .* string.(1:3); "intermediate" .* string.(1:3); "final" .* string.(1:3)],
    buy_min=[rand(Poisson(10),3); zeros(6)],
    buy_max=[rand(Poisson(1000),3); zeros(6)],
    buy_cost=[rand(Gamma(get_gamma(100,10)...),3); zeros(6)],
    sell_min=[zeros(6); rand(Poisson(20),3)],
    sell_max=[zeros(6); rand(Poisson(1100),3)],
    sell_cost=[zeros(6); rand(Gamma(get_gamma(100,20)...),3)]
)

# set Mconv ⊆ M × M
add_parts!(
    production_data, :Mconv, 2,
    m_from_conv=vcat(incident(production_data, ["final2","final3"], :mat_name)...),
    m_to_conv=vcat(incident(production_data, ["final3","raw3"], :mat_name)...),
    conv_yield=rand(Uniform(0.5,0.9),3),
    conv_cost=rand(Gamma(get_gamma(20,5)...),3)
)

# set F
add_parts!(
    production_data, :F, 3,
    fac_name=["A","B","C"],
    cap_min=[10,5,5],
    cap_max=[1000,500,800]
)

# set A
add_parts!(
    production_data, :A
)




# schema for general LP model
@present LPSch(FreeSchema) begin
    # I: constraints, J: variables, C: subset of IxJ, coefficient nonzeros
    (I,J,C)::Ob
    i::Hom(C,I)
    j::Hom(C,J)

    RealType::AttrType

    # min/max values for each constraint
    row_min::Attr(I,RealType)
    row_max::Attr(I,RealType)

    # values for decision variables
    col_profit::Attr(J,RealType)
    col_min::Attr(J,RealType)
    col_optimal::Attr(J,RealType)
    col_max::Attr(J,RealType)

    # (nonzero) coefficient for variable j in constraint i
    coeff_value::Attr(C,RealType)
end

Catlab.to_graphviz(LPSch,graph_attrs=Dict(:dpi=>"72",:size=>"4",:ratio=>"expand"))
