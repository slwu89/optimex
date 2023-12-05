using Catlab

@present GeneralSch(FreeSchema) begin
    (I,J,C)::Ob
    i::Hom(C,I)
    j::Hom(C,J)
    RealType::AttrType
    row_min::Attr(I,RealType)
    row_max::Attr(I,RealType)
    col_profit::Attr(J,RealType)
    col_min::Attr(J,RealType)
    col_optimal::Attr(J,RealType)
    col_max::Attr(J,RealType)
end

Catlab.to_graphviz(GeneralSch,graph_attrs=Dict(:dpi=>"72",:size=>"4",:ratio=>"expand"))

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
    buy_min::Attr(M,RealType)
    buy_opt::Attr(M,RealType)
    buy_max::Attr(M,RealType)
    buy_cost::Attr(M,RealType)
    sell_min::Attr(M,RealType)
    sell_opt::Attr(M,RealType)
    sell_max::Attr(M,RealType)
    sell_cost::Attr(M,RealType)

    cap_min::Attr(F,RealType)
    cap_max::Attr(F,RealType)

    conv_yield::Attr(Mconv,RealType)
    conv_cost::Attr(Mconv,RealType)
    conv_opt::Attr(Mconv,RealType)

    in_min::Attr(Fin,RealType)
    in_opt::Attr(Fin,RealType)
    in_max::Attr(Fin,RealType)

    out_min::Attr(Fout,RealType)
    out_opt::Attr(Fout,RealType)
    out_max::Attr(Fout,RealType)

    act_min::Attr(Fact,RealType)
    act_opt::Attr(Fact,RealType)
    act_max::Attr(Fact,RealType)
    act_cost::Attr(Fact,RealType)
    act_cap_rate::Attr(Fact,RealType)

    act_in_rate::Attr(Ain,RealType)

    act_out_rate::Attr(Aout,RealType)
end

Catlab.to_graphviz(ProductionSch,graph_attrs=Dict(:dpi=>"60",:size=>"8",:ratio=>"expand"))

# we need a synthetic data generator