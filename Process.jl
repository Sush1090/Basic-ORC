

function IsentropicCompression(model::EoSModel,T_in,p_in,z,πc,η)
    @assert πc >=1
    @assert η <= 1
    @assert η >0 
    s_in = Entropy(model,p_in,T_in,z)
    h_in = Enthalpy(model,p_in,T_in,z)

    T_out_isen = Tproperty_S(model,p_in*πc,s_in,z)
    h_out_isen = Enthalpy(model,p_in*πc,T_out_isen,z)
    h_out_actual = h_in + (h_out_isen - h_in)/η
    T_out_actual = Tproperty_H(model,p_in*πc,h_out_actual,z)
    return T_out_actual
end

@register_symbolic IsentropicCompression(model::EoSModel,T_in,p_in,z,πc,η)


function IsentropicExpansion(model::EoSModel,T_in,p_in,z,πc,η)
    @assert πc >=1
    @assert η <= 1
    @assert η >0 
    s_in = Entropy(model,p_in,T_in,z)
    h_in = Enthalpy(model,p_in,T_in,z)

    T_out_isen = Tproperty_S(model,p_in/πc,s_in,z) 
    h_out_isen = Enthalpy(model,p_in/πc,T_out_isen,z)
    h_out_actual = h_in - (-h_out_isen + h_in)*η
    T_out_actual = Tproperty_H(model,p_in/πc,h_out_actual,z)
    return T_out_actual
end
@register_symbolic IsentropicExpansion(model::EoSModel,T_in,p_in,z,πc,η)



function ScaleVar(var::Vector,a,b)
    min_ = minimum(var);max_ = maximum(var)
    var_scaled = (var .- min_)./(max_ - min_)
    return var_scaled
end