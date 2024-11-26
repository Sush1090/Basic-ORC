Enthalpy(model::EoSModel,p,T,z) = Clapeyron.enthalpy(model::EoSModel,p,T,z,phase = "unkown")
@register_symbolic Enthalpy(model::EoSModel,p,T,z)

Entropy(model::EoSModel,p,T,z) = Clapeyron.entropy(model::EoSModel,p,T,z,phase = "unkown")
@register_symbolic Entropy(model::EoSModel,p,T,z)

Tproperty_S(model::EoSModel,p,s,z) = Clapeyron.Tproperty(model::EoSModel,p,s,z,entropy,phase = "unkown")
@register_symbolic Tproperty_S(model::EoSModel,p,s,z)

Tproperty_H(model::EoSModel,p,s,z) = Clapeyron.Tproperty(model::EoSModel,p,s,z,enthalpy,phase = "unkown",verbose = false)
@register_symbolic Tproperty_S(model::EoSModel,p,s,z)

Bubble_temperature(model::EoSModel,p,z) = Clapeyron.bubble_temperature(model,p,z)[1]
@register_symbolic Bubble_temperature(model::EoSModel,p,z)

Dew_temperature(model::EoSModel,p,z) = Clapeyron.dew_temperature(model,p,z)[1]
@register_symbolic Dew_temperature(model::EoSModel,p,z)



function mass_to_moles(model::EoSModel,compositions::AbstractVector,mass)
    @assert isapprox(sum(compositions),1.0,atol = 1e-8)
    mws = Clapeyron.mw(model);
    M = mass*compositions;
    return M./mws
end
@register_symbolic mass_to_moles(model::EoSModel,compositions::AbstractVector,mass)

"""
`1` -> for Gas
`0` -> for Liquid
`-1` -> for Twophase
"""
function ClapyeronPhase_PT(model::EoSModel,p,T,z)
    dt = Dew_temperature(model,p,z)
    bt =  Bubble_temperature(model,p,z)
    if dt <=  T
        return 1;
    end
    if bt<= T <= dt
        return -1;
    end
    if T < bt
        return 0; 
    end
end

function LiquidPhaseChecker(model::EoSModel,p,T,z)
    phase = ClapyeronPhase_PT(model,p,T,z)
    @assert phase == 0 "Phase not Liquid - should be liquid here"
    return phase
end
@register_symbolic LiquidPhaseChecker(model::EoSModel,p,T,z)


function GasPhaseChecker(model::EoSModel,p,T,z)
    phase = ClapyeronPhase_PT(model,p,T,z)
    @assert phase == 1 "Phase not gas - should be gas here"
    return phase
end
@register_symbolic GasPhaseChecker(model::EoSModel,p,T,z)

@connector function fluidPort(;name)
    vars = @variables begin
        p(t), [description = "Pressure (Pa)",input = true]
        T(t), [input = true, description = "Temperature (T)"]
        mdot(t), [input = true, description = "mass flow rate (kg/s)"]
        z(t), [input = true,description  ="Moles vector"]
    end
    ODESystem(Equation[], t, vars, [];name=name)
end

@component function MassSource(;name,composition)
 @named port = fluidPort()
 vars = @variables begin
    s(t)
    T(t)
    p(t)
    h(t)
    z(t)
 end

 para = @parameters begin
    T_source, [description = "Temperature at source (K)"]
    p_source, [description = "pressure at source (Pa)"]
    mdot_source, [description = "mass flow rate at source (kg/s)"]
 end

 eqs = [
    z  ~ mass_to_moles(model,composition,mdot_source)
    T ~ T_source
    p ~ p_source
    h ~ Enthalpy(model,p,T,z)
    s ~ Entropy(model,p,T,z)

    port.T ~ T
    port.p ~ p
    port.mdot ~ mdot_source
    port.z ~ z
 ]

 compose(ODESystem(eqs, t, vars, para;name=name),port)
end

@component function Pump(;name)
    @named inport = fluidPort()
    @named outport = fluidPort()
    vars = @variables begin
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       z_in(t)

       s_out(t)
       T_out(t)
       p_out(t)
       h_out(t)
       z_out(t)

       phase_check_in(t)
       phase_check_out(t)
    end
    para = @parameters begin
        πc = 5,[description = "Pressure ratio (-)"]
        η=1.0, [description = "Isentropic Efficiency (-)"]
     end
     eqs = [

        phase_check_in ~ LiquidPhaseChecker(model,p_in,T_in,z_in)
        z_in  ~ inport.z
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Enthalpy(model,p_in,T_in,z_in)
        s_in ~ Entropy(model,p_in,T_in,z_in)

        z_out ~ z_in
        p_out ~ πc*p_in
        T_out ~ IsentropicCompression(model,T_in,p_in,z_in,πc,η)
        s_out ~ Entropy(model,p_out,T_out,z_out)
        h_out ~ Enthalpy(model,p_out,T_out,z_out)

        outport.T ~ T_out
        outport.p ~ p_out
        outport.mdot ~ inport.mdot
        phase_check_out ~ LiquidPhaseChecker(model,p_out,T_out,z_out)
        outport.z ~ z_out
     ]
    
     compose(ODESystem(eqs, t, vars, para;name=name),inport,outport)
end


@component function Evaporator(;name)
    @named inport = fluidPort()
    @named outport = fluidPort()
    vars = @variables begin
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       z_in(t)

       s_out(t)
       T_out(t)
       p_out(t)
       h_out(t)
       z_out(t)

       phase_check_in(t)
       phase_check_out(t)

    end
    para = @parameters begin
        ΔTsh = 1,[description = "Super heat temperatyre (K)"]
        Δp=0, [description = "Pressure drop at the outlet of evaporator (Pa)"]
     end
     eqs = [
        z_in  ~ inport.z
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Enthalpy(model,p_in,T_in,z_in)
        s_in ~ Entropy(model,p_in,T_in,z_in)
        phase_check_in ~ LiquidPhaseChecker(model,p_in,T_in,z_in)

        z_out ~ z_in
        p_out ~ p_in - Δp
        T_out ~ Dew_temperature(model,p_out,z_out) + ΔTsh
        s_out ~ Entropy(model,p_out,T_out,z_out)
        h_out ~ Enthalpy(model,p_out,T_out,z_out)

        outport.T ~ T_out
        outport.p ~ p_out
        outport.mdot ~ inport.mdot
        phase_check_out ~ GasPhaseChecker(model,p_out,T_out,z_out)
        outport.z ~ z_out
     ]
    
     compose(ODESystem(eqs, t, vars, para;name=name),inport,outport)
end


@component function Expander(;name)
    @named inport = fluidPort()
    @named outport = fluidPort()
    vars = @variables begin
       s_in(t)
       T_in(t)
       p_in(t)
       h_in(t)
       z_in(t)

       s_out(t)
       T_out(t)
       p_out(t)
       h_out(t)
       z_out(t)

       phase_check_in(t)
       phase_check_out(t)
    end
    para = @parameters begin
        πc = 4,[description = "Pressure ratio (-)"]
        η=1.0, [description = "Isentropic Efficiency (-)"]
        pressure_correction = 1, [description = "Pressure addition supplied to expansion to avoid phase transition of values"]
     end
     eqs = [
        z_in  ~ inport.z
        T_in ~ inport.T
        p_in ~ inport.p
        h_in ~ Enthalpy(model,p_in,T_in,z_in)
        s_in ~ Entropy(model,p_in,T_in,z_in)
        phase_check_in ~ GasPhaseChecker(model,p_in,T_in,z_in)

        z_out ~ z_in
        p_out ~ p_in/πc 
        T_out ~ IsentropicExpansion(model,T_in,p_in + pressure_correction,z_in,πc,η)
        s_out ~ Entropy(model,p_out,T_out,z_out)
        h_out ~ Enthalpy(model,p_out,T_out,z_out)

        outport.T ~ T_out
        outport.p ~ p_out
        outport.mdot ~ inport.mdot
        phase_check_out ~ GasPhaseChecker(model,p_out,T_out,z_out)
        outport.z ~ z_out
     ]
    
     compose(ODESystem(eqs, t, vars, para;name=name),inport,outport)
end