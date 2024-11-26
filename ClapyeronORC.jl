using ModelingToolkit, DifferentialEquations, Clapeyron
using ModelingToolkit: t_nounits as t, D_nounits as D


fluids =["pentane","butane"]
global model = cPR(fluids,idealmodel = ReidIdeal)
include("Utils.jl")
include("Process.jl")

composition = [0.5,0.5] # by mass percentage
@named src = MassSource(composition = composition)
@named pump = Pump()
@named evap = Evaporator()
@named expander = Expander()
eqs = [
        connect(src.port,pump.inport)
        connect(pump.outport,evap.inport)
        connect(evap.outport,expander.inport)
    ]
systems=[src,pump,evap,expander] # Define system

@named orc = System(eqs, t, systems=systems)
sys = structural_simplify(orc)
u0 = []
para = [sys.src.p_source => 101325*2, sys.src.T_source => 288,sys.src.mdot_source => 50,
        sys.pump.πc =>4,  sys.expander.πc => sys.pump.πc,
        sys.pump.η =>0.5, sys.expander.η=>0.7]
prob = SteadyStateProblem(sys,u0,para)
sol = solve(prob)


power = sol[expander.h_out][1] - sol[expander.h_in][1]
heat =  sol[evap.h_out][1] - sol[evap.h_in][1]
pump_work = sol[pump.h_out][1] - sol[pump.h_in][1]
η =  (power + pump_work)/heat
@show abs(η)