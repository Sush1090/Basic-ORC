using Clapeyron, Plots
fluids =["pentane","butane"]
model =  cPR(fluids,idealmodel = ReidIdeal)
x = range(0,0.99,length=100)
X = Clapeyron.FractionVector.(x)
P = 101325.0*2


T_bubble = zeros(size(X,1))
T_dew = zeros(size(X,1))

for i in eachindex(X)
    T_bubble[i] = bubble_temperature(model,P,X[i])[1]
    T_dew[i] = dew_temperature(model,P,X[i])[1]
end 

fig1 = plot(x,T_dew,label = "Dew temperature")
plot!(x,T_bubble,label = "Bubble temperature")
xlabel!("Fraction 1st fluid")
ylabel!("Temperature (K)")
title!("$(fluids)")
