# test_propagator
using trajectory_optimization_game_theory 
using Plots 

## ============================================ ##

mu = 398600.4415
r = 6378.0
kep0_p1 = [r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
kep0_p2 = [42164, 0.0, 0.0, 0.0, 0.0, 0.0]
t = (0.0, orbitPeriod(kep0_p2, mu))
prop1 = propagate_2Body(kep2cart(kep0_p1, mu), t, mu)
prop2 = propagate_2Body(kep2cart(kep0_p2, mu), t, mu)

xVals_p1 = []
yVals_p1 = []
zVals_p1 = []
xVals_p2 = []
yVals_p2 = []
zVals_p2 = []
for i in 1:length(prop1.u)
    push!(xVals_p1, prop1.u[i][1])
    push!(yVals_p1, prop1.u[i][2])
    push!(zVals_p1, prop1.u[i][3])
end
for i in 1:length(prop2.u)
    push!(xVals_p2, prop2.u[i][1])
    push!(yVals_p2, prop2.u[i][2])
    push!(zVals_p2, prop2.u[i][3])
end


scatter(xVals_p1, yVals_p1, zVals_p1, label="p1")
scatter!(xVals_p2, yVals_p2, zVals_p2, label="p2")

