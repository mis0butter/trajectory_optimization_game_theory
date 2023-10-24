# test_propagator
using trajectory_optimization_game_theory 
using Plots 

mu = 398600.4415
r = 6378.0
x0 = kep2cart([r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0], mu)
println(x0)
t = (0.0, 90*60.0)
prob = propagate_2Body(x0, t, mu)

xVals = []
yVals = []
zVals = []
for i in 1:length(prob.u)
    push!(xVals, prob.u[i][1])
    push!(yVals, prob.u[i][2])
    push!(zVals, prob.u[i][3])
end

scatter(xVals, yVals, zVals)



