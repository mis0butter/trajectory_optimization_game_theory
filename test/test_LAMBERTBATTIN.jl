using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 
# define IC, target state, and lambert solve 

r1 = [20.0e6, 20.0e6, 0] ;   # [m]
r2 = [-20.0e6, 10.0e6, 0] ;  # [m]
tof = 1.0 * 86400;
mu = 3.986004418e14 ;        # m3/s2

dm      = "pro" 
Dtsec   = tof 

## ============================================ ##
# solve and propagate lambert orbit 

v1, v2 = lambertbattin( r1, r2, mu, dm, tof ) 

x0 = [ r1; v1 ] 
t_lambert, x_lambert = propagate_2Body( x0, tof, mu, 1.0 ) 
x_lambert = mapreduce( permutedims, vcat, x_lambert ) 

fig = plot_orbit( x_lambert )
fig = plot_scatter3d( r2[1], r2[2], r2[3], fig ) 



## ============================================ ##
## ============================================ ##
# test negative propagation 

# xf = x_lambert[end,:] 
vf = v2 ; vf[end] = 100
# xf = [r2; v2]
xf = [r2; vf] 
t_reverse, x_reverse = propagate_2Body( xf, -tof, mu, 1.0 ) 
x_reverse = mapreduce( permutedims, vcat, x_reverse ) 

x_E_0  = x_reverse[end,:] 
t_E, x_E = propagate_2Body( x_E_0, tof, mu, 1.0 ) 
x_E = mapreduce( permutedims, vcat, x_E ) 

## ============================================ ##

fig = plot_orbit( rv_lambert )
fig = plot_orbit( x_reverse, fig )  
fig = plot_orbit( x_E, fig )  

# plot target 
scatter!( r2[1], r2[2], r2[3]; marker = :circle, markersize = 10, color = :black ) 
text!( r2[1], r2[2], r2[3]; text = "target", color = :gray, offset = (0,-10), align = (:center, :bottom) )




## ============================================ ##
# test seebatt and seebattk 

x        = 0.265264639249882 
x_matlab = 5.32072850287158 
x_julia  = seebatt(x) 
println( "err = ", x_matlab - x_julia ) 

x        = 2.13801399029609 
x_matlab = 6.97629585037618 
x_julia  = seebatt(x )
println( "err = ", x_matlab - x_julia ) 

u        = 0.694215270078835 
u_matlab = 0.302248110562253 
u_julia  = seebattk(u) 
println( "err = ", u_matlab - u_julia ) 

u        = 0.42135922481563 
u_matlab = 0.313748053815816 
u_julia  = seebattk(u) 
println( "err = ", u_matlab - u_julia ) 

## ============================================ ##
# test lambert 

tof     = 1.0 * 86400;
mu      = 3.986004418e14 ;        # m3/s2
dm      = "pro" 
Dtsec   = tof 

r1 = [20.0e6, 20.0e6, 0] ;   # [m]
r2 = [-20.0e6, 10.0e6, 0] ;  # [m]

v1_matlab = [   778.548926356586
                4306.37485999949
                0 ]
v2_matlab = [   2230.87531343386
                -4643.26359035983
                0 ]
v1, v2 = lambertbattin( r1, r2, mu, dm, tof ) 

## ============================================ ##







