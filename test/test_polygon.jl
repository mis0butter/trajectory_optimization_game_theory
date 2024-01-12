using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

## ============================================ ##
# init params 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180 ]
rv_0_P = kep2cart(kep0_P, mu) 
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 120.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 

# tof for pursuer to catch up to evader 
tof = 2000 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

## ============================================ ##
# ... it's just vector addition to find the vertices of the polygon 

function polygon_vertices( rv_f_E ) 

    # center of polygon 
    r_f = rv_f_E[1:3]  ;    axis_2 = -r_f / norm(r_f)
    v_f = rv_f_E[4:6]  ;    axis_1 = v_f / norm(v_f) 

    # define vector normal to orbit plane 
    axis_3 = cross( axis_1, axis_2 )  ; axis_3 = axis_3 / norm(axis_3)     

    # top vertex: move up from r_f along axis 3 
    r_top = r_f + axis_3 * r/10 

    # top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * r/10 + sind(60) * axis_2 * r/10 
    r_topin  = r_f + vec

    # bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * r/10 + sind(60) * axis_2 * r/10 
    r_botin  = r_f + vec 

    # bottom vertex: move down from r_f along axis 3 
    r_bot = r_f - axis_3 * r/10 

    # bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * r/10 - sind(60) * axis_2 * r/10 
    r_botout = r_f + vec 

    # top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
    vec       = cosd(60) * axis_3 * r/10 - sind(60) * axis_2 * r/10 
    r_topout  = r_f + vec 

    out = ( top = r_top, topin = r_topin, botin = r_botin, bot = r_bot, botout = r_botout, topout = r_topout ) 
    return out 
end 

## ============================================ ##


# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 

# let's propagate the evader SC 
rv_f_E = rv_E[end,:] 

# center of polygon 
r_f = rv_f_E[1:3]  ; axis_1 = v_f / norm(v_f)   
v_f = rv_f_E[4:6]  ; axis_2 = -r_f / norm(r_f)

# define vector normal to orbit plane 
axis_3 = cross( axis_1, axis_2 )  ; axis_3 = axis_3 / norm(axis_3) 

# ok, let's plot this so that it all looks right 
fig = plot_vector3d( [ r_f ] , [ axis_1 * r ] , fig, r/100, :black, "1" ) 
fig = plot_vector3d( [ r_f ] , [ axis_2 * r ] , fig, r/100, :black, "2" ) 
fig = plot_vector3d( [ r_f ] , [ axis_3 * r ] , fig, r/100, :black, "3" ) 

# vertices of polygon along axis 2-3 plane 
# ok, let's define the distance of vertices of polygon from center: how about r / 100 ? 

# top vertex: move up from r_f along axis 3 
r_top = r_f + axis_3 * r/10 
# fig   = plot_scatter3d( r_top[1], r_top[2], r_top[3], fig, :circle ) 
l_topin = [ r_top[1], r_top[2], r_top[3] ] 

# top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
vec      = cosd(60) * axis_3 * r/10 + sind(60) * axis_2 * r/10 
r_topin  = r_f + vec
# fig      = plot_scatter3d( r_topin[1], r_topin[2], r_topin[3], fig, :circle )  

mat = [ r_top' ; r_topin' ] 
fig = plot_line3d( mat, fig ) 

# bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
vec      = - cosd(60) * axis_3 * r/10 + sind(60) * axis_2 * r/10 
r_botin  = r_f + vec 
# fig      = plot_scatter3d( r_botin[1], r_botin[2], r_botin[3], fig, :circle ) 

mat = [ r_topin' ; r_botin' ] 
fig = plot_line3d( mat, fig ) 

# bottom vertex: move down from r_f along axis 3 
r_bot = r_f - axis_3 * r/10 
# fig   = plot_scatter3d( r_bot[1], r_bot[2], r_bot[3], fig, :circle ) 

mat = [ r_botin' ; r_bot' ] 
fig = plot_line3d( mat, fig ) 

# bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
vec      = - cosd(60) * axis_3 * r/10 - sind(60) * axis_2 * r/10 
r_botout = r_f + vec 
# fig      = plot_scatter3d( r_botout[1], r_botout[2], r_botout[3], fig, :circle ) 

mat = [ r_bot' ; r_botout' ] 
fig = plot_line3d( mat, fig ) 

# top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
vec       = cosd(60) * axis_3 * r/10 - sind(60) * axis_2 * r/10 
r_topout  = r_f + vec 
# fig       = plot_scatter3d( r_topout[1], r_topout[2], r_topout[3], fig, :circle ) 

mat = [ r_botout' ; r_topout' ] 
fig = plot_line3d( mat, fig ) 

mat = [ r_topout' ; r_top' ] 
fig = plot_line3d( mat, fig ) 



