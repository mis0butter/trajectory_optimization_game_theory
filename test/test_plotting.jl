using GLMakie 

## ============================================ ##
# https://juliadatascience.io/glmakie

function peaks(; n=49)
    x = range(-3, 3, n)
    y = range(-3, 3, n)
    a = 3 * (1 .- x') .^ 2 .* exp.(-(x' .^ 2) .- (y .+ 1) .^ 2)
    b = 10 * (x' / 5 .- x' .^ 3 .- y .^ 5) .* exp.(-x' .^ 2 .- y .^ 2)
    c = 1 / 3 * exp.(-(x' .+ 1) .^ 2 .- y .^ 2)
    return (x, y, a .- b .- c)
end

# peaks plot 

x, y, z = peaks()
x2, y2, z2 = peaks(; n=15)
fig = Figure(resolution=(1200, 400))
axs = [Axis3(fig[1, i]; aspect=(1, 1, 1)) for i = 1:3]
hm = GLMakie.surface!(axs[1], x, y, z)
GLMakie.wireframe!(axs[2], x2, y2, z2)
GLMakie.contour3d!(axs[3], x, y, z; levels=20)
Colorbar(fig[1, 4], hm, height=Relative(0.5))
fig 

## ============================================ ##

fig = plot_surface( x, y, z ) 
fig = plot_contour3d( x2, y2, z2 ) 

## ============================================ ##

xyz = [0 0 0] 
uvw = [1 2 3] 

fig = plot_vector3d( xyz, uvw ) 

## ============================================ ##

xyz = [ 1 2 3 ; 4 5 6 ] 

l = size(xyz, 1) 

xyz_vec = [] 
for i = 1 : l 
    push!( xyz_vec, xyz[i,:] ) 
end 

xyz_vec = [ xyz[i,:] for i = 1 : size(xyz, 1) ] 

## ============================================ ##

xyz = [ (1,-2,-3), (1,2,3) ] 
ps  = [ Point3f(x,y,z) for (x,y,z) in xyz ] 

uvw = [ (1,2,3), (-1,-2,-3) ] 
ns  = [ Vec3f(x,y,z) for (x,y,z) in uvw ] 

xyz = [ zeros(3) for i in 1:3 ] 
uvw = [ [1,0,0] , [0,1,0] , [0,0,1] ] 


fig = arrows(
    ps, ns, fxaa=true, # turn on anti-aliasing
    linecolor = :gray, arrowcolor = :black,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
    align = :center, axis=(type=Axis3,)
) 
Auto() 

ps = [ Point3f( 1,-2,3 ) ] 
ns = [ Vec3f( -1,2,3 ) ] 

arrows!(  
    ps, ns, fxaa=true, # turn on anti-aliasing
    linecolor = :red, arrowcolor = :red,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
    align = :center, 
)

fig = plot_vector3d( xyz, uvw ) 

# GLMakie.scatter!( 0,0,0; marker = :circle, markersize = 40, color = :red ) 
# fig = plot_scatter3d( 0,    0,   0, fig ) 
# fig = plot_scatter3d( 10, -20, -30, fig, :circle, :green ) 
# fig = plot_scatter3d( 10,  20,  30, fig, :xcross, :blue )  

## ============================================ ##
 
xyz = [ zeros(3) for i in 1:3 ] 
uvw = [ [2,0,0] , [0,1,0] , [0,0,1] ] 

fig = plot_vector3d( [ xyz[1] ] , [ uvw[1] ] ) 
fig = plot_vector3d( [ xyz[2] ] , [ uvw[2] ], fig, :red ) 
fig = plot_vector3d( [ xyz[3] ] , [ uvw[3] ], fig, :green ) 
# fig


    