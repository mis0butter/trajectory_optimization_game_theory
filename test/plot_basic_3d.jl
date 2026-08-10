using CairoMakie 
using LinearAlgebra: norm 

## ==================================================================== 

r   = 6378.0
xyz = [ zeros(3) for i in 1:3 ] 
uvw = r .* [ [1,0,0] , [0,1,0] , [0,0,1] ] 

width = norm(uvw[1])/100 
color = :black 

ps = [ Point3f(x,y,z) for (x,y,z) in xyz ] 
ns = [ Vec3f(u,v,w) for (u,v,w) in uvw ] 

fig = Figure() 
ax = Axis3(fig[1,1])
arrows!(ps, ns, fxaa = true, align = :origin)
c = 1.5
# xlims!(ax, -c*r, c*r) 
# ylims!(ax, -c*r, c*r) 
# zlims!(ax, -c*r, c*r) 

fig 

## ================================== 

# r   = 6378.0 
r   = 1
xyz = [ zeros(3) for i in 1:3 ] 
uvw = r .* [ [1,0,0] , [0,1,0] , [0,0,1] ] 

ps = [ Point3f(x,y,z) for (x,y,z) in xyz ] 
ns = [ Vec3f(u,v,w) for (u,v,w) in uvw ] 

width = r/100

fig = Figure() 
ax = Axis3(fig[1,1])
arrows!(
    ps, ns, fxaa = true, 
    linewidth = width, 
    align = :origin, 
) 
# :origin, :head, :lineend, :tailend, :headstart or :center  

c = 3.0 
xlims!(ax, -c*r, c*r) 
ylims!(ax, -c*r, c*r) 
zlims!(ax, -c*r, c*r) 

ax.azimuth = rand() * 2pi 
ax.elevation = rand() * 2pi 

fig 

## ==================================================================== 

ps = [Point3f(x, y, z) for x in -5:2:5 for y in -5:2:5 for z in -5:2:5]
ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)

## ==================================================================== 

ps = [ Point3f(0,0,0) ] 
ns = [ Vec3f(1,0,0) ]  

arrows(ps, ns, fxaa = true, align = :origin)
