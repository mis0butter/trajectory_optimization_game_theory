using GLMakie 
GLMakie.activate!()

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



# fig


    