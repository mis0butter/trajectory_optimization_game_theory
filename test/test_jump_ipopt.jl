using JuMP 
using Ipopt 
import Random
import Statistics
import Test 
using GLMakie 

using trajectory_optimization_game_theory

## ============================================ ##
# rosenbrock function

function rosenbrock(x, y)
    out = (1.0 - x)^2 + 100.0 * (y - x^2)^2 
  return out
end

# given rosenbrock another method 
rosenbrock(x) = rosenbrock(x[1], x[2]) 

# define another example obj fn 
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# minimize 
x_min = min_ipopt( obj_fn, 3 ) 
x_min = min_ipopt( rosenbrock, 2 ) 

## ============================================ ##

# test solution 
x_min, y_min = min_rosenbrock() 
z_min = rosenbrock( x_min, y_min ) 

x, y = collect(-1.5:0.1:1.5), collect(-1.5:0.1:1.5) 
z    = rosenbrock.(x,y') 

# plot 
fig  = plot_surface(x, y, z) 
fig  = plot_scatter3d(x_min, y_min, z_min, fig ) 
ax = fig.content[1] 
ax.title = "Rosenbrock function" 

display(fig)  


## ============================================ ##
## ============================================ ##

function eq_constraint_1( i, x, h, t )
    out = x[i+1] - x[i] - 0.5 * h * ( sin(t[i+1]) + sin(t[i]) ) 
    return out 
end 

function eq_constraint_2( i, u, h, t ) 
    t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i]
end 

function ineq_constraint_1( i, x, h, t )
    return i, x, h, t 
end 

function obj( u, h, t, alpha ) 
    # sum(
    #     0.5 * h * (u[i+1]^2 + u[i]^2) +
    #     0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
    # )
    i   = 1 
    out = 
        0.5 * h * (u[i+1]^2 + u[i]^2) +
        0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) 

    for i in 2:N 
        out += 
        0.5 * h * (u[i+1]^2 + u[i]^2) +
        0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) 
    end 
    return out 
end 

## ============================================ ##

# set-up 
N = 1000
h = 1 / N
alpha = 350
model = Model(Ipopt.Optimizer)

# def vars constraints 
@variables(model, begin
    -1    <= t[1:(N+1)] <= 1
    -0.05 <= x[1:(N+1)] <= 0.05
    u[1:(N+1)]
end)

# def obj 
@objective(
    model,
    Min,
    obj( u, h, t, alpha ), 
)

# def constraints 
@constraint(
    model,
    [i = 1:N],
    eq_constraint_1( i, x, h, t ) == 0 
)
@constraint(
    model,
    [i = 1:N],
    eq_constraint_2( i, u, h, t ) == 0 
) 

# optimize using JuMP 
JuMP.optimize!(model)

println("""
termination_status = $(termination_status(model))
primal_status      = $(primal_status(model))
objective_value    = $(objective_value(model)) 
""")

# set solution 
x = value.(x) 
t = value.(t) 
u = value.(u) 
    

## ============================================ ##

function example_clnlbeam()
    N = 1000
    h = 1 / N
    alpha = 350
    model = Model(Ipopt.Optimizer)
    @variables(model, begin
        -1 <= t[1:(N+1)] <= 1
        -0.05 <= x[1:(N+1)] <= 0.05
        u[1:(N+1)]
    end)
    @objective(
        model,
        Min,
        sum(
            0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 1:N
        ),
    )
    @constraint(
        model,
        [i = 1:N],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0,
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    objective_value    = $(objective_value(model))
    """)
    
    x = value.(x) 
    t = value.(t) 
    u = value.(u) 
    return x, t, u 
end

example_clnlbeam()




