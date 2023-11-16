using JuMP 
using Ipopt 
import Random
import Statistics
import Test 

## ============================================ ##

function eq_constraint_1( i, x, h, t )
    
    out = x[i+1] - x[i] - 0.5 * h * ( sin(t[i+1]) + sin(t[i]) ) 

    return out 
end 

function ineq_constraint_1( i, x, h, t )
    
    return i, x, h, t 

end 

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
        # x[i+1] - x[i] - 0.5 * h * ( sin(t[i+1]) + sin(t[i] ) ) == 0,
        eq_constraint_1( i, x, h, t ) == 0 
    )
    @constraint(
        model,
        [i = 1:N],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0,
    )
    # @constraint(
    #     model,
    #     [i = 1:N],
    #     ineq_constraint_1( i, x, h, t ) >= 0,
    # )
    out = JuMP.optimize!(model)
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

## ============================================ ##

x, t, u = example_clnlbeam()




