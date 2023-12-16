using JuMP 
using Ipopt 

## ============================================ ##

" Minimize objective function using JuMP "
function min_ipopt( obj_fn, N ) 

    # set model 
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # define variables 
    @variable(model, x[1:N])
    
    # define objective 
    @objective(model, Min, obj_fn(x)) 
    
    # optimize 
    optimize!(model)

    return value.(x) 
end 

export min_ipopt 


