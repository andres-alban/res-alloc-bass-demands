using JuMP, Ipopt

@doc raw"""
    sumBassLocal(p,q,m,c0,alpha,F,T,Pi,xstart)

This function locally maximizes a sum of Bass functions with a hot start at `xstart`
and returns the maximizer and maximum. It uses JuMP's interface with the IpOpt solver.
"""
function sumBassLocal(p,q,m,c0,alpha,F,T,Pi::Number,xstart)
    n = length(p) # Number of sites
    model = Model(optimizer_with_attributes(Ipopt.Optimizer,"print_level" => 0))
    @variable(model,x[1:n])
    set_start_value.(x,xstart)
    @constraint(model, F .<= x .<= T)
    # for i in 1:n
    #     set_start_value(x[i],node.x[i]) # Starting point at the maximizer found by SigmidalProgramming.jl
    #     @constraint(model, F[i] <= x[i] <= T[i])
    # end
    C = ones(1,n)
    @constraint(model,sum(C[i]*x[i] for i=1:n) == Pi)
    SumBass(x...) = sum(alpha[i]*Bass(x[i],p[i],q[i],m[i],c0[i]) for i=1:length(x))
    function SumBass_derivative(g,x...)
        for i in 1:length(x)
            g[i] = alpha[i]*Bass_derivative(x[i],p[i],q[i],m[i],c0[i])
            nothing
        end
    end
    register(model,:SumBass,n,SumBass,SumBass_derivative)
    @NLobjective(model,Max,SumBass(x...))
    # register(model,:Bass,5,Bass,autodiff=true)
    # @NLobjective(model,Max,sum(alpha[i]*Bass(x[i],p[i],q[i],m[i],c0[i]) for i=1:n))
    optimize!(model)
    maximizer = value.(x)
    maximum = objective_value(model)
    # if termination_status(model) == MOI.LOCALLY_SOLVED
    #     return value.(x), objective_value(model)
    # else
    #     error("The model was not solved correctly.")
    # end
    return maximizer,maximum
end
