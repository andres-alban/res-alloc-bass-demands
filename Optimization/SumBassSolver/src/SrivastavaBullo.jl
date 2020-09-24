# using Optim

@doc raw"""
    sumBassSB(p,q,m,c0,alpha,F,T,Pi)

This function implements Srivastava and Bullo's algorithm to maximize a sum of
Bass functions and returns the maximizer and maximum.
"""
function sumBassSB(p,q,m,c0,alpha,F,T,Pi::Real)
    # max FLP
    ub = maximum(Bass_inflection.(p,q,m,c0))
    # res = optimize(lambda -> -FLP(lambda)[2], 0.0, ub)
    # There are some issues with the optimizer because the objective is discontinuous with the maximizer at one of the discontinuities
    # return FLP(Optim.minimizer(res))
    # Instead make a grid search
    lambdas = collect(range(0,ub,length=1000))
    FLPs = [FLP(lambdas[i],p,q,m,c0,alpha,F,T,Pi)[2] for i in 1:length(lambdas)]
    lambdastar = lambdas[partialsortperm(FLPs,1:50,rev=true)]
    lambdas .= range(minimum(lambdastar),maximum(lambdastar),length=1000)
    FLPs .= [FLP(lambdas[i],p,q,m,c0,alpha,F,T,Pi)[2] for i in 1:length(lambdas)]
    # maximizer = FLP(lambdas[argmax(FLPs)],p,q,m,c0,alpha,F,T,Pi)[1]
    # maximum = sum(alpha .* Bass.(maximizer,p,q,m,c0))
    return FLP(lambdas[argmax(FLPs)],p,q,m,c0,alpha,F,T,Pi)
end

@doc raw"""
    FLP(lambda,p,q,m,c0,alpha,F,T,Pi)

This function evaluates the linear program relaxation of the maximization of the
sum of Bass functions.
"""
function FLP(lambda,p,q,m,c0,alpha,F,T,Pi::Real)
    n = length(p)
    Pi_eff = Pi - sum(F)
    f = Vector{Function}(undef,n)
    # df = Vector{Function}(undef,n)
    # df_inv = Vector{Function}(undef,n)
    # z = Vector{Float64}(undef,n)
    for i in 1:n
        f[i] = t -> alpha[i]*Bass(t,p[i],q[i],m[i],c0[i])
        # df_inv[i] = lambda -> Bass_derivative_inverse(lambda/alpha[i],p[i],q[i],m[i],c0[i])[1]
    end

    df_invs = [Bass_derivative_inverse(lambda/alpha[i],p[i],q[i],m[i],c0[i])[1] for i in 1:n]
    df_invs .= [F[i] < df_invs[i] < T[i] ? df_invs[i] : 0.0 for i in 1:n]

    LPindex = Vector{Float64}(undef,n)
    xLP = Vector{Float64}(undef,n)
    gain_rate = Vector{Float64}(undef,3)
    for i in 1:n
        gain_rate .= [df_invs[i] > 0 ? f[i](df_invs[i])/df_invs[i] : 0, F[i]>0 ? f[i](F[i])/F[i] : 0, T[i] > 0 ? f[i](T[i])/T[i] : 0]
        xLPind = argmax(gain_rate)
        xLP[i] = [df_invs[i],F[i],T[i]][xLPind]
        LPindex[i] = gain_rate[xLPind]
    end
    LPperm = sortperm(LPindex,rev = true)
    cum_capacity = 0.0
    k = n
    for i in 1:n
        cum_capacity += xLP[LPperm[i]] - F[LPperm[i]]
        if cum_capacity >= Pi_eff
            k = i
            break
        end
    end
    if k < n
        xLP[LPperm[(k+1):n]] .= F[(k+1):n]
    end
    xLP[LPperm[k]] -= max(cum_capacity - Pi_eff,0.0) # if the capacity is not
    return xLP, sum(alpha[i].*Bass(xLP[i],p[i],q[i],m[i],c0[i]) for i in 1:n)
end
