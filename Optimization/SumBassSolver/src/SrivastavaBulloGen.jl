@doc raw"""
    sumSigmoidSolverSB(f,z,F,T,Pi)

This function implements Srivastava and Bullo's algorithm to maximize a sum of
general sigoidal functions and returns the maximizer and maximum. `f` is an array
of sigmoidal functions, `df` its derivatives and `z` is the array of inflection points.
"""
function sumSigmoidSolverSB(f,df,z,F,T,Pi,df_inv=nothing)
    ub = maximum([-optimize(x->-df[i](x),F[i],T[i]).minimum for i in 1:length(df)])
    lambdas = collect(range(0.,ub,length=1000))
    FLPs = [FLP_sigmoid(lambdas[i],f,df,df_inv,z,F,T,Pi)[2] for i in 1:length(lambdas)]
    lambdastar = lambdas[partialsortperm(FLPs,1:50,rev=true)]
    lambdas .= range(minimum(lambdastar),maximum(lambdastar),length=1000)
    FLPs .= [FLP_sigmoid(lambdas[i],f,df,df_inv,z,F,T,Pi)[2] for i in 1:length(lambdas)]
    return FLP_sigmoid(lambdas[argmax(FLPs)],f,df,df_inv,z,F,T,Pi)
end

@doc raw"""
    FLP_sigmoid(lambda,f,df,df_inv,z,F,T,Pi)

This function evaluates the linear program relaxation of the maximization of the
sum of sigmoidal functions.
"""
function FLP_sigmoid(lambda,f,df,df_inv,z,F,T,Pi::Real)
    n = length(f)
    Pi_eff = Pi - sum(F)

    df_invs = Vector{Float64}(undef,n)
    if isnothing(df_inv)
        for i in 1:n
            try
                df_invs[i] = find_zero(x-> df[i](x)-lambda,(z[i],T[i]))
            catch
                df_invs[i] = 0.
            end
        end
    else
        df_invs = [df_inv[i](lambda) for i in 1:n]
        df_invs .= [F[i] < df_invs[i] < T[i] ? df_invs[i] : 0.0 for i in 1:n]
    end

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
    xLP[LPperm[k]] -= max(cum_capacity - Pi_eff,0.0) # if the capacity is not reached
    return xLP, sum(f[i](xLP[i]) for i in 1:n)
end