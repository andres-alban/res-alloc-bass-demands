# Solver for different starting points
using Roots

@doc raw"""
    sumBassDiffStartSolver(p,q,m,c0,alpha,F,T,Pi)

This function maximizes the sum of Bass functions with same p, q, and m, but
heterogeneous 'starting points' c0. The function returns the maximizer
and maximum.
"""
function sumBassDiffStartSolver(p::Real,q::Real,m::Real,c0,alpha::Real,F::Real,T::Real,Pi::Real)
    n=length(c0)
    if n*F > Pi
        throw(DomainError(Pi,"the problem is unfeasible because the capacity is too small."))
    end
    x0 = Bass_inverse.(c0,p,q,m,0)
    z0 = Bass_inflection(p,q,m,0)
    C1 = findall(x0 .+ F .>= z0)
    C2 = findall(x0 .+ F .< z0)
    ind2 = sortperm(x0[C2],rev=true)
    ind1 = sortperm(x0[C1])
    ind = vcat(C2[ind2],C1[ind1])

    function P(l,lambda,p,q,m,F,T,x0,lenC2)
        xhat = Bass_derivative_inverse(lambda,p,q,m,0)[1]
        result = sum(min.(xhat .- x0[1:l],T)) + F*(lenC2 - l)
        if lenC2<n
            result += sum(min.(max.(xhat .- x0[(lenC2+1):end],F),T))
        end
        return result
    end

    lenC2 = length(C2)
    lambda_ub = Bass_derivative(z0,p,q,m,0)
    lambdas = Vector{Float64}(undef,lenC2)
    Psi = Vector{Float64}(undef,lenC2+1)
    for l in 1:lenC2
        try
            lambdas[l] = find_zero(y -> P(l,y,p,q,m,F,T,x0[ind],lenC2) - Pi,(0,lambda_ub))
            xhat = Bass_derivative_inverse(lambdas[l],p,q,m,0)[1]
            Psi[l] = alpha*sum(Bass(min(xhat - x0[ind[i]],T),p,q,m,c0[ind[i]]) for i in 1:l)
            if l<lenC2
                Psi[l] += alpha*sum(Bass(F,p,q,m,c0[ind[i]]) for i in (l+1):lenC2)
            end
            if lenC2<n
                Psi[l] += alpha*sum([Bass(min(max(xhat - x0[ind[i]],F),T),p,q,m,c0[ind[i]]) for i in (lenC2+1):n])
            end
        catch
            lambdas[l] = 0
            Psi[l] = 0
        end
    end

    Psi[lenC2+1] = alpha*Bass(min(Pi,T),p,q,m,c0[ind[1]]) + alpha*sum([Bass(F,p,q,m,c0[ind[i]]) for i in 2:n])
    lstar = argmax(Psi)
    x = Vector{Float64}(undef,n)
    if lstar == lenC2+1
        x .= F
        x[ind[1]] = min(Pi,T)
    else
        xhat = Bass_derivative_inverse(lambdas[lstar],p,q,m,0)[1]
        x[ind[1:lstar]] = [min(xhat - x0[ind[i]],T) for i in 1:lstar]
        x[ind[(lstar+1):lenC2]] .= F
        if lenC2<n
            x[ind[(lenC2+1):n]] = [min(max(xhat - x0[ind[i]],F),T) for i in (lenC2+1):n]
        end
    end
    return x, sum(alpha .* Bass.(x,p,q,m,c0))
end
