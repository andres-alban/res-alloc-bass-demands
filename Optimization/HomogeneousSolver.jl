# Solver for homogeneous sites
@doc raw"""
    sumBassHomogeneousSolver(p,q,m,c0,alpha,F,T,n,Pi)

This function maximizes the sum of homogeneous Bass functions and returns the maximizer
and maximum. Note that any permutation of the maximizer is also a maximizer.
"""
function sumBassHomogeneousSolver(p::Real,q::Real,m::Real,c0::Real,alpha::Real,F::Real,T::Real,n::Int,Pi::Real)
    if n*F > Pi
        throw(DomainError(Pi,"the problem is unfeasible because the capacity is too small."))
    end
    if n*T <= Pi # Trivial solution
        x = fill(float(T),n)
        v = alpha*n*Bass(T,p,q,m,c0)
        return x,v
    end

    # Case 1
    vstar = 0.0
    xstar = 0
    for i in 1:n # This loop ensures that case does not violate the constraint x<T
        if F+(Pi-F*n)/i < T
            objective = [j * (Bass(F + (Pi-F*n)/j,p,q,m,c0) - Bass(F,p,q,m,c0)) for j in i:n]
            xstar = argmax(objective) + i - 1
            vstar = alpha*(xstar*Bass(F + (Pi-F*n)/xstar,p,q,m,c0) + (n-xstar)*Bass(F,p,q,m,c0))
            break
        end
    end

    # Case 2
    xstar2 = floor(Int,(Pi-n*F)/(T-F))
    vstar2 = alpha*(xstar2*Bass(T,p,q,m,c0) + Bass(Pi-xstar2*T-(n-xstar2-1)*F,p,q,m,c0) +
        (n-xstar2-1)*Bass(F,p,q,m,c0))

    x = Vector{Float64}(undef,n)
    if vstar >= vstar2
        x[1:xstar] .= F + (Pi-F*n)/xstar
        x[(xstar+1):n] .= F
        return x,vstar
    else
        x[1:xstar2] .= T
        x[xstar2+1] = Pi-xstar2*T-(n-xstar2-1)*F
        x[xstar2+2:n] .= F
        return x,vstar2
    end

end
