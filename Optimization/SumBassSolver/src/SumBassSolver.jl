module SumBassSolver
using DataStructures, Ipopt, JuMP, Roots, SigmoidalProgramming, Optim

export Bass, Bass_derivative, Bass_inverse, Bass_derivative_inverse, Bass_inflection
export sumBassSolver
export sumBassHomogeneousSolver, sumBassDiffStartSolver
export sumSigmoidSolverSB

include("Bass.jl")
include("UdellBoyd.jl")
include("SrivastavaBullo.jl")
include("LocalOpt.jl")
include("HomogeneousSolver.jl")
include("DiffstartSolver.jl")
include("SrivastavaBulloGen.jl")

@doc raw"""
    sumBassSolver(p,q,m,c0,alpha,F,T,n,Pi;alg="UB", localopt=false)

This function maximizes the sum Bass functions and returns the maximizer
and maximum. `alg` is one of UB (Udell and Boyd's algorithm) or
SB (Srivastava and Bullo's algorithm). `localopt` indicates making a local search
after the algorithms terminate.
"""
function sumBassSolver(p,q,m,c0,alpha,F,T,Pi::Real;alg = "SB",localopt=false)
    n = length(p) # Number of sites
    if sum(F) > Pi
        throw(DomainError(Pi,"the problem is unfeasible because the capacity is too small."))
    end
    if n==0
        throw(DomainError(p,"inputs need to be of size at least 1"))
    end
    if length(q) != n || length(m) != n || length(c0) != n || length(alpha) != n || length(F) != n || length(T) != n
        throw(DomainError("p,q,m,c0,alpha,F,T","inputs must be vectors of the same length."))
    end
    if sum(c0 .> m) > 0
        throw(DomainError("c0,m","existing clients needs to be smaller than the market potential."))
    end

    if alg == "UB"
        maximizer,maximum = sumBassUB(p,q,m,c0,alpha,F,T,Pi)
    elseif alg == "SB"
        maximizer,maximum = sumBassSB(p,q,m,c0,alpha,F,T,Pi)
    else
        throw(DomainError(alg,"alg is one of UB (Udell and Boyd's algorithm) or SB (Srivastava and Bullo's algorithm)."))
    end

    if localopt
        maximizer,maximum = sumBassLocal(p,q,m,c0,alpha,F,T,Pi,maximizer)
    end

    return maximizer,maximum
end


end

