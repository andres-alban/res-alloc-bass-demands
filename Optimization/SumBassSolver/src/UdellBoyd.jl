
@doc raw"""
    sumBassUB(p,q,m,c0,alpha,F,T,Pi)

This function implements Udell an Boyd's algorithm to maximize a sum of Bass functions
and returns the maximizer and maximum.
"""
function sumBassUB(p,q,m,c0,alpha,F,T,Pi::Real)
    n = length(p)
    # Generate model for solution
    f = Vector{Function}(undef,n)
    df = Vector{Function}(undef,n)
    z = Vector{Float64}(undef,n)
    for i in 1:n
        f[i] = t -> alpha[i]*Bass(t,p[i],q[i],m[i],c0[i])
        df[i] = t -> alpha[i]*Bass_derivative(t,p[i],q[i],m[i],c0[i])
        z[i] = Bass_inflection(p[i],q[i],m[i],c0[i])
    end

    C = ones(1,n)
    d = [Pi]

    A = zeros(0,n)
    b = zeros(0)

    problem = LinearSP(f, df, z, A, b, C, d)

    pq, bestnodes, lbs, ubs = solve_sp(F, T, problem,TOL=1e-6)

    node = dequeue!(pq)

    maximizer = Vector{Float64}(undef,n)
    maximizer .= node.x
    maximum = sum(alpha[i]*Bass(maximizer[i],p[i],q[i],m[i],c0[i]) for i in 1:n)
    return maximizer, maximum
end

    # Maximizer
    # node.x
    # Lower bound of maximum
    # node.lb
    # Upper bound of maximum
    # node.ub
