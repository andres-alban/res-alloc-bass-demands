using SumBassSolver
using Random

f = [x-> 1/(1 + exp(-x)), x->2/(1 + exp(-x/2))]
df = [x->exp(-x)/(1 + exp(-x))^2, x-> exp(-x/2)/(1 + exp(-x/2))^2]
z = [0,0]
F = [-1,-1]
T = [2,2]
Pi = 2

sumSigmoidSolverSB(f,df,z,F,T,Pi)

Random.seed!(121)
# Generate random data. For implementation input proper data
n=50
p = 0.003*rand(n)
q = 0.3*rand(n)
m = 1000 .+ 500*rand(n)
c0 = m .* rand(n) ./ 2
alpha = rand(n)
F = zeros(n)
T = fill(100.0,n)
Pi = 2500

f = Vector{Function}(undef,n)
df = Vector{Function}(undef,n)
z = Vector{Float64}(undef,n)
df_inv = Vector{Function}(undef,n)
for i in 1:n
    f[i] = t -> alpha[i] * Bass(t,p[i],q[i],m[i],c0[i])
    df[i] = t -> alpha[i] * Bass_derivative(t,p[i],q[i],m[i],c0[i])
    z[i] = Bass_inflection(p[i],q[i],m[i],c0[i])
    df_inv[i] = t -> Bass_derivative_inverse(t/alpha[i],p[i],q[i],m[i],c0[i])[1]
end

x1,v1 = sumSigmoidSolverSB(f,df,z,F,T,Pi)
x2,v2 = sumSigmoidSolverSB(f,df,z,F,T,Pi,df_inv)
x3,v3 = @time sumBassSolver(p,q,m,c0,alpha,F,T,Pi,alg="SB",localopt=true)

