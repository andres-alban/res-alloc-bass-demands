include("../SumBassSolver.jl")
using .SumBassSolver
using Random

Random.seed!(121)
# Generate random data. For implementation input proper data
p = 0.003*rand(50)
q = 0.3*rand(50)
m = 1000 .+ 500*rand(50)
c0 = m .* rand(50) ./ 2
alpha = rand(50)
F = zeros(50)
T = fill(100.0,50)
Pi = 2500

x1,v1 = @time sumBassSolver(p,q,m,c0,alpha,F,T,Pi)
x2,v2 = @time sumBassSolver(p,q,m,c0,alpha,F,T,Pi,localopt=true)
x3,v3 = @time sumBassSolver(p,q,m,c0,alpha,F,T,Pi,alg="SB")
x4,v4 = @time sumBassSolver(p,q,m,c0,alpha,F,T,Pi,alg="SB",localopt=true)

vv1 = sum(alpha .* Bass.(x1,p,q,m,c0))
vv2 = sum(alpha .* Bass.(x2,p,q,m,c0))
vv3 = sum(alpha .* Bass.(x3,p,q,m,c0))
vv4 = sum(alpha .* Bass.(x4,p,q,m,c0))
sum(abs.(x1-x4))

lambdas = Bass_derivative.(x,p,q,m,c0).*alpha

i = findall(x -> abs(x-lambdas[1]) > 0.01,lambdas) # Sites at the boundary
x[i]

F=fill(10,50)
Pi = 400
x1,v1 = sumBassSolver(p,q,m,c0,alpha,F,T,Pi,localopt=true)

Pi = 600
x1,v1 = sumBassSolver(p,q[1],m,c0,alpha,F,T,Pi,localopt=true)
