using Random
include("../SumBassSolver.jl")
using .SumBassSolver

Random.seed!(121)
n=50
p=0.001
q=0.1
m=100
c0=rand(Float64,n).*m./2
alpha = 1
F = 10
T = 100
Pi = 2500

x,v = sumBassDiffStartSolver(p,q,m,c0,alpha,F,T,Pi)

println(alpha .* Bass_derivative.(x,p,q,m,c0))

pp=fill(p,n)
qq=fill(q,n)
mm=fill(m,n)
alphas=fill(alpha,n)
FF=fill(F,n)
TT=fill(T,n)

x1,v1 = sumBassSolver(pp,qq,mm,c0,alphas,FF,TT,Pi)
x2,v2 = sumBassSolver(pp,qq,mm,c0,alphas,FF,TT,Pi,alg="SB")
x3,v3 = sumBassSolver(pp,qq,mm,c0,alphas,FF,TT,Pi,alg="SB",localopt=true)
