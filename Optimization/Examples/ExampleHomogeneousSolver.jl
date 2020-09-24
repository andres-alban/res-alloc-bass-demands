using SumBassSolver

p=0.001
q=0.1
m=100
c0=10
alpha = 1
F = 10
T = 50
n = 50
Pi = 600

x1,v1 = sumBassHomogeneousSolver(p,q,m,c0,alpha,F,T,n,Pi)

pp=fill(p,n)
qq=fill(q,n)
mm=fill(m,n)
cc0 = fill(c0,n)
alphas=fill(alpha,n)
FF=fill(F,n)
TT=fill(T,n)

# The general algorithm does find the right solution but takes really long
# The local optimization finishes in a local optimum
x2,v2 = @time sumBassSolver(pp,qq,mm,cc0,alphas,FF,TT,Pi)
x3,v3 = @time sumBassSolver(pp,qq,mm,cc0,alphas,FF,TT,Pi,alg="SB")
x4,v4 = @time sumBassSolver(pp,qq,mm,cc0,alphas,FF,TT,Pi,alg="SB",localopt=true)
