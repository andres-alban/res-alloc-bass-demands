using SumBassSolver
using Plots

p=fill(0.001,2)
q=[0.1,0.11]
m=fill(100,2)
c0=fill(0,2)
alpha = fill(1,2)
F = fill(0,2)
T = fill(100,2)
Pi=90

sumBassSolver(p,q,m,c0,alpha,F,T,Pi)
sumBassSolver(p,q,m,c0,alpha,F,T,Pi,localopt=true)
sumBassSolver(p,q,m,c0,alpha,F,T,Pi,alg="SB")
sumBassSolver(p,q,m,c0,alpha,F,T,Pi,alg="SB",localopt=true)

Pis = 50:100
v = Vector{Float64}(undef,length(Pis))
x = Array{Float64}(undef,2,length(Pis))

@time for i in 1:length(Pis)
    x[:,i],v[i] = sumBassSolver(p,q,m,c0,alpha,F,T,Pis[i],localopt=true)
end

pl1 = plot(Pis,v);
pl2 = plot(Pis,x[1,:],label = "plot1");
plot!(Pis,x[2,:],label = "plot2");
plot!(Pis,fill(Bass_inflection(p[1],q[1],m[1],c0[1]),length(Pis)),label="inflection1");
plot!(Pis,fill(Bass_inflection(p[2],q[2],m[2],c0[2]),length(Pis)),label="inflection2");
pl2
