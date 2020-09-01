
    # To allow for p=0 when c0>0, theta is defined different compared to the paper
    # The problem is that if p=0, we would divide by zero, even when p would cancel out in the final expression
    # if theta' is the definition of the paper, then theta = theta'/p

function Bass(t,p,q,m,c0)::Float64
    theta = (m - c0)/(q*c0 + p*m)
    m* (1 - p*theta*exp(-(p + q)*t)) / (1 + q*theta*exp(-(p + q)*t))
end

function Bass_derivative(t,p,q,m,c0)
    if isinf(t)
        return 0.0
    end
    exp((p+q)*t)*m*(-c0 + m)*(p+q)^2*(m*p+c0*q)/((-c0+m)*q + exp((p+q)*t)*(m*p+c0*q))^2
end

function Bass_inflection(p,q,m,c0)
    theta = (m - c0)/(q*c0 + p*m)
    t = max(-1/(p + q) * log(1/q/theta),0.0)
end

function Bass_inverse(Y,p,q,m,c0)
    theta = (m - c0)/(q*c0 + p*m)
    y = Y/m
    t = log(theta*(p+y*q)/(1-y))/(p+q)
end

function Bass_derivative_inverse(lambda,p,q,m,c0)
    if lambda == 0.0
        return Inf,0.0
    end
    infl = Bass_inflection(p,q,m,c0)
    if Bass_derivative(infl,p,q,m,c0) <= lambda
        t1 = infl
        t2 = infl
        return t1,t2
    else
        theta = (m - c0)/(q*c0 + p*m)
        y=lambda/m
        a = q^2*y
        b = 2*q*y - (p+q)^2
        c = y
        x1 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
        x2 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
        t1 = -log(x1/theta)/(p+q)
        t2 = max(-log(x2/theta)/(p+q),0.0)
        return t1, t2
    end
end
