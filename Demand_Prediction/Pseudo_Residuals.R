Power = function(a,b){
  a^b
}

penalty_deriv = function(t,rho,eps,k,m=NULL){
  if(is.null(m)){
    m = length(t)
  }
  -sum(rho*ifelse(t<0,
                  0,ifelse(t < (eps/(m*rho))^(1/k),
                           (3*k*m^2*rho^2*t^(3*k-1))/(eps^2) - (4*k*m^3*rho^3*t^(4*k-1))/(2*eps^3),k*t^(k-1))))
}

residual_1 = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,y,v,penalty_params,m){
  LF = -2*(b1 + F1*(2*b2 + 3*b3*F1))*mf*((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))
                                         /(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))
                                         - (Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                 F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) - ((b0 + F1*(b1 + F1*(b2 + b3*F1)))*(F2 + F3)*F4*mf*
                                                                                                                                                                                                                                                      (-((F3*(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf +  (F2*(F4 - (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*v))*
                                                                                                                                                                                                                                                            Power(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf +  (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/ Power(exp(1),(F2 + F3)*(-1 + v)),2))/Power(exp(1),(F2 + F3)*v)) - (F2*Power(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                    (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*(-1 + v)),2)*(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*v)))/
                                                                                                                                                                                                                                                         Power(exp(1),(F2 + F3)*v) + (F2*(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf +  (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/
                                                                                                                                                                                                                                                                                            Power(exp(1),(F2 + F3)*(-1 + v)))* Power(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf +  (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*v),2))/
                                                                                                                                                                                                                                                         Power(exp(1),(F2 + F3)*(-1 + v)) + (F3*Power(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),
                                                                                                                                                                                                                                                                                                      2)*(Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)))/Power(exp(1),2*(F2 + F3)*(-1 + 2*v))
                                                                                                                                                                                                                                                      ))/((F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf)*Power(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*(-1 + v)),2)
                                                                                                                                                                                                                                                          *Power(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + (F3*(-F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/Power(exp(1),(F2 + F3)*v),2)))*
    (-((b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf* (-((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*
                                                                                               (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))) + (Power(exp(1),(F2 + F3)*v)*F3*F4 + 
                                                                                                                                                                                                                                                                                                             F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))))+ y)
  LF = sum(LF)
  
  if(is.null(penalty_params)){
    penalty = 0
  } else{
    t = 1/120 - F1
    t = c(t, F1 - 7/120)
    penalty = penalty_deriv(t,penalty_params$rho,penalty_params$eps,penalty_params$k,m)
  }
  
  LF + penalty
}
residual_2 = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,y,v,penalty_params,m){
  LF = (2*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf*(F4 - (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)*((F3*(Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*
                                                                                                                                              mf))*((b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf*(1 + F2*(-1 + v)) + F3*F4*(-1 + v)))/Power(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                     F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2) + (F3*F4*(-1 + F2*(-1 + v)) + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*Power(F2,2)*mf*
                                                                                                                                                                                                                                                                                                                                   (-1 + v))/(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                              (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) - ((b0 + F1*(b1 + F1*(b2 + b3*F1)))*Power(F2,2)*mf*v + F3*F4*(-1 + F2*v))/
                                                                                            (Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) - 
                                                                                            (F3*(Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))*
                                                                                               (b2*Power(F1,2)*mf + b3*Power(F1,3)*mf + F3*F4*v + b2*Power(F1,2)*F2*mf*v + b3*Power(F1,3)*F2*mf*v + b1*F1*mf*(1 + F2*v) + b0*(mf + F2*mf*v)))/
                                                                                            Power(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2))*
          (-((b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf* (-((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*
                                                                                                     (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                              F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))) +  (Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))
                                                   /(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))
          )) + y))/(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf)
  LF = sum(LF)
  
  if(is.null(penalty_params)){
    penalty = 0
  } else{
    t = 0 - F2
    penalty = penalty_deriv(t,penalty_params$rho,penalty_params$eps,penalty_params$k,m)
  }
  
  LF + penalty
}
residual_3 = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,y,v,penalty_params,m){
  LF = (2*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf*(F4 - (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)*((F2*(F4 + F3*F4*(-1 + v) + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf*(-1 + v)))/
                                                                                            (Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) + 
                                                                                            ((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))* mf))*(-((b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf) + 
                                                                                                                                                                                                                                  Power(F3,2)*F4*(-1 + v) + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*F3*mf*(-1 + v)) )/Power(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                          F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2) - (F2*(F4 + F3*F4*v + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf*v))/
                                                                                            (Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf +  F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) - 
                                                                                            ((Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))*(-(b2*Power(F1,2)*F2*mf) - b3*Power(F1,3)*F2*mf + Power(F3,2)*F4*v + 
                                                                                                                                                                                                                   b2*Power(F1,2)*F2*F3*mf*v + b3*Power(F1,3)*F2*F3*mf*v + b0*F2*mf*(-1 + F3*v) + b1*F1*F2*mf*(-1 + F3*v)))/ Power(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                                     F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2))*(-((b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf* (-((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     (Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))) + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   (Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)))) + y))/(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf)
  LF = sum(LF)
  
  if(is.null(penalty_params)){
    penalty = 0
  } else{
    t = 0 - F3
    penalty = penalty_deriv(t,penalty_params$rho,penalty_params$eps,penalty_params$k,m)
  }
  
  LF + penalty
}
residual_4 = function(mf,b0,b1,b2,b3,F1,F2,F3,F4,y,v,penalty_params,m){
  LF = (2*Power(b0 + F1*(b1 + F1*(b2 + b3*F1)),2)*(F2 + F3)*Power(mf,2)* (-(F2/(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                  F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))) +  F2/(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                   F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)) -  (F3*(Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + 
                                                                                                                                                                                                                                                                             F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))* mf)))/Power(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                                                                              F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2) + (F3*(Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)))/
                                                                            Power(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf),2))*(-((b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf*
                                                                                                                                                                                                                                             (-((Power(exp(1),(F2 + F3)*(-1 + v))*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))/(Power(exp(1),(F2 + F3)*(-1 + v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + 
                                                                                                                                                                                                                                                                                                                                                                                    F3*((-1 + Power(exp(1),(F2 + F3)*(-1 + v)))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))) + (Power(exp(1),(F2 + F3)*v)*F3*F4 + F2*(F4 + (-1 + Power(exp(1),(F2 + F3)*v))*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf))
                                                                                                                                                                                                                                              /(Power(exp(1),(F2 + F3)*v)*(b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf + F3*((-1 + Power(exp(1),(F2 + F3)*v))*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*mf)))) + y))/(F3*F4 + (b0 + F1*(b1 + F1*(b2 + b3*F1)))*F2*mf)
  LF = sum(LF)
  
  if(is.null(penalty_params)){
    penalty = 0
  } else{
    t = 0 - F4
    t = c(t, F4 - mf*(b0 + b1 * F1 + b2 * F1^2 + b3 * F1^3))
    penalty = penalty_deriv(t,penalty_params$rho,penalty_params$eps,penalty_params$k,m)
  }
  
  LF + penalty
}