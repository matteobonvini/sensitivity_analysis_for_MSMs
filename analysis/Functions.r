library(lpSolve)
library(quantreg)
source("Fit.r")
LP_sub = function(B,Y,A,piA,piAX,gamma,index){
     ### This function bounds beta[index]
     ### under propensity ratio sensitivity for one value gamma
     ### using linear programming.
     ###
     ### Linear marginal structural model: Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = design matrix
     ### Y = outcome
     ### A = treatment
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score

     z = B[,index]
     V = B[,-index]
     f = resid(lm(z ~ 0 +V,weights=piA/piAX))
     n = length(Y)

     k = ncol(B)
     C = mean(B[,index]*f*piA/piAX)
     Q = (1/n)*(Y*f*piA/(C*piAX))
cat("Hello",sum(Q),"\n")
     H = matrix(0,k,n)
     ff = matrix(f,k,n,byrow=TRUE)
     H = t(B)*ff*piA/piAX
     H = rbind(H,rep(1/n,n))
     rhs = t(B) %*% as.matrix(f)
     rhs = c(rhs,1)

     I = diag(n)
     HH = rbind(H,I,I)
     rhs = c(rhs,rep(gamma,n),rep(1/gamma,n))
     con = c(rep("=",k+1),rep("<=",n),rep(">=",n))
     tmp = lp("max",Q,HH,con,rhs)
     up = tmp$objval
     alpha_up = tmp$solution
     tmp = lp("min",Q,HH,con,rhs)
     low = tmp$objval
     alpha_low = tmp$solution

     return(list(low=low,up=up,alpha_low=alpha_low,alpha_up=alpha_up))
     }





LP = function(B,Y,A,piA,piAX,gamma,index){
     ### This function bounds beta[index]
     ### under propensity ratio sensitivity for a vector of values of gamma
     ### using linear programming.
     ###
     ### Linear marginal structural model: Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = design matrix
     ### Y = outcome
     ### A = treatment
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score

     ng = length(gamma)
     up  = rep(0,ng)
     low = rep(0,ng)
     for(i in 1:ng){
          tmp = LP_sub(B,Y,A,piA,piAX,gamma[i],index)
          low[i] = tmp$low
          up[i]  = tmp$up
          }
      return(list(low=low,up=up))
      }
     




###Ratio_sub = function(B,Y,A,W,gamma,index){
###     ### This function bounds beta[index] using fractional linear programming
###     ### under propensity ratio sensitivity for one value gamma
###     ###
###     ### Linear marginal structural model: Y(a) = sum_j beta_j b_j(a)
###     ### B = [b1 ... bk] = design matrix
###     ### Y = outcome
###     ### A = treatment
###     ### W = piA/piAX
###     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
###     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score
###
###
###
###     n = length(Y)
###     k = ncol(B)
###
###     b1 = B[,index]
###     b2 = as.matrix(B[,-index])
###     r = resid(lm(b1 ~ 0 +b2,weights=W))
###     ry = resid(lm(Y ~ 0 +b2,weights=W))
###     f = ry*r*W
###     g = r^2*W
###
###
###
###
###     ### create constraint matrix H
###     H = matrix(0,n,k-1)
###     for(j in 1:(k-1)){
###          H[,j] = b2[,j]*r*W
###          }
###     H = t(H)
###     H = rbind(H,rep(1/n,n))
###     I = diag(n)
###     H = rbind(H,I,I)
###
###     rhs = c(rep(0,k-1),1,rep(gamma,n),rep(1/gamma,n))
###     con = c(rep("=",k),rep("<=",n),rep(">=",n))
###
###     S = rep(0,nt)
###     for(i in 1:nt){
###          h = f - tgrid[i]*g
###          tmp = lp("max",h,H,con,rhs)
###          S[i] = tmp$objval
###          v_up = tmp$solution
###          }
###     up = max(tgrid[S >= 0])
###
###     S = rep(0,nt)
###     for(i in 1:nt){
###          h = f - tgrid[i]*g
###          tmp = lp("min",h,H,con,rhs)
###          S[i] = tmp$objval
###          v_low = tmp$solution
###          }
###     low = max(tgrid[S >= 0])
###
###
###     return(list(low=low,up=up))
###     }
###

Ratio = function(B,Y,A,W,gamma,index){
     ### This function bounds beta[index] using fractional linear programming
     ### under propensity ratio sensitivity for one value gamma
     ###
     ### Linear marginal structural model: Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = design matrix
     ### Y = outcome
     ### A = treatment
     ### W = piA/piAX
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score

     aux = gamma
     ng = length(aux)
     up = rep(0,ng)
     low = rep(0,ng)
     for(j in 1:ng){
          tmp = Ratio_sub(B,Y,A,W,aux[j],index)
          low[j] = tmp$low
          up[j] = tmp$up
          }
     return(list(low=low,up=up))
     }



Homotopy = function(B,Y,W,gamma,index){
     ### This function bounds beta[index]
     ### under propensity ratio sensitivity for a vector of values gamma
     ### using homotopy: i.e. find local optimum for gamma_j
     ### starting at solution for gamma_{j-1}.
     ###
     ### Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = desifgn matrix
     ### Y = outcome
     ### A = treatment
     ### W = piA/piAX
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score
  


     ng = length(gamma)
     n = length(Y)
     up = rep(0,ng)
     low = rep(0,ng)


     ##upper bound
     V = diag(W)
     G = solve(t(B) %*%  V %*% B)
     beta = (G %*% t(B) %*% V %*% as.matrix(Y))
     up[1] = beta[index]
     d = rep(0,n)
     for(i in 1:n){
          ri = as.matrix(B[i,])
          tmp = (ri %*% t(ri)) %*% as.matrix(beta)
          si = Y[i]*ri
          d[i] = ((G %*% (si - tmp))[index])/W[i]
          }


     for(j in 2:ng){
            t = gamma[j]/(1+gamma[j])
            q = quantile(d,t)
            alpha = rep(1/gamma[j],n)
            alpha[d >= q] = gamma[j]
            V = diag(alpha*W)
            G = solve(t(B) %*%  V %*% B)
            beta = (G %*% t(B) %*% V %*% as.matrix(Y))
            up[j] = beta[index]
            for(i in 1:n){
                 ri = as.matrix(B[i,])
                 tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                 si = Y[i]*ri
                 d[i] = ((G %*% (si - tmp))[index])/W[i]
                 }
            }

 


     ##lower bound
     V = diag(W)
     G = solve(t(B) %*%  V %*% B)
     beta = (G %*% t(B) %*% V %*% as.matrix(Y))
     low[1] = beta[index]
     d = rep(0,n)
     for(i in 1:n){
          ri = as.matrix(B[i,])
          si = Y[i]*ri
          tmp = (ri %*% t(ri)) %*% as.matrix(beta)
          d[i] = ((G %*% (si - tmp))[index])/W[i]
          }
     for(j in 2:ng){
            t = 1/(1+gamma[j])
            q = quantile(d,t)
            alpha = rep(1/gamma[j],n)
            alpha[d <= q] = gamma[j]
            V = diag(alpha*W)
            G = solve(t(B) %*%  V %*% B)
            beta = (G %*% t(B) %*% V %*% as.matrix(Y))
            low[j] = beta[index]
            for(i in 1:n){
                 ri = as.matrix(B[i,])
                 si = Y[i]*ri
                 tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                 d[i] = ((G %*% (si - tmp))[index])/W[i]
                 }
            }

 
       return(list(low=low,up=up))
       }



Exact = function(B,Y,piA,piAX,gamma,index){
     B = as.matrix(B)
     n = length(Y)
     if(n > 20)stop("Can't get exact solution for n this big")
     
     W = piA/piAX
     ng  = length(gamma)
     up  = rep(0,ng)
     low = rep(0,ng)

     for(j in 1:ng){
          m = round( n/(1+gamma[j]))
          if(m < 1)m = 1
          if(m > n)m = n
          S = combn(1:n,m)
          N = ncol(S)
          b = rep(0,N)
          for(i in 1:N){
               I = S[,i]
               alpha = rep(1/gamma[j],n)
               alpha[I] = gamma[j]
               V = diag(alpha*W)
               q = solve(t(B) %*% V %*% B) %*% t(B) %*% V %*% as.matrix(Y)
               b[i] = q[index]
               }
          low[j] = min(b)
          up[j] = max(b)
          }
     return(list(low=low,up=up))
     }




Basic = function(f,gamma){
     ### bound E[f(Z)*alpha(Z)]
     ### subject to E[alpha(Z)] = 1

     f = rev(sort(f))
     n = length(f)
     mup  = round(n/(1+gamma))
     mlow = round(n*gamma/(1+gamma))
     alpha = rep(1/gamma,n)
     alpha[1:mup] = gamma
     alphaup = alpha
     up = mean(f*alpha)
     alpha = rep(gamma,n)
     alpha[1:mlow] = 1/gamma
     low = mean(f*alpha)
     alphalow = alpha
     return(list(low=low,up=up))
     }



Sherman = function(B,Y,V,A,index){
     ### compute (B^T V_i B)^{-1} B^T V_i Y 
     ### for each i for which A[i] = 1
     ### where V_{ii} is flipped to 1/V_{ii}
     n = length(Y)
     AA = solve(t(B) %*%  V %*% B)
     ident = diag(rep(1,n)) 

     old = (AA %*% t(B) %*% V %*% as.matrix(Y))[index]
     b = rep(old,n)

     for(i in 1:n){
       if(A[i] == 1){
            Delta = (1/V[i,i]) - V[i,i]
            ri = as.matrix(B[i,])
            ei = as.matrix(ident[i,])
            one = Delta * (AA %*% ri %*% t(ri) %*% AA)
            two = c(1 + Delta* (t(ri) %*% AA %*% ri))
            three = V + Delta*(ei %*% t(ei))
            tmp = (AA - one/two) %*% t(B) %*% three %*% as.matrix(Y)
            b[i] = tmp[index]
            }
          }

     return(b)
     }


Coordinate_Wise_Sub = function(B,Y,piA,piAX,gamma,index,m=3,which){
     ### Cooordinate ascent/descent
     ### m = number of times to cycle through coordinates
     ### which = "max" or "min"
     ### One value of gamma
     B = as.matrix(B)
     n = length(Y)
     W = piA/piAX
     V = diag(W)
     d = ncol(B)


     beta = solve(t(B) %*% V %*% B) %*% t(B) %*% V %*% as.matrix(Y)
     b = beta[index]

     ### random shuffle
     I = sample(1:n)
     B = as.matrix(B[I,])
     Y = Y[I]
     W = W[I]  
     V = diag(W)
     alpha = rep(1,n)
     for(j in 1:m){
          for(i in 1:n){
               V1 = V
               V2 = V
               V1[i,i] = gamma
               V2[i,i] = 1/gamma
     
               b1 = solve(t(B) %*% V1 %*% B) %*% t(B) %*% V1 %*% as.matrix(Y)      
               b2 = solve(t(B) %*% V2 %*% B) %*% t(B) %*% V2 %*% as.matrix(Y)      
               if(which == "max"){
                    if(b1[index] >= b2[index])V = V1 else V = V2
                    }
               if(which == "min"){
                    if(b1[index] <= b2[index])V = V1 else V = V2
                    }
               }
          }
     beta = solve(t(B) %*% V %*% B) %*% t(B) %*% V %*% as.matrix(Y)
     return(beta[index])
     }






Coordinate_Wise = function(B,Y,piA,piAX,gamma,index,N,m=3){
     ### Coordinate ascent/descent
     ### Vector of gamma values
     ### m = number of cycles through coordinates
     ### N = number of random starting positions
     ng  = length(gamma)
     Up  = rep(0,ng)
     Low = rep(0,ng)


     for(j in 1:ng){
          up = rep(0,N)
          low = rep(0,N)
          for(i in 1:N){
               up[i]  = Coordinate_Wise_Sub(B,Y,piA,piAX,gamma[j],index,m,which="max")
               low[i] = Coordinate_Wise_Sub(B,Y,piA,piAX,gamma[j],index,m,which="min")
               }
          Up[j]  = max(up)
          Low[j] = min(low)
          }
     return(list(low=Low,up=Up))
     }








Outcome_Exact = function(B,A,X,Y,delta,index,piA,piAX,piAXU){
     ### Outcome sensitivity
     ### exact bound
     ### requires piAXU = pi(A|U,X)
     n = length(Y)
     if(n > 15)stop("Cannot compute exact bound with n > 15")
     nd = length(delta)
     W = diag(piA/piAX)
     WW = diag(piA/piAXU)
     E  = expand.grid(replicate(n, 0:1, simplify = FALSE))
     E[E == 0 ]=-1

     tmp = lm(Y ~ A + X)
     mu = tmp$fitted


     G = solve(t(B) %*% W %*% B)
     beta = G %*% t(B) %*% W %*% as.matrix(mu)

     offset = G %*% t(B) %*% WW %*% t(E)

     offset = offset[index,]
     up = beta[index] + delta*max(offset)
     low = beta[index] + delta*min(offset)


     return(list(low=low,up=up))
     }



Outcome_Bound = function(B,A,X,Y,delta,index,piA,piAX){
     ### outcome sensitivity model
     ### |mu(u,x,a) - mu(x,a)| <= delta
     ### Y(a) = sum_j b_j(a)
     ### B = [b1 ... bk]
     ### Bound beta[index]
 
     n = length(Y)
     nd = length(delta)
     W = diag(piA/piAX)

##     tmp = lm(Y ~ A + X)
tmp = lm(Y~ 0+ B + X)
     mu = tmp$fitted

     G = solve(t(B) %*% W %*% B)
     beta = G %*% t(B) %*% W %*% as.matrix(mu)

     r = (G %*% t(B))[index,]


     bbar = apply(abs(B),2,mean)   



     low = beta[index] - delta*sum(abs(r))
     up  = beta[index] + delta*sum(abs(r))

     return(list(low=low,up=up))
     }




Subset_sub = function(B,Y,A,piA,piAX,epsilon,index,a,b){
     ### (1-epsilon) P_0 + \epsilon P_1
     ### where P_0 is unconfounded
     ### bound beta_1(P_0)
     ### This is the subroutine for one value of epsilon

     nt = 100
     tgrid = seq(a,b,length=nt)
     n = length(Y)
     k = ncol(B)

     b1 = B[,index]
     b2 = as.matrix(B[,-index])
     r = resid(lm(b1 ~ 0 +b2,weights=piA/piAX))
     ry = resid(lm(Y ~ 0 +b2,weights=piA/piAX))
     f = ry*r*piA/piAX
     g = r^2*piA/piAX

     S = rep(0,nt)
     for(i in 1:nt){
          h = f - tgrid[i]*g
          q = quantile(h,epsilon)
          S[i] = mean( h*(h>q))
          }
     i = which.min(abs(S))
     up = tgrid[i]

     S = rep(0,nt)
     for(i in 1:nt){
          h = f - tgrid[i]*g
          q = quantile(h,1-epsilon)
          S[i] = mean( h*(h<q))
          }
     i = which.min(abs(S))
     low = tgrid[i]

     return(list(low=low,up=up))
     }

     

Subset = function(B,Y,A,piA,piAX,epsilon,index,a,b){
     ### (1-epsilon) P_0 + \epsilon P_1
     ### where P_0 is unconfounded
     ### bound beta_1(P_0)
     neps = length(epsilon)
     up  = rep(0,neps)
     low = rep(0,neps)
     for(i in 1:neps){
          tmp = Subset_sub(B,Y,A,piA,piAX,epsilon[i],index,a,b)
          low[i] = tmp$low
          up[i]  = tmp$up
          }
      return(list(low=low,up=up))
      }
     




Subset_sub_propensity = function(B,Y,A,piA,piAX,epsilon,gamma,index){
     ### (1-epsilon) P_0 + \epsilon P_1
     ### P1 in gamma propensity model
     ### where P_0 is unconfounded
     ### bound beta_1(P)
     ### This is the subroutine for one value of epsilon
     n = length(Y)
     k = ncol(B)

     z = B[,index]
     V = B[,-index]
     f = resid(lm(z ~ 0 +V))
     C = mean(B[,index]*f)

     xi = Y*f*piA/(C*piAX)
     xi = sort(xi)
     m1 = round(n*epsilon)
     if(m1 < 1)m1 = 1
     if(m1 > n)m1 = n
     m2 = round(n*(1-epsilon))
     if(m2 < 1)m2 = 1
     if(m2 > n)m2 = n

     low = sum(xi[1:m2])/(n*(1-epsilon))
     up  = sum(xi[m1:n])/(n*(1-epsilon))

     return(list(low=low,up=up))
     }



##############

Homotopy_subset_propensity = function(B,Y,piA,piAX,gamma,epsilon,index){
     ### This function bounds beta[index]
     ### under propensity ratio sensitivity for a vector of values gamma
     ### using homotopy: i.e. find local optimum for gamma_j
     ### starting at solution for gamma_{j-1}.
     ###
     ### Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = desifgn matrix
     ### Y = outcome
     ### A = treatment
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score
     ### Mode:
     ### (1-epsilon) P_0 + \epsilon P_1
     ### P1 in gamma propensity model
     ### where P_0 is unconfounded
     ### bound beta_1(P)
  

     W = piA/piAX
     ng = length(gamma)
     neps = length(epsilon)
     n = length(Y)
     up = matrix(0,ng,neps)
     low = matrix(0,ng,neps)

     for(s in 1:neps){
          ##upper bound
          V = diag(W)
          G = solve(t(B) %*%  V %*% B)
          beta = (G %*% t(B) %*% V %*% as.matrix(Y))
          up[1,s] = beta[index]
          d = rep(0,n)
          for(i in 1:n){
               ri = as.matrix(B[i,])
               tmp = (ri %*% t(ri)) %*% as.matrix(beta)
               si = Y[i]*ri
               d[i] = ((G %*% (si - tmp))[index])/W[i]
               }
          for(j in 2:ng){
                 t = epsilon[s]*gamma[j]/(1+gamma[j])
                 q1 = quantile(d,t)
                 t = 1-epsilon[s]/(1+gamma[j])
                 q2 = quantile(d,t)
                 alpha = rep(1,n)
                 alpha[d > q2] = gamma[j]
                 alpha[d < q1] = 1/gamma[j]
                 V = diag(alpha*W)
                 G = solve(t(B) %*%  V %*% B)
                 beta = (G %*% t(B) %*% V %*% as.matrix(Y))
                 up[j,s] = beta[index]
                 for(i in 1:n){
                      ri = as.matrix(B[i,])
                      tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                      si = Y[i]*ri
                      d[i] = ((G %*% (si - tmp))[index])/W[i]
                      }
                 }
     
      
     
     
          ##lower bound
          V = diag(W)
          G = solve(t(B) %*%  V %*% B)
          beta = (G %*% t(B) %*% V %*% as.matrix(Y))
          low[1,s] = beta[index]
          d = rep(0,n)
          for(i in 1:n){
               ri = as.matrix(B[i,])
               si = Y[i]*ri
               tmp = (ri %*% t(ri)) %*% as.matrix(beta)
               d[i] = ((G %*% (si - tmp))[index])/W[i]
               }
          for(j in 2:ng){
                 t = epsilon[s]/(1+gamma[j])
                 q1 = quantile(d,t)
                 t = 1-epsilon[s]*gamma[j]/(1+gamma[j])
                 q2 = quantile(d,t)
                 alpha = rep(1,n)
                 alpha[d > q2] = 1/gamma[j]
                 alpha[d < q1] = gamma[j]
                 V = diag(alpha*W)
                 G = solve(t(B) %*%  V %*% B)
                 beta = (G %*% t(B) %*% V %*% as.matrix(Y))
                 low[j,s] = beta[index]
                 for(i in 1:n){
                      ri = as.matrix(B[i,])
                      si = Y[i]*ri
                      tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                      d[i] = ((G %*% (si - tmp))[index])/W[i]
                      }
                 }
     
           }
       return(list(low=low,up=up))
       }





Subset_Outcome = function(B,A,X,Y,delta,epsilon,index,piA,piAX){
     ### outcome sensitivity model
     ### (1-epsilon)P_0 + epsilon P_1
     ### |mu_1(u,x,a) - mu_1(x,a)| <= delta
     ### Y(a) = sum_j b_j(a)
     ### B = [b1 ... bk]
     ### Bound beta[index]
 
     nd = length(delta)
     neps = length(epsilon)
     low = matrix(0,nd,neps)
     up  = matrix(0,nd,neps)
     n = length(Y)

     W = diag(piA/piAX)

     tmp = lm(Y ~ A + X)
     mu = tmp$fitted

     G = solve(t(B) %*% W %*% B)
     beta = G %*% t(B) %*% W %*% as.matrix(mu)

     r = (G %*% t(B))[index,]


     bbar = apply(abs(B),2,mean)   

     for(i in 1:neps){

          low[,i] = (1-epsilon[i])*beta[index] - epsilon[i]*delta*sum(abs(r))
          up[,i]  = (1-epsilon[i])*beta[index] + epsilon[i]*delta*sum(abs(r))
          }

     return(list(low=low,up=up))
     }






Ratio_sub = function(B,Y,A,W,gamma,index){
     ### This function bounds beta[index] using fractional linear programming
     ### under propensity ratio sensitivity for one value gamma
     ###
     ### Linear marginal structural model: Y(a) = sum_j beta_j b_j(a)
     ### B = [b1 ... bk] = design matrix
     ### Y = outcome
     ### A = treatment
     ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
     ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score



     n = length(Y)
     k = ncol(B)

     b1 = B[,index]
     b2 = as.matrix(B[,-index])
     r = resid(lm(b1 ~ 0 +b2,weights=W))
     ry = resid(lm(Y ~ 0 +b2,weights=W))
     f = ry*r*W
     g = r^2*W


     ### create constraint matrix H
     H = matrix(0,n,k-1)
     for(j in 1:(k-1)){
          H[,j] = b2[,j]*r*W
          }
     H = t(H)
     H = rbind(H,rep(1/n,n))
     H = cbind(H,c(rep(0,k-1),-1))
     I = diag(n)
     H = rbind(H,cbind(I,rep(-gamma,n)))
     H = rbind(H,cbind(-I,rep(1/gamma,n)))
     H = rbind(H,c(rep(0,n),1))
     H = rbind(H,c(g,0))

     rhs = c(rep(0,k),rep(0,n),rep(0,n),0,1)
     con = c(rep("=",k),rep("<=",n),rep("<=",n),">=","=")

     h = c(f,0)
     tmp = lp(direction="max",objective.in=h,const.mat=H,const.dir=con,const.rhs=rhs)
     up = tmp$objval
     tmp = lp(direction="min",objective.in=h,const.mat=H,const.dir=con,const.rhs=rhs)
     low = tmp$objval

     return(list(low=low,up=up))
     }




Fit_Bounds = function(A,Y,State,gamma,degree,delta=4){
     ## Fit the marginal structural model
     ## Y = weekly deaths
     ## log Y =  beta M_t + sum_{j=1}^k beta_j basis_j(t)
     ## k = degree
     ## NOTE: degree is a vector, one for each state
     ## get bounds for each state

     n = nrow(Y)
     m = ncol(Y)
     ng = length(gamma)

     YY = Y
     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:m
     Ypast = YY[,1:(ncol(YY)-delta-1)]
     Ypast = cbind(rep(0,n),Ypast)
     beta = rep(0,n)
     se = rep(0,n)
     one = rep(1,m)

     up  = matrix(0,n,ng)
     low = matrix(0,n,ng)
     up_local = matrix(0,n,ng)
     low_local = matrix(0,n,ng)


     for(i in 1:n){
          a = A[i,]
          a = a-a[1]
          M = cumsum(a)
          y = Y[i,]
          L = log(y+1)
          ypast = Ypast[i,]
          w = Get_Weights_One_State(a,ypast)
          basis = poly(time,degree[i])
          tmp = lm(L ~ M + basis,weights=w)
          B = cbind(M,basis)
          one = rep(1,nrow(B))
          B = cbind(one,B)        

          tmp = Homotopy(B,L,w,gamma,2)
          up[i,] = tmp$up
          low[i,] = tmp$low

          tmp = Local(B,L,w,gamma,2)
          up_local[i,] = tmp$up
          low_local[i,] = tmp$low

          }



     return(list(low=low,up=up,low_local=low_local,up_local=up_local))
     }


          


Local = function(B,Y,W,gamma,index){
     ### bound for gamma not too far from 1
 
     n = length(Y)
     d = ncol(B)
     M = matrix(0,d,d)
     for(i in 1:n){
         b = as.matrix(B[i,])
         M = M + (W[i]/n)*(b %*% t(b))
         }
 
     alpha = t(B) %*%  as.matrix(Y*W)/n
     beta = solve(M) %*% alpha
 
     deriv = rep(0,n)
     Minv = solve(M)
     for(i in 1:n){
          b = as.matrix(B[i,])
          r = Y[i] - sum(b*beta)
          tmp = W[i]*r*b
          deriv[i] = (Minv %*% tmp)[index]
          }
     deriv = mean(abs(deriv))
 
     low = beta[index] - log(gamma)*deriv
     up  = beta[index] + log(gamma)*deriv
 
     return(list(low=low,up=up))
     }






Homotopy_g = function(B,Y,W,gamma){
     ### This function bounds g(a_t,beta) = B^T b(a)
     ### where a = Observed A shifted back 2 weeks


  
     small = .00001
     W[W == 0] = small

     ng = length(gamma)
     n = length(Y)
     up = rep(0,ng)
     low = rep(0,ng)

     Up = matrix(0,n,ng)
     Low = matrix(0,n,ng)

     a = B[,2]
     q = length(a)
     a = c(a[3:q],a[q],a[q])


for(time in 1:n){
     ##upper bound
     V = diag(W)
     G = solve(t(B) %*%  V %*% B)
     beta = (G %*% t(B) %*% V %*% as.matrix(Y))


     c = matrix(B[time,],ncol=1)
     c[2] = a[time] ##shift A by minus 2 weeks


     Up[time,1] = sum(c*beta)
     k = length(beta)
     D = matrix(0,n,k)
     d = rep(0,n)

     for(i in 1:n){
          ri = as.matrix(B[i,])
          tmp = (ri %*% t(ri)) %*% as.matrix(beta)
          si = Y[i]*ri
##          d[i,] = ((G %*% (si - tmp))[index])/W[i]
          D[i,] = ((G %*% (si - tmp)))/W[i]
          d[i] = sum(D[i,]*c)
          }




     for(j in 2:ng){
            t = gamma[j]/(1+gamma[j])
            q = quantile(d,t)
            alpha = rep(1/gamma[j],n)
            alpha[d >= q] = gamma[j]
            V = diag(alpha*W)
            G = solve(t(B) %*%  V %*% B)
            beta = (G %*% t(B) %*% V %*% as.matrix(Y))
            Up[time,j] = sum(c*beta)
            for(i in 1:n){
                 ri = as.matrix(B[i,])
                 tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                 si = Y[i]*ri
                 D[i,] = ((G %*% (si - tmp)))/W[i]
                 d[i] = sum(D[i,]*c)
                 }
            }

 


     ##lower bound
     V = diag(W)
     G = solve(t(B) %*%  V %*% B)
     beta = (G %*% t(B) %*% V %*% as.matrix(Y))
     c = matrix(B[time,],ncol=1)
     c[2] = a[time] ##shift A by minus 2 weeks
     Low[time,1] = sum(c*beta)
     d = rep(0,n)
     D = matrix(0,n,k)



     for(i in 1:n){
          ri = as.matrix(B[i,])
          si = Y[i]*ri
          tmp = (ri %*% t(ri)) %*% as.matrix(beta)
          D[i,] = ((G %*% (si - tmp)))/W[i]
          d[i] = sum(D[i,]*c)
          }

     for(j in 2:ng){
            t = 1/(1+gamma[j])
            q = quantile(d,t)
            alpha = rep(1/gamma[j],n)
            alpha[d <= q] = gamma[j]
            V = diag(alpha*W)
            G = solve(t(B) %*%  V %*% B)
            beta = (G %*% t(B) %*% V %*% as.matrix(Y))
            Low[time,j] = sum(c*beta)
            for(i in 1:n){
                 ri = as.matrix(B[i,])
                 si = Y[i]*ri
                 tmp = (ri %*% t(ri)) %*% as.matrix(beta)
                 D[i,] = ((G %*% (si - tmp)))/W[i]
                 d[i] = sum(D[i,]*c)
                 }
            }


       }
 
       return(list(Low=Low,Up=Up))
       }





Fit_Bounds_g = function(A,Y,State,gamma,degree,delta=4){
     ## Fit the marginal structural model
     ## Y = weekly deaths
     ## log Y =  beta M_t + sum_{j=1}^k beta_j basis_j(t)
     ## k = degree
     ## NOTE: degree is a vector, one for each state
     ## get bounds on g(a_t) for each state

     n = nrow(Y)
     m = ncol(Y)
     ng = length(gamma)

     YY = Y
     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:m
     Ypast = YY[,1:(ncol(YY)-delta-1)]
     Ypast = cbind(rep(0,n),Ypast)
     beta = rep(0,n)
     se = rep(0,n)
     one = rep(1,m)

     Up  = list()
     Low = list()
     up_local = matrix(0,n,ng)
     low_local = matrix(0,n,ng)

     Deaths = Y

     for(i in 1:n){
          a = A[i,]
          a = a-a[1]
          M = cumsum(a)
          y = Y[i,]
          L = log(y+1)
          ypast = Ypast[i,]
          w = Get_Weights_One_State(a,ypast)
          basis = poly(time,degree[i])
          tmp = lm(L ~ M + basis,weights=w)
          B = cbind(M,basis)
          one = rep(1,nrow(B))
          B = cbind(one,B)        

          tmp = Homotopy_g(B,L,w,gamma)

          Up[[i]]  = tmp$Low
          Low[[i]] = tmp$Up

#          tmp = Local(B,L,w,gamma,2)
#          up_local[i,] = tmp$up
#          low_local[i,] = tmp$low

          }



     return(list(Low=Low,Up=Up,low_local=low_local,up_local=up_local,Deaths=Deaths))
     }


