#### Functions

library(lubridate)
library(dplyr)
library(tidyr)
library(mgcv)
library(MASS)
#library(locfit)
library(cgam)
library(mvtnorm)
library(splines)
library(sandwich)
library(car)


we = function(a){
     ### convert vector a to weekly averages
     n = length(a)
     b = floor(n/7)
     nn = b*7
     a = a[1:nn]
     a = matrix(a,b,7,byrow=TRUE)
     a = apply(a,1,sum,na.rm=TRUE)
     return(a)
     }







Fit = function(A,Y,State,weights=TRUE,nboot=1000,gamma=.1,degree,delta=4){
     ## Fit the marginal structural model
     ## Y = weekly deaths
     ## log Y =  beta M_t + sum_{j=1}^k beta_j basis_j(t)
     ## k = degree
     ## NOTE: degree is a vector, one for each state


     n = nrow(Y)
     m = ncol(Y)

     YY = Y
     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:m

     Ypast = YY[,1:(ncol(YY)-delta-1)]
     Ypast = cbind(rep(0,n),Ypast)



     beta = rep(0,n)
     se = rep(0,n)
     

     E = matrix(0,n,m)
     Low = matrix(0,n,m)
     Up = matrix(0,n,m)
     Weights = matrix(0,n,m)

     Counterfactual1    = matrix(0,n,m)   ### One week earlier
     CounterfactualUp1  = matrix(0,n,m)
     CounterfactualLow1 = matrix(0,n,m)
     Counterfactual2    = matrix(0,n,m)   ### Two weeks earlier
     CounterfactualUp2  = matrix(0,n,m)
     CounterfactualLow2 = matrix(0,n,m)
     Counterfactual3    = matrix(0,n,m)   ### Stay Vigilant
     CounterfactualUp3  = matrix(0,n,m)
     CounterfactualLow3 = matrix(0,n,m)

     Smooth = matrix(0,n,m)

     one = rep(1,m)


     Total = rep(0,n)
     TotalUp = rep(0,n)
     TotalLow = rep(0,n)

     Ysum = rep(0,n)


     Lout = NULL
     FIT = NULL
     RESIDUALS = NULL


     for(i in 1:n){
          a = A[i,]
          a = a-a[1]
          M = cumsum(a)
          y = Y[i,]
          L = log(y+1)
          ypast = Ypast[i,]


          w = Get_Weights_One_State(a,ypast)
          Weights[i,] = w



          basis = poly(time,degree[i])
          tmp = lm(L ~ M + basis,weights=w)
          Lout = rbind(Lout,L)
          FIT = rbind(FIT,tmp$fitted)
          RESIDUALS = rbind(RESIDUALS,resid(tmp))
          beta[i] = tmp$coef[2]
          acov =  NeweyWest(tmp,verbose=TRUE,prewhite=FALSE,adjust=TRUE)  ###this adjusts for serial correlation

          se[i] = sqrt(acov[2,2])
          

          ### smooth time component
          b = tmp$coef[3:(ncol(basis)+2)]
          b = matrix(b,ncol=1)
          Smooth[i,] = basis %*% b



          ### shift one intervention
          astar = a
          astar = c(astar[-1],astar[m])
          Mstar = cumsum(astar)
          X = cbind(rep(1,m),Mstar,basis)
          B = as.matrix(tmp$coef)
          Counterfactual1[i,] = exp(X %*% B)
          S = rep(0,m)
          for(j in 1:m){
               x = as.matrix(X[j,])
               S[j] = sqrt(t(x) %*% acov %*% x)
               }
          CounterfactualUp1[i,]  = Counterfactual1[i,]*exp(2*S)
          CounterfactualLow1[i,] = Counterfactual1[i,]*exp(-2*S)


          ### shift two intervention
          astar = a
          astar = c(astar[-(1:2)],astar[c(m-1,m)])
          Mstar = cumsum(astar)
          X = cbind(rep(1,m),Mstar,basis)
          B = as.matrix(tmp$coef)
          Counterfactual2[i,] = exp(X %*% B)
          S = rep(0,m)
          for(j in 1:m){
               x = as.matrix(X[j,])
               S[j] = sqrt(t(x) %*% acov %*% x)
               }
          CounterfactualUp2[i,]  = Counterfactual2[i,]*exp(2*S)
          CounterfactualLow2[i,] = Counterfactual2[i,]*exp(-2*S)


          ###stay vigilant
          r = c(rep(1,5),rep(2,m-5))
          astar = rep(a,r)
          astar = astar[1:m]
          Mstar = cumsum(astar)
          X = cbind(rep(1,m),Mstar,basis)
          B = as.matrix(tmp$coef)
          Counterfactual3[i,] = exp(X %*% B)
          S = rep(0,m)
          for(j in 1:m){
               x = as.matrix(X[j,])
               S[j] = sqrt(t(x) %*% acov %*% x)
               }
          CounterfactualUp3[i,]  = Counterfactual3[i,]*exp(2*S)
          CounterfactualLow3[i,] = Counterfactual3[i,]*exp(-2*S)



### Total saved under stay vigilant
            Bstar = rmvnorm(nboot,B,acov)
            thetastar = rep(0,nboot)
            for(j in 1:nboot){
                 o = as.matrix(Bstar[j,])
                 thetastar[j] = sum(exp(X %*% o)) - sum(y)
                 }
             Total[i] = sum(Counterfactual3[i,]) - sum(y)
             TotalUp[i] = quantile(thetastar,.975)
             TotalLow[i] = quantile(thetastar,.025)
             Ysum[i] = sum(y)

          }



     return(list(beta=beta,se=se,Up=Up,Low=Low,E=E,time=time,
            Counterfactual1=Counterfactual1,CounterfactualUp1=CounterfactualUp1,CounterfactualLow1=CounterfactualLow1,
            Counterfactual2=Counterfactual2,CounterfactualUp2=CounterfactualUp2,CounterfactualLow2=CounterfactualLow2,
            Counterfactual3=Counterfactual3,CounterfactualUp3=CounterfactualUp3,CounterfactualLow3=CounterfactualLow3,
            Smooth=Smooth,Total=Total,TotalUp=TotalUp,TotalLow=TotalLow,Y=Y,Ysum=Ysum,
            Lout=Lout,FIT=FIT,RESIDUALS=RESIDUALS))
     }


          




Get_Weights = function(A,Y,linear=TRUE,h=5){
     ### markov moment matching method

     ### uses, A, A^2, log(Y) and (log(Y))^2; 


     n = nrow(Y)
     m = ncol(Y)
     Weights = matrix(0,n,m)
     one = rep(1,m-2) 
     one = matrix(one,ncol=1)
     D = rep(0,5)
     D[1] = m-2
     D = matrix(D,ncol=1)
     for(i in 1:n){
          
          a = A[i,]
          L = log(Y[i,]+1)
          an = a[3:m]
          an2 = an^2
          ap = a[2:(m-1)]
          Ln = L[2:(m-1)]
          Ln2 = Ln^2
          Lp = L[1:(m-2)]
          
          if(linear == TRUE){
               mu  = fitted(rlm(an ~ ap))
               mu2 = fitted(rlm(an2 ~ ap))
               nu  = fitted(rlm(Ln ~ Lp+ap))
               nu2 = fitted(rlm(Ln2 ~ Lp+ap))
               }

          if(linear == FALSE){  ###use nonparametric fit
               mu = fitted(locfit(an ~ lp(ap,deg=1,h=h)))
               mu2 = fitted(locfit(an2 ~ lp(ap,deg=1,h=h)))
               nu = fitted(locfit(Ln ~ lp(Lp,ap,deg=1,h=h)))
               nu2 = fitted(locfit(Ln2 ~ lp(Lp,ap,deg=1,h=h)))
               }


          H = matrix(1,m-2,5)
          H[,2] = (an - mu)*(Ln-nu)
          H[,3] = (an2 - mu2)*(Ln-nu)
          H[,4] = (an - mu)*(Ln2-nu2)
          H[,5] = (an2 - mu2)*(Ln2-nu2)

          W = one - H %*% solve(t(H) %*% H) %*% (t(H) %*% one - D)
          W[W  < 0] = 0 
          W = c(1,1,W)
          W = m*W/sum(W)
    
          Weights[i,] = W
          }


     return(Weights)
     }




Get_Weights_One_State = function(a,y,linear=TRUE,h=10){
     ### markov moment matching method for one state

     ### uses, A, A^2, log(Y) and (log(Y))^2; 


     m = length(y)
     one = rep(1,m-2) 
     one = matrix(one,ncol=1)
     D = rep(0,5)
     D[1] = m-2
     D = matrix(D,ncol=1)
     L = log(y+1)
     an = a[3:m]
     an2 = an^2
     ap = a[2:(m-1)]
     Ln = L[2:(m-1)]
     Ln2 = Ln^2
     Lp = L[1:(m-2)]
     

     if(linear == TRUE){
          mu  = fitted(lm(an ~ ap))
          mu2 = fitted(lm(an2 ~ ap))
          nu  = fitted(lm(Ln ~ Lp+ap))
          nu2 = fitted(lm(Ln2 ~ Lp+ap))
          }



     if(linear == FALSE){  ###use nonparametric fit
          mu = fitted(locfit(an ~ lp(ap,deg=1,h=h)))
          mu2 = fitted(locfit(an2 ~ lp(ap,deg=1,h=h)))
          nu = fitted(locfit(Ln ~ lp(Lp,ap,deg=1,h=h)))
          nu2 = fitted(locfit(Ln2 ~ lp(Lp,ap,deg=1,h=h)))
          }


     H = matrix(1,m-2,5)
     H[,2] = (an - mu)*(Ln-nu)
     H[,3] = (an2 - mu2)*(Ln-nu)
     H[,4] = (an - mu)*(Ln2-nu2)
     H[,5] = (an2 - mu2)*(Ln2-nu2)

     W = one - H %*% solve(t(H) %*% H) %*% (t(H) %*% one - D)
     W[W  < 0] = 0 
     W = c(1,1,W)
     W = m*W/sum(W)

     Weights = W



     return(Weights)
     }


Blip = function(A,Y,State,degree,delta=1){
     ### Parametric  blip outcome regression
     ### L_t = L_{t-1} + beta A_t + f(t) + epsilon_t 
     ### To be consistent with MSM with degree k, we use othogonal models of order degree k-1

     degree = degree - 1


     n = nrow(Y)
     m = ncol(Y)
     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:(m-1)
     beta = rep(0,n)
     se = rep(0,n)

     aic = rep(0,n)

     if(degree > 0)basis = poly(time,degree)

     Fit = matrix(0,n,m-1)

     for(i in 1:n){
          a = A[i,]
          a = a - a[1]
          M = cumsum(a)
          y = Y[i,]
          y[y<0]=0
          L = log(y+1)


          an = a[-1]
          ap = a[1:(m-1)]
          Ln = L[-1]
          Lp = L[1:(m-1)]


          D = Ln - Lp
          if (degree > 0)tmp = lm(D ~ an + basis)
          if (degree == 0)tmp = lm(D ~ an)
          beta[i] = tmp$coef[2]
          se[i] = summary(tmp)$coef[2,2]
          aic[i] = AIC(tmp)

          Fit[i,] = Lp + fitted(tmp)
          }

     Y = Y[,-1]
     return(list(beta=beta,se=se,time=time,Fit=Fit,Y=Y,aic=aic))
     }









Sensitivity = function(A,Y,State,degree,delta=4){
     ## Y = weekly deaths
     ## log Y = s(t) + beta M_t 
     ## Find delta s.t. unoberved confouding of size Delta*se makes effect non-significant
     n = nrow(Y)
     m = ncol(Y)
     Ypast = Y[,1:(ncol(Y)-delta-1)]
     Ypast = cbind(rep(0,n),Ypast)

     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:m

     UP =  matrix(0,n,m)     
     LOW = matrix(0,n,m)     
     WEIGHTS = matrix(0,n,m)

     beta = rep(0,n)
     se.beta = rep(0,n)
     Delta = rep(0,n)



     for(i in 1:n){
          a = A[i,]
          a = a-a[1]
          M = cumsum(a)
          y = Y[i,]
          L = log(y+1)
          ypast = Ypast[i,]

          w = Get_Weights_One_State(a,ypast)
          WEIGHTS[i,] = w
          df = data.frame(L=L,M=M,time=time)

          basis = poly(time,degree[i])
          tmp = lm(L ~ M + basis,weights=w)
          beta[i] = tmp$coef[2]
          se.beta[i] = summary(tmp)$coef[2,2]

          if(beta[i] > 2*se.beta[i])Delta[i] = beta[i]/se.beta[i] - 2
          if(beta[i] < -2*se.beta[i])Delta[i] = -beta[i]/se.beta[i] - 2


          }


          
     return(list(beta=beta,se.beta=se.beta,Delta=Delta))
     }

          



Sensitivity_Gamma = function(A,Y,h=4,State,G,grid,degree,niter=100,delta=4){
     ## Y = weekly deaths
     ## log T = s(t) + beta M_t where s(t) = sum_j basis_j(t)
     ## Sensitivity analysis
     n = nrow(Y)
     m = ncol(Y)
     AA = A
     YY = Y
     Y = Y[,(delta+1):m]
     A = A[,1:(m-delta)] 
     m = ncol(Y)
     time = 1:m

     Ypast = YY[,1:(ncol(YY)-delta-1)]
     Ypast = cbind(rep(0,n),Ypast)


     WEIGHTS = matrix(0,n,m)

     beta = rep(0,n)
     se.beta = rep(0,n)
     change = rep(0,n)
     ww = rep(1,m)



     up = rep(0,n)
     low = rep(0,n)

     for(i in 1:n){
          basis = poly(time,degree[i])
          a = A[i,]
          a = a-a[1]
          M = cumsum(a)
          y = Y[i,]
          L = log(y+1)


          ypast = Ypast[i,]

          w = Get_Weights_One_State(a,ypast)
          WEIGHTS[i,] = w




          tmp = lm(L ~ M + basis,weights=w)
          beta[i] = tmp$coef[2]
          se.beta[i] = summary(tmp)$coef[2,2]

          X = cbind(rep(1,m),M,basis)
          index = 2 
          tmp = Gamma(X,L,G,w,index,niter)
          low[i] = tmp[1]
          up[i] = tmp[2]
          }


          
     return(list(beta=beta,se.beta=se.beta,low=low,up=up))
     }

          



Gamma = function(X,y,G,w,index,niter=100){
     ### max and min of  beta.hat[index]
     ### this is used in Sensitivity_Gamma
     n = length(y)
     V = diag(w)
     Y = matrix(y,ncol=1)
     Vmax = V
     Vmin = V

     for(j in 1:niter){
          for(i in 1:n){
               v1 = Vmax
               v1[i,i] = G*w[i]
               v2 = Vmax
               v2[i,i] = w[i]/G
               b1 = solve(t(X) %*% v1 %*% X) %*% t(X) %*% v1 %*% Y
               b1 = b1[index]
               b2 = solve(t(X) %*% v2 %*% X) %*% t(X) %*% v2 %*% Y
               b2 = b2[index]
               if(b1 > b2)Vmax[i,i] = G*w[i] 
               if(b1 < b2)Vmax[i,i] = w[i]/G
               if(b1 == b2)Vmax[i,i] = w[i]


               v1 = Vmin
               v1[i,i] = G*w[i]
               v2 = Vmin
               v2[i,i] = w[i]/G
               b1 = solve(t(X) %*% v1 %*% X) %*% t(X) %*% v1 %*% Y
               b1 = b1[index]
               b2 = solve(t(X) %*% v2 %*% X) %*% t(X) %*% v2 %*% Y
               b2 = b2[index]
               if(b1 < b2)Vmin[i,i] = G*w[i] 
               if(b1 > b2)Vmin[i,i] = w[i]/G
               if(b1 == b2)Vmin[i,i] = w[i]


               }
          }


          up =  (solve(t(X) %*% Vmax %*% X) %*% t(X) %*% Vmax %*% Y)[index]
          low = (solve(t(X) %*% Vmin %*% X) %*% t(X) %*% Vmin %*% Y)[index]


     return(c(low,up))
     }




