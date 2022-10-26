setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Local = function(B,Y,piA,piAX,gamma,index){
     ### local bound for gamma not too far from 1
     ### using F(v) and E[V]=1
 
     n = length(Y)
     d = ncol(B)
     W = piA/piAX
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






Local_refined = function(B,Y,piA,piAX,gamma,index){
     ### bound for gamma not too far from 1
 
     n = length(Y)
     d = ncol(B)
     ng = length(gamma)
     W = piA/piAX
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

     low = rep(0,ng)
     up  = rep(0,ng)
     for(i in 1:ng){
          q = quantile(deriv,gamma[i]/(1+gamma[i]))
          L = rep(-log(gamma[i]),n)
          L[deriv > q] = log(gamma[i])
          up[i] = beta[index] + mean(L*deriv)

          q = quantile(deriv,1/(1+gamma[i]))
          L = rep(-log(gamma[i]),n)
          L[deriv < q] = log(gamma[i])
          low[i] = beta[index] + mean(L*deriv)
          }


 
     return(list(low=low,up=up))
     }





LocalG = function(B,Y,piA,piAX,gamma,index){
     ### bound for gamma not too far from 1
     ### using G(v); not really needed 

     n = length(Y)
     d = ncol(B)
     W = piA/piAX
     M = matrix(0,d,d)
     for(i in 1:n){
         b = as.matrix(B[i,])
         M = M + (W[i]/n)*(b %*% t(b))
         }
 
     alpha = t(B) %*%  as.matrix(Y*W)/n
     beta = solve(M) %*% alpha
 

     Minv = solve(M)

     deriv = rep(0,n)
     for(i in 1:n){
          b = as.matrix(B[i,])
          deriv[i] = Y[i]*W[i]*c((Minv %*% b)[index,])
          }


     deriv = mean(abs(deriv))
 
     low = beta[index] - log(gamma)*deriv
     up  = beta[index] + log(gamma)*deriv
     return(list(low=low,up=up))
     }






Homotopy = function(B,Y,piA,piAX,gamma,index){
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



     W = piA/piAX
     ng = length(gamma)
     n = length(Y)
     up = rep(0,ng)
     low = rep(0,ng)
     Vout = matrix(1,n,ng)

     ##upper bound
     V = diag(W)
     G = solve(t(B) %*%  V %*% B)
     beta = (G %*% t(B) %*% V %*% as.matrix(Y))
     up[1] = beta[index]
     d = rep(0,n)


     for(i in 1:n){
          ri = as.matrix(B[i,])
          tmp = sum((ri *beta))
          d[i] = ( (G %*% ri)[index]) * (Y[i] - tmp)*W[i]
          }
     for(j in 2:ng){
            t = gamma[j]/(1+gamma[j])
            q = quantile(d,t)
            alpha = rep(1/gamma[j],n)
            alpha[d >= q] = gamma[j]
            Vout[,j] = alpha
            V = diag(alpha*W)
            G = solve(t(B) %*%  V %*% B)
            beta = (G %*% t(B) %*% V %*% as.matrix(Y))
            up[j] = beta[index]
            for(i in 1:n){
               ri = as.matrix(B[i,])
               tmp = sum((ri*beta))
               d[i] = ((G %*% ri)[index]* (Y[i] - tmp))*W[i]
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
          d[i] = ((G %*% (si - tmp))[index])*W[i]
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
                 d[i] = ((G %*% (si - tmp))[index])*W[i]
                 }
            }

 
       return(list(low=low,up=up,Vout=Vout))
       }





Linear = function(B,Y,piA,piAX,gamma,index){
     ### bound on beta
     ### exact formula
     ### using G(v) and E[V]=1
     B = as.matrix(B)
     n = length(Y)
     ng = length(gamma)
     W = piA/piAX
     d = ncol(B)
     e = rep(0,d)
     e[index] = 1
     e = as.matrix(e)

     M = matrix(0,d,d)
     for(i in 1:n){
         b = as.matrix(B[i,])
         M = M + (W[i]/n)*(b %*% t(b))
         }
     Minv = solve(M)

     f = rep(0,n)
     for(i in 1:n){
          b = as.matrix(B[i,]) 
          f[i] = Y[i]*W[i]*c(t(e) %*% Minv %*% b)
          }
     Vout = matrix(1,n,ng)
     up  = rep(0,ng)
     low = rep(0,ng)
     for(i in 1:ng){
          q = quantile(f,gamma[i]/(1+gamma[i]))
          v = rep(1/gamma[i],n)
          v[f>q] = gamma[i]
          up[i] = mean(f*v)
          Vout[,i] = v
          q = quantile(f,1/(1+gamma[i]))
          v = rep(1/gamma[i],n)
          v[f<q] = gamma[i]
          low[i] = mean(f*v)

          }


     return(list(low=low,up=up,Vout=Vout))
     }




     

set.seed(977)

# pdf("F_Versus_G.pdf",width=6,height=2)
png("./results/F_Versus_G.png",width=15,height=5, units = "cm", res = 313)
par(mfrow=c(1,2))
par(mar=c(4.5,5.1,2,1),mex=.7)
par(cex.lab=.9,cex.axis=.7)

##Example 1: F(v) bound is tighter
n = 1000
beta = 3
X = rnorm(n)
A = X + rnorm(n)
Y = beta*A + 2*X + rnorm(n)
piA = dnorm(A,0,2)
piAX = dnorm(A,X,1)
W = piA/piAX
B = cbind(rep(1,n),A)
gamma = seq(1,3,length=20)
index = 2


out = Homotopy(B,Y,piA,piAX,gamma,index)
matplot(gamma,cbind(out$low,out$up),type="l",lty=c(1,1),col=c(1,1),
 xlab=expression(gamma),ylab=expression(beta),lwd=c(2,2), ylim=c(0,8))
mtext("(a)",side=3,line=.5)
tmp = Linear(B,Y,piA,piAX,gamma,index)
lines(gamma,tmp$low,lwd=2,col="darkcyan")
lines(gamma,tmp$up,lwd=2,col="darkcyan")
tmp = Local(B,Y,piA,piAX,gamma,index)
lines(gamma,tmp$low,lwd=2,col="red",lty=3)
lines(gamma,tmp$up,lwd=2,col="red",lty=3)



##Example 2: G(v) bound is tighter
n = 1000
gamma = seq(1,3,length=20)
piA = rep(1,n)
piAX = rep(1,n)
index = 1
U = seq(0.5,1,length=n)
A = 3-U
Y = 5*U
pi = .5
Y = pi*Y + (1-pi)*rnorm(n)
B = as.matrix(A)



out = Homotopy(B,Y,piA,piAX,gamma,index)
matplot(gamma,cbind(out$low,out$up),type="l",lty=c(1,1),col=c(1,1),
 xlab=expression(gamma),ylab="",lwd=c(2,2))
mtext("(b)",side=3,line=.5)
tmp = Linear(B,Y,piA,piAX,gamma,index)
lines(gamma,tmp$low,lwd=2,col="darkcyan")
lines(gamma,tmp$up,lwd=2,col="darkcyan")
tmp = Local(B,Y,piA,piAX,gamma,index)
lines(gamma,tmp$low,lwd=2,col="red",lty=3)
lines(gamma,tmp$up,lwd=2,col="red",lty=3)

dev.off()


