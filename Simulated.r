setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##Simulated example

rm(list=ls())
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

set.seed(1232)


png("./results/Simulated.png",width=15,height=5, units = "cm", res = 313)
par(mfrow=c(1,2))
par(mar=c(4.5,5.1,2,1),mex=.7)
par(cex.lab=.9,cex.axis=.7)

n = 100
beta = 3
A = rnorm(n)
Y = beta*A + rnorm(n)
W = rep(1,n)
piA = dnorm(A)
piAX = dnorm(A)
B = cbind(rep(1,n),A)
gamma = seq(1,3,length=20)
index = 2


out1 = Homotopy(B,Y,piA,piAX,gamma,index)
matplot(gamma,cbind(out1$low,out1$up),type="l",lty=c(1,1),col=c(1,1),xlab=expression(gamma),ylab=expression(beta),lwd=c(2,2),
ylim=c(2,4))
mtext("(a)",side=3,line=.5)
out2 = Local(B,Y,piA,piAX,gamma,index)
lines(gamma,out2$low,col="red",lwd=2,lty=3)
lines(gamma,out2$up,col="red",lwd=2,lty=3)

##Now do outcome sensitivity
m = 20
delta = seq(0,5,length=m)
tmp = lm(Y ~ A)
mu.hat = fitted(tmp)
d = ncol(B)
D = matrix(0,2,2)
C = matrix(0,2,1)
for(i in 1:n){
     b = as.matrix(B[i,])
     D = D + (W[i]/n)*(b %*% t(b))
     C = C + (W[i]/n)*mu.hat[i]*b
     }

Dinv = solve(D)
beta = (Dinv %*% C)[2]

s = Dinv[2,]
up = rep(0,m)
low = rep(0,m)
for(i in 1:n){
     tmp = sum(s*B[i,]*piA[i])
     up = up   + (delta*tmp*(tmp >0) - delta*tmp*(tmp <0))/n
     low = low + (delta*tmp*(tmp <0) - delta*tmp*(tmp >0))/n
     }
up = up + beta 
low = low + beta
matplot(delta,cbind(low,up),type="l",lty=c(1,1),col=c(1,1),lwd=2,xlab=expression(delta),ylab="",ylim=c(2,4))
mtext("(b)",side=3,line=.5)
abline(h=0,lty=2)
dev.off()
