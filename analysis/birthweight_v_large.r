setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(nnet)
library(splines)
rm(list=ls())

data = read.csv("./data/birthweight_data_sensitivity_msm.csv",header=TRUE)

Propensity_sub = function(T,Xc,Xd,drop=0){
        ## Propensity score
        ## can also drop a variable
        n = length(T)
        d.Xc = ncol(Xc)
        knots = c(3,1,3,1,3,1,1)
        if((drop <= 7)&(drop >0)){
                knots = knots[-drop]
                Xc = Xc[,-drop]
                d.Xc = d.Xc - 1
        }
        if(drop > 7){
                Xd = Xd[,-drop]
        }
        R = Xd
        for (j in 1:d.Xc){
                R=cbind(R, bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
        }
        gps.fit = nnet::multinom(T ~ R, maxit=5000)
        gps.hat = gps.fit$fitted.values
        
        
        
        piAX = rep(0,n)
        for(i in 1:n){
                piAX[i] = gps.hat[i,(T[i]+1)]
        }
        return(piAX)
}

###############################################
##### Start analysis
###############################################

data = data[,-1]
Y = as.vector(data[,1])
T = as.vector(data[,2])
Xd = as.matrix(data[,3:46])
Xc = as.matrix(data[,47:53])
n = length(Y)
d.Xd = ncol(Xd)
d.Xc = ncol(Xc)
T.levels = sort(unique(T))
d.T = length(T.levels)

piAX = Propensity_sub(T,Xc,Xd,drop=0)
tmp = table(T)/n
piA = rep(0,n)
for(i in 1:n){
        piA[i] = tmp[T[i]+1]
}
W = piA/piAX

n = length(Y)
d = rep(0,n)
WW = diag(W)
A = T
one = rep(1,n)
B = cbind(one,A,A^2)
d = solve((t(B) %*% WW %*% B)) %*% t(B) %*% WW
gamma = seq(1,1.2,length=10)

Homotopy_subset_propensity = function(B,Y,piA,piAX,gamma,epsilon,index){
        ### This function bounds beta[index]
        ### under propensity ratio sensitivity for a vector of values gamma
        ### using homotopy: i.e. find local optimum for gamma_j
        ### starting at solution for gamma_{j-1}.
        ###
        ### Y(a) = sum_j beta_j b_j(a)
        ### B = [b1 ... bk] = design matrix
        ### Y = outcome
        ### A = treatment
        ### piA = (pi(A_1),...,pi(A_n)) marginal density of treatment
        ### piAX = (pi(A_1|X_1),...,pi(A_n|X_n)) propensity score
        
        
        W = piA/piAX
        ng = length(gamma)
        n = length(Y)
        up = rep(0,ng)
        low = rep(0,ng)
        Ident = diag(1,n,n)
        
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
                alpha = (1-epsilon)*rep(1,n) + epsilon*alpha
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
                alpha = (1-epsilon)*rep(1,n) + epsilon*alpha
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

### Subset analysis
png("./results/smokingeps.png", units = "cm", height=5, width=15, res=313)

par(mfrow = c(1, 2), mar = c(4.5, 5.1, 2, 1), mex = 0.7, cex.lab = 0.9, cex.axis = 0.7)

out1 = Homotopy(B, Y, piA, piAX, gamma, index = 2)
matplot(gamma, cbind(out1$low, out1$up), type="l", col = c(1,1), 
        xlab = expression(gamma), ylab = expression(beta[1]), lty=c(1,1), 
        lwd = c(3,3))
abline(h = out1$up[1], lty = 2)

out2 = Homotopy(B,Y,piA,piAX,gamma,index = 3)
matplot(gamma, cbind(out2$low, out2$up), type="l", col = c(1,1), 
        xlab = expression(gamma), ylab = expression(beta[2]), lty = c(1,1),
        lwd = c(3,3))
abline(h = out2$up[1], lty = 2)

out3 = Homotopy_subset_propensity(B, Y, piA, piAX, gamma, epsilon = 0.1, index = 2)
matplot(gamma, cbind(out1$low, out1$up), type = "l", col = c(1,1), xlab = expression(gamma),
        ylab = expression(beta[1]), lty = c(1,1), lwd = c(2,2))
mtext("(a)",side = 3, line = 0.5)
abline(h = out1$up[1], lty = 2)
lines(gamma,out3$up, lwd = 2,col = "lightgray")
lines(gamma,out3$low, lwd = 2,col = "lightgray")

out4 = Homotopy_subset_propensity(B, Y, piA, piAX, gamma, epsilon = 0.5, index = 2)
lines(gamma, out4$up, lwd = 2,col = "grey50")
lines(gamma, out4$low, lwd = 2,col = "grey50")


out5 = Homotopy_subset_propensity(B, Y, piA, piAX, gamma, epsilon = 0.1, index = 3)
matplot(gamma, cbind(out2$low,out2$up), type = "l", col = c(1,1), 
        xlab = expression(gamma), ylab = expression(beta[2]), lty = c(1,1),
        lwd = c(2,2))
mtext("(b)", side = 3,line = 0.5)
abline(h = out2$up[1],lty = 2)
lines(gamma, out5$up, lwd = 2,col = "lightgray")
lines(gamma, out5$low, lwd = 2,col = "lightgray")

out6 = Homotopy_subset_propensity(B, Y, piA, piAX, gamma, epsilon = 0.5, index = 3)
lines(gamma, out6$up, lwd = 2, col = "grey50")
lines(gamma, out6$low, lwd = 2, col = "grey50")

dev.off()