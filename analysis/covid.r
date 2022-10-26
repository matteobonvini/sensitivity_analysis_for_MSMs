setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##### MAIN PROGRAM
##outcome is death

## clear out old variables
rm(list=ls())

source("Functions.r")




#######Get data


df = read.csv(file="./data/covid_data_sensitivity_msm.csv",head=TRUE)
n = nrow(df)
state = df$state
remove = c("gu","pr","mp","as","vi")
n = nrow(df)
`%notin%` <- Negate(`%in%`)
I = (1:n)[state %notin% remove]
df = df[I,]
n = nrow(df)
date = df$date
date = as.Date(date)
I = (1:n)[date > "2020-02-15"]
df = df[I,]



state = df$state
date = df$date
deaths = df$deaths
deaths[deaths < 0] = 0
home = df$completely_home_prop
bars = df$bars_prop
restaurants = df$restaurants_prop
State = unique(state)
population = df$chr_population
population[is.na(population)] = 0
n = length(deaths)



### Weekly data
Y = NULL
A = NULL   ##stay at home
R = NULL   ##restautants
B = NULL   ##bars
P = NULL   ##population
for(i in 1:50){
     I = (1:n)[state == State[i]]
     Y = rbind(Y,we(deaths[I]))
     B = rbind(B,we(bars[I]))
     A = rbind(A,we(home[I]))
     R = rbind(R,we(restaurants[I]))
     P = rbind(P,we(population[I]))
     }
Population = P[,1]/7
A = A/7


delta=4


##############################################################################################################




r = nrow(A)
o = order(Population)
A = A[o,]
Y = Y[o,]
State = State[o]
Population = Population[o]
Lp = log(Population)



### Blip Model
aic = matrix(0,50,4)

out = Blip(A,Y,State,degree=1,delta)
B1 = out$beta
S1 = out$se
aic[,1] = out$aic

out = Blip(A,Y,State,degree=2,delta)
B2 = out$beta
S2 = out$se
out_Fit = out$Fit
out_time = out$time
out_Y = out$Y
aic[,2] = out$aic

out = Blip(A,Y,State,degree=3,delta)
B3 = out$beta
S3 = out$se
aic[,3] = out$aic

out = Blip(A,Y,State,degree=4,delta)
B4 = out$beta
S4 = out$se
aic[,4] = out$aic

Degree = apply(aic,1,which.min)



### Fit the MSM
gamma = seq(1,3,length=20)
degree = Degree
out = Fit_Bounds(A,Y,State,gamma,degree,delta=4)

png("./results/Bounds_Beta_Covid.png", units = "cm", height=5, width=15, res=313)
par(mfrow=c(1,4))
par(mar=c(4.5,5.1,2,1),mex=.8)
par(cex.lab=1.7,cex.axis=1.5)
I = c(35,47,49,50)
for(i in I){
     matplot(gamma,cbind(out$low[i,],out$up[i,]),type="l",ylim=c(-10,1),xlab=expression(gamma),col=c(1,1),lty=c(1,1),lwd=2,
         ylab="",cex.axis=.9)
     title(State[i],cex.main=1.0)
     if(i == 35)mtext(expression(beta),side=2,line=2.5,cex=.8)
     abline(h=0,lty=2)
     lines(gamma,out$up_local[i,],col=2,lwd=2,lty=3)
     lines(gamma,out$low_local[i,],col=2,lwd=2,lty=3)
     }
dev.off()

### Bounds on g(a_t)
png("./results/Bounds_g.png", units = "cm", height=10, width=18, res=313)
par(mfrow=c(2,2))
par(oma=c(3,3,0,0), mar=c(3,3,2,1),mex=.8)
par(cex.lab=1.7,cex.axis=1.5)
gamma = c(1,1.5,2,3)
degree = Degree
out = Fit_Bounds_g(A,Y,State,gamma,degree,delta=4)
I = c(35,47,49,50)
for(i in I){
  matplot(1:40,cbind(exp(out$Low[[i]]),exp(out$Up[[i]]),out$Deaths[i,]),
   xlab="",ylab="",lty=1,type="n",cex.axis=.9,cex.lab=.9, main=State[i])
  polygon(c(1:40,40:1),exp(cbind(out$Low[[i]][,4],rev(out$Up[[i]][,4]))),col="grey100")
  polygon(c(1:40,40:1),exp(cbind(out$Low[[i]][,3],rev(out$Up[[i]][,3]))),col="grey80")
  polygon(c(1:40,40:1),exp(cbind(out$Low[[i]][,2],rev(out$Up[[i]][,2]))),col="grey60")
  lines(1:40,exp(cbind(out$Low[[i]][,1])))
  points(1:40,out$Deaths[i,],pch=19)
  }
mtext(side=1,line=1,"Weeks",outer=T,cex.lab=.9)
mtext(side=2,line=1,"Deaths",outer=T,cex.lab=.9)
dev.off()

