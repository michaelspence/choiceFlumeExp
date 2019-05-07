### getting the chlorine levels for all time
## for 1
sampled1d <- chl.DATA[chl.DATA$ExptNo == "D1" & chl.DATA$Treatment == "dosed" ,]
sampled1n <- chl.DATA[chl.DATA$ExptNo == "D1" & chl.DATA$Treatment == "none" ,]
xdata1d <- sampled1d[complete.cases(sampled1d),]$Exp.time.decmins * 60  
ydata1d <- sampled1d[complete.cases(sampled1d),]$CHLOR    
xdata1n <- c(0,max(sampled1n$Exp.time.decmins)*60)
ydata1n <- c(0,0)

## for 2
sampled2d <- chl.DATA[chl.DATA$ExptNo == "D2" & chl.DATA$Treatment == "dosed" ,]
sampled2n <- chl.DATA[chl.DATA$ExptNo == "D2" & chl.DATA$Treatment == "none" ,]
xdata2d <- sampled2d[complete.cases(sampled2d),]$Exp.time.decmins * 60  
ydata2d <- sampled2d[complete.cases(sampled2d),]$CHLOR    
xdata2n <- c(0,max(sampled2n$Exp.time.decmins)*60)
ydata2n <- c(0,0)

## for 3
sampled3d <- chl.DATA[chl.DATA$ExptNo == "D3" & chl.DATA$Treatment == "dosed" ,]
sampled3n <- chl.DATA[chl.DATA$ExptNo == "D3" & chl.DATA$Treatment == "none" ,]
xdata3d <- sampled3d[complete.cases(sampled3d),]$Exp.time.decmins * 60  
ydata3d <- sampled3d[complete.cases(sampled3d),]$CHLOR    
xdata3n <- c(0,max(sampled3n$Exp.time.decmins)*60)
ydata3n <- c(0,0)

## for 4
sampled4d <- chl.DATA[chl.DATA$ExptNo == "D4" & chl.DATA$Treatment == "dosed" ,]
sampled4n <- chl.DATA[chl.DATA$ExptNo == "D4" & chl.DATA$Treatment == "none" ,]
xdata4d <- sampled4d[complete.cases(sampled4d),]$Exp.time.decmins * 60  
ydata4d <- sampled4d[complete.cases(sampled4d),]$CHLOR    
xdata4n <- c(0,max(sampled4n$Exp.time.decmins)*60)
ydata4n <- c(0,0)

## for 5
sampled5d <- chl.DATA[chl.DATA$ExptNo == "D5" & chl.DATA$Treatment == "dosed" ,]
sampled5n <- chl.DATA[chl.DATA$ExptNo == "D5" & chl.DATA$Treatment == "none" ,]
xdata5d <- sampled5d[complete.cases(sampled5d),]$Exp.time.decmins * 60  
ydata5d <- sampled5d[complete.cases(sampled5d),]$CHLOR    
xdata5n <- c(0,max(sampled5n$Exp.time.decmins)*60)
ydata5n <- c(0,0)

## for 6
sampled6d <- chl.DATA[chl.DATA$ExptNo == "D6" & chl.DATA$Treatment == "dosed" ,]
sampled6n <- chl.DATA[chl.DATA$ExptNo == "D6" & chl.DATA$Treatment == "none" ,]
xdata6d <- sampled6d[complete.cases(sampled6d),]$Exp.time.decmins * 60  
ydata6d <- sampled6d[complete.cases(sampled6d),]$CHLOR    
xdata6n <- c(0,max(sampled6n$Exp.time.decmins)*60)
ydata6n <- c(0,0)

## for 7
sampled7d <- chl.DATA[chl.DATA$ExptNo == "D7" & chl.DATA$Treatment == "dosed" ,]
sampled7n <- chl.DATA[chl.DATA$ExptNo == "D7" & chl.DATA$Treatment == "none" ,]
xdata7d <- sampled7d[complete.cases(sampled7d),]$Exp.time.decmins * 60  
ydata7d <- sampled7d[complete.cases(sampled7d),]$CHLOR    
xdata7n <- c(0,max(sampled5n$Exp.time.decmins)*60)
ydata7n <- c(0,0)

## for 8
sampled8d <- chl.DATA[chl.DATA$ExptNo == "D8" & chl.DATA$Treatment == "dosed" ,]
sampled8n <- chl.DATA[chl.DATA$ExptNo == "D8" & chl.DATA$Treatment == "none" ,]
xdata8d <- sampled8d[complete.cases(sampled8d),]$Exp.time.decmins * 60  
ydata8d <- sampled8d[complete.cases(sampled8d),]$CHLOR   


chlor <- cbind(ydata1d,ydata2d,ydata3d,ydata4d,ydata5d,ydata6d,ydata7d,ydata8d)
x <- c(0,300,600,900,1200,1500,1800)
plot(x,seq(0,0.1,length.out = length(x)),type='n',ylab='Chlorine',xlab='Time')
for(i in 1:8){
  points(x,chlor[,i],col=i)
}

## load the package for gaussian processes
library(DiceKriging)
?DiceKriging

dat <- data.frame(cbind(x,chlor))
# ml <- km(design=data.frame(X1=dat[,1]),response = dat[,2],nugget=0.0001)
# #nugget.estim=FALSE, noise.var=NULL
# x.n <- seq(0,1800,length.out = 1000)
# test <- predict(ml,newdata=data.frame(X1=x.n),"UK")
# plot(x,seq(0,0.1,length.out = length(x)),type='n',ylab='Chlorine',xlab='Time')
# points(x,chlor[,1],col=1)
# lines(x.n,test$mean,lty=2)
# lines(x.n,test$upper95,lty=3)
# lines(x.n,test$lower95,lty=3)

##### going to fit the model in stan and then use the fitted model with DiceKriging to simulate possible values.

# library(rstan)
# ## first test with one data set
# dat.stan <- list(N=length(x),M=8,y=t(chlor),sigma=0.0001)
# fit<- stan("gp.stan",data=dat.stan,iter = 2000,chains = 1)
# fit.ex<- stan("exponential.stan",data=dat.stan,iter = 2000,chains = 1)
# 
# ## stan doesn't really work. Possibly a trial and error using 
# ml <- km(design=data.frame(X1=dat[,1]),response = dat[,2],nugget=0.0001)
## create an object
mu <- 0.04;vars <- 0.01^2; covs <- 2000;er <- 0.0025^2
plot.test <- function(x,mu,vars,covs,er,xn = seq(0,1800,length.out = 1000))
{
  ml1 <- km(design=data.frame(X1=dat[,1]),response = dat[,x+1],nugget=er,coef.trend = mu,coef.var = vars,coef.cov = covs)
  #ml1 <- km(design=data.frame(X1=dat[,1]),response = dat[,x+1],coef.trend = mu,coef.var = vars,coef.cov = covs)
  test <- predict(ml1,newdata=data.frame(X1=xn),"UK")
  plot(dat[,1],seq(0,0.1,length.out = length(dat[,1])),type='n',ylab='Chlorine',xlab='Time')
  points(dat[,1],chlor[,x],pch=20)
  lines(xn,test$mean,lty=2)
  lines(xn,test$upper95,lty=3)
  lines(xn,test$lower95,lty=3)
  return(ml1)
}
#par(mfrow=c(4,2))
#blank<-sapply(1:8,plot.test,mu = mu,vars = vars, covs =covs,er=er)


### simulating paths

sim1path <- function(ml,mu,covs,vars,xn = seq(0,1800,length.out = 1000))
{
  post <- simulate(ml,nsim = 1,cond=T)
  ml1 <- km(design=data.frame(X1=dat[,1]),response = t(post),coef.trend = mu,coef.var = vars,coef.cov = covs)  

#plot(dat[,1],seq(0,0.1,length.out = length(dat[,1])),type='n',ylab='Chlorine',xlab='Time')
# points(dat[,1],chlor[,x],pch=20)
# lines(x.n,test$mean,lty=2)
# lines(x.n,test$upper95,lty=3)
# lines(x.n,test$lower95,lty=3)
  sim.dat <- simulate(ml1,nsim = 1,newdata=data.frame(X1=xn),cond=T,nugget.sim=1e-15)
  #lines(x.n,sim.dat[1,],lty=4)
  return(sim.dat)
}
par(mfrow=c(4,2))
for(i in 1:8)
{
  ml <- plot.test(i,mu = mu,vars = vars, covs =covs,er=er)
  sim<-sim1path(ml,mu,covs,vars,xn = seq(0,1800,length.out = 901))
  lines(seq(0,1800,length.out = 901),sim)
}

