#### random effects for year of experiment!!
## build a better data collection
chl.DATA <- read.csv("chl.DATA.csv")
chl.DATA$seconds <- chl.DATA$Exp.time.decmins*60

xdatad <- seq(0,1800,length.out = 901)
xdatan <- c(0,1800)
ydatan <- c(0,0)

dats <- sapply(paste("D",1:8,sep=""),function(Dn){
  list(dosed = chl.DATA[chl.DATA$ExptNo == Dn & chl.DATA$Treatment == "dosed" ,],
       notdosed = data.frame(cbind(chl.DATA[chl.DATA$ExptNo == Dn & chl.DATA$Treatment == "none" ,]$seconds[1:181],matrix(chl.DATA[chl.DATA$ExptNo == Dn & chl.DATA$Treatment == "none" ,]$Dosed.chamber,nrow=181))))},simplify = F)

getting_switch<- function(se){
  res<-matrix(1,length(se),2)
  res[1,2] <- se[1] + 1
  for (i in 2:length(se))
  {
    res[i,1] <- se[i-1] + 1 
    res[i,2] <- se[i] + 1
  }
  return(res)
}

cons <- lapply(dats,function(x){
  if(ncol(x$notdosed)>2){
    return(apply(x$notdosed[,-1],2,getting_switch))}
  return(getting_switch(x$notdosed[,-1]))})

## sorting out the other cons
cons[-1]<-lapply(cons[-1],function(x){return(list(matrix(x[,1],ncol=2),matrix(x[,2],ncol=2),matrix(x[,3],ncol=2)))})
cons[[1]] <- list(cons[[1]])

### getting the transition matrix
ords <- matrix(c(1:5,2:6,1:6,2:6,1:5,1:6),ncol=2)
gettingTran <- function(chlor,beta1,beta2,beta3,kappa=1,finding_max=F)
{
  ret <- matrix(0,6,6)
  temp1 <- beta2 * chlor
  temp2 <- exp(c((beta1 + beta3 + temp1),(beta1[(5:1)] - temp1))) / kappa
  temp3 <- c(temp2[1],temp2[2] + temp2[6],temp2[3] + temp2[7],temp2[4] + temp2[8],temp2[5] + temp2[9],temp2[10])
  if(finding_max==T)
  {
    return(max(temp3))
  }
  ret[ords] <- c(temp2,1-temp3)
  #diag(ret) <- 1 - rowSums(ret)
  return(ret)
}
### 

### For 1) I will sample a path and then do linear interpolation
### ## the linear interpolation
chlor_fun <- function(xnew, x, y){
  #if (any(x==xnew)){return(y[which(x==xnew)])}
  #n0 <- max(which(x < xnew))
  n0 <- floor(xnew/2) + 1 ## a very dirty fix
  n1 <- n0 + 1
  x_0 <- x[n0] ; y0 <- y[n0]
  x_1 <- x[n1] ; y1 <- y[n1]
  return(((y1 - y0)/(x_1-x_0)) * (xnew - x_0) + y0)
}

prop_beta <- function(x,beta,COVS)
{
  beta[x,] + as.numeric(t(chol(COVS[,,x])) %*% as.matrix(rnorm(5,0,1)))
}

library(Matrix)

gettingTran1 <- function(chlor,beta1,beta2,beta3,kappa=1,finding_max=F)
{
  ret <- matrix(0,6,6)
  temp1 <- beta2 * chlor
  diag(ret[-6,-1]) <- exp(beta1 + beta3 + temp1)/kappa
  diag(ret[-1,-6]) <- exp(beta1[(5:1)] - temp1)/kappa
  if(finding_max==T)
  {
    return(max(rowSums(ret)))
  }
  diag(ret) <- 0 - rowSums(ret)
  return(ret)
}

getting_like <- function(s_t,beta1,beta2,beta3,sampled,seconds,xdata,ydata,new_r=T,rpot_switches=NA)
{
  t <- seconds[s_t-1]; ends <- seconds[s_t]; sta_state <- sampled[s_t-1]; end_state <- sampled[s_t] ### FIn the start time, start state, end time and end state of the time frame
  ### it works - now lets turn this into a function
  if (length(xdata) == 2)
  {
    kappa <- gettingTran(0,beta1,beta2,beta3,finding_max = T)
  }
  else{
    ## find_ min and max chlorine level in the time window
    min_t <- which(abs(xdata - t) < 0.05); max_t <- which(abs(xdata - ends) < 0.05);
    min_ct <- min(ydata[min_t : max_t]) ; max_ct <- max(ydata[min_t : max_t]);
    kappa <- max(gettingTran(min_ct,beta1,beta2,beta3,finding_max = T),gettingTran(max_ct,beta1,beta2,beta3,finding_max = T)) ## this is only a temporary fix - probably work for the whole data but will not in general
  }
  #browser()
  if (new_r==T){
    #if(s_t==36 & sta_state==5 & end_state==0){browser()}
    rpot_switches <- runif(1)
    num <- qpois(rpot_switches,(ends - t) * kappa)
    rpot_switches <- c(rpot_switches,runif(num))
    #pot_switch <- qunif(rpot_switch[-1],t,ends)
  }
  else{
    rpot_switches<-rpot_switches[[s_t-1]]
    num <- qpois(rpot_switches[1],(ends - t) * kappa)
    if (length(rpot_switches) < (num+1))
    {
      #browser()
      rpot_switches <- c(rpot_switches,runif(num - length(rpot_switches) + 1))
    }
  }
  pot_switch <- qunif(rpot_switches[2:(num+1)],t,ends)
  
  p.state <- t(as.matrix(rep(0,6))); p.state[sta_state+1] <- 1
  if(num==0){return(c(p.state[end_state+1],dpois(num,(ends - t) * kappa),list(rpot_switches)))}
  chl_ps <- sapply(pot_switch,chlor_fun,x=xdata,y=ydata)
  ## have to check to see if there are enough switches - not going to do this here as it seems that have loads of runs.
  ## getting the matricies
  tran.mat<- sapply(chl_ps,gettingTran,beta1=beta1,beta2=beta2,beta3=beta3,kappa=kappa,simplify=F)
  
  p.state <- Reduce(`%*%`,tran.mat,p.state)
  return(c(p.state[end_state+1],dpois(num,(ends - t) * kappa),list(rpot_switches)))
}


get_ll<- function(x,dats,beta1,beta2,beta3,xdatad,ydatad,xdata,ydata,mu,sigma,log_like,dosed=T,new_mov=T,pot_switches,siggy,cons)
{
  pot_switches <- pot_switches[[x]]
  if (dosed==F)
  {
    undosed_ll <- log_like[[1]][1,x]
  }
  #browser()
  temp<-sapply(2:181,getting_like,beta1=beta1[x,],beta2=beta2[x],beta3=beta3[x],sampled=dats[[x]]$dosed$Dosed.chamber,seconds=dats[[x]]$dosed$seconds,xdata=xdatad,ydata=ydatad[,x],new_r=new_mov,rpot_switches=pot_switches)
  pot_switches<- temp[3,]
  #pot_switches[[length(pot_switches)]] <- temp[3,]
  dosed_ll <- sum(log(unlist(temp[1,])))
  lpsd <- (log(unlist(temp[1,])))
  #dosed_ll<- sum(log(sapply(2:181,getting_like,beta1=beta1[x,],beta2=beta2,beta3=beta3[x],sampled=dats[[x]]$dosed$Dosed.chamber,seconds=dats[[x]]$dosed$seconds,xdata=xdatad,ydata=ydatad[,x])))
  if (dosed==T)
  {
    tran1<- log(expm(10 * gettingTran1(0,beta1=beta1[x,],beta2=beta2[x],beta3=beta3[x])))
    ll <- unlist(lapply(cons[[x]],function(y,tran1){sum(tran1[y])},tran1=tran1))
    undosed_ll <- sum(ll)
    
    
    # test <- sapply(1:(ncol(dats[[x]]$notdosed)-1),function(input,beta1,beta2,beta3,xdata,ydata,x,pot_switches)
    # {
    #   sapply(2:181,getting_like,beta1=beta1,beta2=beta2[x],beta3=beta3,sampled=dats[[x]]$notdosed[,input+1],seconds=(0:180)*10,xdata,ydata,new_r=new_mov,rpot_switches=pot_switches[[input]])
    # },beta1=beta1[x,],beta2=beta2,beta3=beta3[x],xdata=xdatan,ydata=ydatan,x=x,pot_switches=pot_switches)
    # 
    # sorting_mul <- apply(test,2,function(test)
    # {
    #   ll <- sum(log(unlist(matrix(test,nrow=3)[1,])))
    #   pll <- (log(unlist(matrix(test,nrow=3)[1,])))
    #   ps <- matrix(test,nrow=3)[3,]
    #   return(c(ll,list(pll),list(ps)))
    # })
    # pot_switches[1:(length(pot_switches)-1)] <- lapply(sorting_mul,function(x){x[[3]]})
    # undosed_ll<-unlist(lapply(sorting_mul,function(x){x[[1]]}))
    # lpsu  <-(lapply(sorting_mul,function(x){x[[2]]}))
  }
  ### priors
  #browser()
  return(c(undosed_ll,dosed_ll,list(lpsd),list(pot_switches)))
}

sorting_ll <- function(cl,dats,beta1,beta2,beta3,xdatad,ydatad,xdata,ydata,mu,sigma,dosed=T,new_mov=T,pot_switches,siggy,log_like=matrix(0,2,8),cons)
{
  log_like <- sapply(1:8,get_ll,dats=dats,beta1=beta1,beta2=beta2,beta3=beta3,xdatad=xdatad,ydatad=ydatad,xdata=xdata,ydata=ydata,mu=mu,sigma=sigma,log_like=log_like,dosed=dosed,siggy=siggy,pot_switches=pot_switches,new_mov=new_mov,cons)
  ll <- matrix(c(unlist(log_like[1,]),unlist(log_like[2,])),nrow=2,byrow=T)
  #log_like <- parSapply(cl,1:8,get_ll,dats=dats,beta1=beta1,beta2=beta2,beta3=beta3,xdatad=xdatad,ydatad=ydatad,xdata=xdata,ydata=ydata,mu=mu,sigma=sigma,log_like=log_like,dosed=dosed,siggy=siggy,pot_switches=pot_switches,new_mov=new_mov)
  
  #lps <- c(list(c(log_like[[1]][[3]],log_like[[1]][4])),lapply(log_like[-1],function(x){c(x[[5]],x[6])}))
  lps <- log_like[3,]
  pot_switches1 <- log_like[4,]
  #pot_switches1 <- c(list(log_like[[1]][5:6]),lapply(log_like[-1],function(x){x[7:10]}))
  return(list(ll = ll, lps =lps, pot_switches=pot_switches1))
}

priors<- function(beta1,mu,sigma,siggy){
  sum(dnorm(beta1,mu,sigma,log=T)) + sum(sapply(beta1,function(x){dnorm(max(x,0),0,siggy,log=T)}))  }

## update the hyperparameters - prior inv-gamma(2,2)
update.hyper<- function(beta,siggy=0.25)
{
  as <- 2 + length(beta)/2
  bs <- 2 + 1/2 * sum((beta - mean(beta))^2)
  sigma <- sqrt(1/rgamma(1,as,bs))
  picked <- F
  i <- 1
  while (picked == F)
  {
    mu <- rnorm(1,mean(beta),sigma/sqrt(length(beta)))
    if (runif(1) < dnorm(max(0,mu),0,siggy)/dnorm(0,0,siggy))
    {
      return(rbind(mu,sigma))
    }
    if (i > 10)
    {
      return(rbind(0,sigma))
    }
    i <- i +1
  }
}

source("gaussian_process.r")
ml <- sapply(1:8,plot.test,mu = 0.04,vars = vars, covs =covs,er=er)
ydatad <- matrix(unlist(lapply(ml,sim1path,mu=0.04,covs=covs,vars=vars,xn = seq(0,1800,length.out = 901))),ncol=8)
mu <- c(rep(-0.5,6),0)
sigma <- rep(0.5,7)
beta <- matrix(rnorm(40,mu[1:5],sigma[1:5]),8,5)
betaf <- rnorm(8,mu[6],sigma[6])
beta_c <- rnorm(8,mu[7],sigma[7])
siggy <- 1
#ntimes<-c(rep(5,5),rep(2,3))
#ntimes1<-c(rep(5,5),rep(2,3))
beta_s <- function(beta){
  beta[c(1,2,3,4,5,6,7,8),]
}

beta_f <- function(beta){
  beta[c(1,2,3,4,5,6,7,8)]
}

beta_cc <- function(beta){
  beta[c(1,2,3,4,5,6,7,8)]
}

library(parallel)

cl <- makeCluster(2)
clusterExport(cl, as.character(as.list(ls())))
clusterEvalQ(cl,library(DiceKriging))



log_like <- sorting_ll(cl,dats,beta_s(beta),beta_cc(beta_c),beta_f(betaf),xdatad,ydatad,xdatan,ydatan,mu,sigma,dosed=T,new_mov=T,pot_switches=c(list(rep(list(c()),2)),rep(list(rep(list(c()),4)),7)),siggy=siggy,cons=cons)
## could do this line in parallel by each experiment having it's own core.
library(abind)
s_ll <- array(log_like$ll,dim = c(2,8,1))
betas <- array(beta,dim = c(8,5,1)) ;beta_cs <- beta_c; betafs <- betaf
mus <- mu; sigmas <- sigma
ydatads <- array(ydatad,dim=c(901,8,1))
#load("one_per_year_COVS.Rdata")
#bcF <- c(0.05,0.05)
#bcP <- 0.75
#load("turning_pars.Rdata")
################################################################################################################


chgs<- matrix(0,nrow=2,ncol=8)

start_t<-Sys.time()

for (i in 1: 20000){
  ###### updating beta
  beta_ <- t(sapply(1:8,prop_beta,beta=beta,COVS=COVS))
  betaf_ <- rnorm(8,betaf,bcF)
  
  ll_<-sorting_ll(cl,dats,beta_s(beta_),beta_c,beta_f(betaf_),xdatad,ydatad,xdatan,ydatan,mu,sigma,dosed=T,new_mov=F,pot_switches=log_like$pot_switches,siggy=siggy,log_like=log_like[1:2],cons)
  
  
  tmp <- colSums(ll_$ll-log_like$ll)# + ll_$lps-log_like$lps)
  
  accepts <- which(log(runif(8)) <= tmp + apply(cbind(beta_,betaf_),1,priors,mu=mu,sigma=sigma,siggy=siggy) -  apply(cbind(beta,betaf),1,priors,mu=mu,sigma=sigma,siggy=siggy))
  if (length(accepts) > 0)
  {
    beta[accepts,] <- beta_[accepts,]
    betaf[accepts] <- betaf_[accepts]
    log_like$ll[,accepts] <- ll_$ll[,accepts]
    log_like$lps[accepts] <- ll_$lps[accepts]
    log_like$pot_switches[accepts] <- ll_$pot_switches[accepts]
  }
  ##### update beta_c - not at the moment (do this with the first step???)
  beta_c_ <- rnorm(8,beta_c,bcP)
  ll_<-sorting_ll(cl,dats,beta_s(beta),beta_c_,beta_f(betaf),xdatad,ydatad,xdatan,ydatan,mu,sigma,dosed=F,new_mov=F,pot_switches=log_like$pot_switches,siggy=siggy,log_like=log_like[1:2],cons)
  
  tmp <- colSums(ll_$ll-log_like$ll)# + ll_$lps-log_like$lps)
  
  accepts <- which(log(runif(8)) <= tmp + dnorm(beta_c_,mu[7],sigma[7],log=T) -  dnorm(beta_c,mu[7],sigma[7],log=T))
  if (length(accepts) > 0)
  {
    beta_c[accepts] <- beta_c_[accepts]
    log_like$ll[,accepts] <- ll_$ll[,accepts]
    log_like$lps[accepts] <- ll_$lps[accepts]
    log_like$pot_switches[accepts] <- ll_$pot_switches[accepts]
  }
  
  ydatad_ <- matrix(unlist(parLapply(cl,ml,sim1path,mu=0.04,covs=covs,vars=vars,xn = seq(0,1800,length.out = 901))),ncol=8)
  
  ll_<-sorting_ll(cl,dats,beta_s(beta),beta_c,beta_f(betaf),xdatad,ydatad_,xdatan,ydatan,mu,sigma,dosed=F,new_mov=F,pot_switches=log_like$pot_switches,siggy=siggy,log_like=log_like[1:2],cons)
  accepts <- (log(runif(16)) <= ll_$ll - log_like$ll)
  if (sum(accepts) > 0)
  {
    chg<-which(apply(accepts,2,function(x){any(x==T)})==T)
    for ( ii in chg)
    {
      # if(accepts[1,ii]==T)
      # {
      #   chgs[1,ii] <- chgs[1,ii] + 1
      #   log_like$pot_switches[[ii]][1:(length(log_like$pot_switches[[ii]])-1)] <- ll_$pot_switches[[ii]][1:(length(log_like$pot_switches[[ii]])-1)]
      #   log_like$ll[1,ii] <- ll_$ll[1,ii]
      # }
      if(accepts[2,ii]==T)
      {
        chgs[2,ii] <- chgs[2,ii] + 1
        log_like$pot_switches[[ii]] <- ll_$pot_switches[[ii]]
        log_like$ll[2,ii] <- ll_$ll[2,ii]
        log_like$lps[[ii]]<- ll_$lps[[ii]]
        ydatad[,ii] <- ydatad_[,ii]
      }
    }
  }
  #####
  
  ll_<-sorting_ll(cl,dats,beta_s(beta),beta_c,beta_f(betaf),xdatad,ydatad,xdatan,ydatan,mu,sigma,dosed=T,new_mov=T,pot_switches=log_like$pot_switches,siggy=siggy,log_like=log_like[1:2],cons)
  
  for (ii in 1:8)
  {
    accept <- which(log(runif(180)) < (ll_$lps[[ii]] - log_like$lps[[ii]] ))
    log_like$lps[[ii]][accept] <- ll_$lps[[ii]][accept]
    log_like$pot_switches[[ii]][accept] <- ll_$pot_switches[[ii]][accept]
    log_like$ll[2,ii] <- sum(unlist(log_like$lps[[ii]]))
  }
  
  
  ##### update the hyperparameters
  tmp<-apply(cbind(beta,betaf),2,update.hyper,siggy=siggy)
  mu[1:5] <- tmp[1,1:5] ; sigma[1:5] <- tmp[2,1:5];
  mu[6] <- tmp[1,6] ; sigma[6] <- tmp[2,6];
  tmp <- update.hyper(beta_c,siggy=1000)
  mu[7] <- tmp[1]; sigma[7] <- tmp[2]
  ### collecting the results
  betas <- abind(betas,beta); sigmas <- rbind(sigmas,sigma) ; mus <- rbind(mus,mu); betafs <- rbind(betafs,betaf)
  s_ll <- abind(s_ll,log_like$ll); ydatads <- abind(ydatads,ydatad)
  beta_cs <- rbind(beta_cs,beta_c)
  if (i/10 == floor(i/10)){
    if ((i / 10)%%10 == 0){
      plot.ts(t(betas[1,,]))
    }
    if ((i / 10)%%10 == 1){
      plot.ts(t(betas[2,,]))
    }
    if ((i / 10)%%10 == 2){
      plot.ts(t(betas[3,,]))
    }
    if ((i / 10)%%10 == 3){
      plot.ts(t(betas[4,,]))
    }
    if ((i / 10)%%10 == 4){
      plot.ts(t(betas[5,,]))
    }
    if ((i / 10)%%10 == 5){
      plot.ts(t(betas[6,,]))
    }
    if ((i / 10)%%10 == 6){
      plot.ts(t(betas[7,,]))
    }
    if ((i / 10)%%10 == 7){
      plot.ts(t(betas[8,,]))
    }
    if ((i / 10)%%10 == 8){
      plot.ts(beta_cs)
    }
    if ((i / 10)%%10 == 9){
      plot.ts(betafs)
    }
    print(colMeans(apply(cbind(beta_cs,betafs),2,diff)!=0))
    print(c(i,difftime(Sys.time() ,start_t,units="secs") / i))
    print(chgs)
  }
  #print(c(i,max(beta)))
}
