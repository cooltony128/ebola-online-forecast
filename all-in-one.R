rm(list=ls())

#### Sarita code: Note - change wd so it matches with computer
##setwd("C:/Users/coolt/Documents/EBOLA") #set to working folder
data<-read.csv("2019-01-07_data.csv", header =T)


##############################################################################
colnames(data)<-c("date", "case")#rename
data<-data[!is.na(data$case),] #drop any rows with no cases
data<-data[data$case!=0,]

data$date<-strptime(data$date, format="%m/%d/%Y", tz="UTC")
data$date<-as.numeric(data$date)/86400 #convert to days. 86400 seconds in a day
data<-data[order(data$date),] #sort by time. crucial for zeroing out time
data$date<-data$date-data$date[1] #make first enty t=0

#initialize point process dataset
pp.data<-data.frame(matrix(ncol = 1, nrow = sum(data$case))); colnames(pp.data)<-"t"

set.seed(1111)
counter<-1#initialize
for(i in 1:nrow(data)) {#look at each row of "data"...
  for (j in 1:data$case[i]){#for as many cases there are in that row...
    pp.data[counter,1]<-runif(1,data$date[i]-1,data$date[i]) #jitter back as much as one day uniformly and add to pp.data
    counter<-counter+1
  }
}#instead of replacing each point in pp.data through a loop. probably more efficient to just grow a vector of
#... different lengths runifs with c(). But whatever this isn't huge data

pp.data<-data.frame(pp.data[order(pp.data$t),]); colnames(pp.data)<-"t" #sort by time. I don't know why sorting here removes data.frame format
pp.data$t<-pp.data$t-pp.data$t[1] #make first time be zero (pp.data$t[1] will some negative number)


write.csv(pp.data, file="ppdata.csv", row.names = F)
#################################################################################


NonparamETAS <- function(tpp, nbins, tot_time, num_iter, verbose = T) {
  supDist <- function (x, y) return (max (abs (x - y)))
  
  # Function that returns a nonparametric estimate of the background intensity mu and 
  # triggering function g(t-ti) for an ETAS model of a self exciting temporal point process.  
  # The triggering function g is estimated by a histogram.
  #
  # Args:
  #   tpp: a vector containing a realization of a self exciting point process
  #   nbins: the number of bins to be used in the histogram estimate of g(t-ti)
  #   tot_time: the total time over which the point process data was observed
  #   num_iter: the number of interations for the algorithm to run
  #   verbose: If TRUE (default), prints the sup of the distance between probabilty matrices 
  #            at each iteration
  #
  # Returns:
  #   mu: list element; estimate of background intensity
  #   g: list element; vector of length nbins where g[i] is the height of the histogram over 
  #      the interval (delta_t*(i-1), delta_t*(i)], where delta_t is the bin width  
  #
  N <- length(tpp)
  delta_t <- (tot_time + 1)/nbins  # bin width of histogram estimator g 
  g <- rep(0, nbins)  # intialize triggering function
  
  tdiffsmtx <- outer(tpp, tpp, FUN = "-")  # N by N matrix of times differences for tpp 
  tdiffs <- tdiffsmtx[lower.tri(tdiffsmtx, diag = F)]	#just the lower triangular part as a vector
  A <- list(0)  
  
  for(i in 1:nbins) {
    A[[i]] <- which(tdiffs > delta_t*(i-1) & tdiffs <= delta_t*i)
    if(length(A[[i]]) == 1) {
      if(is.na(A[[i]])) A[[i]] = 0
    } 
  }  
  
  bin_ind <- c(rep(0, length(unlist(A)))) #length here is # of pairs in all of A. should be same length as tdiffs
  for(i in 1:length(A)){ #length of A is just nbins
    bin_ind[A[[i]]] <- i #A[[i]] contains the address of the time differences falling into the i^th bin
  }
  
  # initialize matrix of probabilities
  Pmtx_new <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Pmtx_new[i,] <- rep(c(1/i), N)
  }
  P_new <- Pmtx_new[lower.tri(Pmtx_new)] #note: diag = F is the default
  
  # start convergence algorithm 
  for(k in 1:num_iter) {
    Pmtx_old <- Pmtx_new
    P_old <- P_new		
    
    # M-step.  
    g <- c() #g already initialized above
    for(j in 1:nbins) {
      g[j] <- sum(P_old[A[[j]]]) / (delta_t * sum(P_old))   
    }
    
    kp <- sum(P_old)/N #I added this
    
    # background rate is estimatated p_{jj}
    mu <- sum(diag(Pmtx_old)) / tot_time
    
    # take the probabilities from the M-step and put them back in their right place in the g matrix using the bin indicator 
    gmtx <- matrix(0, nrow = N, ncol = N)
    gmtx[lower.tri(gmtx)] <- g[bin_ind]*kp #I added kp
    
    # set the diagonal elements of the g matrix to the new background rate from the E-step 
    diag(gmtx) <- mu
    
    # E-step
    Pmtx_new <- matrix(0, nrow = N, ncol = N)  # Initialize new p matrix
    for(i in 1:N) {                            # For each event 
      Pmtx_new[i,] <- gmtx[i,] / sum(gmtx[i,]) # At each earthquake i, 
    }
    
    #cat('Est aftershocks:',sum(P_old)*delta_t)
    nafter <- sum(P_old)*delta_t
    P_new <- Pmtx_new[lower.tri(Pmtx_new)]
    
    if(k > 2 & k %%10 == 0) {
      sup <- supDist(Pmtx_old, Pmtx_new)
      if(verbose) {
        cat (
          "Iteration: ", k,
          "SupDist: ", formatC(sup, digits = 8, width = 12, forma = "f"),
          "\n")
      }
      if( sup < 1e-5 ) return(list(mu = mu, g = g, nbins = nbins, kp=kp, delta_t = delta_t, nafter = nafter,Pmtx_new=Pmtx_new)) 	
    }
  }
  return(list(mu = mu, g = g, nbins = nbins, kp=kp, delta_t = delta_t, Pmtx_new=Pmtx_new))	
}
## Fit the recursive model. 
## lambda = mu + SUM Ki beta exp(-beta (t-ti)), 
## where Ki = c / (lambda_i)^p. 
## theta = (mu, c, beta, p).

dat1<-read.csv("ppdata.csv", header=T)#jittered ppdata

date_start <- as.Date("2018-05-03") # define the beginning of the outbreak, May 3 2018, to be consistent with the paper

# Info on Variables:
# T: the number of days of the outbreak
# t: the times of each ebola case (random case distribution over each day)
# n: the total number of cases for the entire outbreak

# n is the total number of cases, so find the number of rows of the dataset
n <- length(dat1$t)

# number of days of the outbreak
T <- trunc(dat1$t[n]) + 1

t <- dat1$t # extract the t column

# Put into usable format for rest of the code:
z<-c()
z$n <- n
z$t <- t



## Before running this, you should already have a list t of case times, and 
## T should be the end of your time window, in days. The start of the time window = 0. 
## For example, just to play around with some fake data of length 578, you could just do 
##T = 300
#T=200







## This function computes the negative loglikelihood for the recursive model and 
## outputs a list with lambda and the negative loglikelihood. 
## It also prints 2 checks on the fit. Both should be about one. 
neglogl = function(theta){
  
  mu = theta[1]
  c = theta[2]
  beta = theta[3]
  p = theta[4]
  lam = rep(mu,n)
  ks = rep(1/mu,n)
  z = 0
  eps = 0.000001
  
  if(mu < eps) return(99999)
  if(c < eps) return (99999)
  if(beta < eps) return(99999)
  
  for(i in 2:n){
    
    lam[i] = mu + sum(ks[1:(i-1)]*beta * exp(-beta*(t[i]-t[1:(i-1)])))
    
    if(lam[i] < eps) z = 2
    else ks[i] = c*lam[i]^(-p)
    
  }
  
  
  if(z > 1) return(99999)
  
  sumloglam = sum(log(lam))
 
  intlam = mu*T*1 + sum(ks)

  cat(theta,". ")
  cat("\n loglik is ",sumloglam-intlam,".\n Int/n = ",intlam/n,".")
  cat("sum of (1/lam)/T is ",sum(1/lam)/T,".\n\n")
  
  xlam = list(nlog = intlam-sumloglam, lam=lam)
  
}


## The next two lines need to be done in terminal using R. Make sure measles.c is in your 
## working R directory.  
neglogc = function(theta1){
  
  t1 <<- t1 + 1
  
  if(min(c(theta1[1],theta1[2],theta1[3])) < 0.00000001) return(9e20)
  else
    a3 = neglogl(theta1)$nlog 
  
  if(p1 > 1){
    cat("\n iter = ", t1, " ")
  }
  
  a3
}



getlam = function(theta1){
  
  a1 = neglogl(theta1)$lam
  a1
  
}


t1 = 0
lam = rep(n/T,n)
p1 = 2
eps = .00000001

theta = c(.001,.03,.01,.5)

b2 = optim(theta,neglogc,hessian=T) #it stopped at iteration 511

theta = b2$par




#sqrt(diag(solve(b2$hess))) ## SEs 
## Sometimes you get Error in solve.default(b2$hess): system is computationally singular. 
## You have to get SEs by simulation in this case. I have code for that elsewhere. 


## This is to plot the resulting lambda. 
x = 1:T
n1 = length(x)
mu = theta[1]
c = theta[2]
beta = theta[3]
p = theta[4]
a1 = rep(mu,n1) ## a1[i] will be lambda at time x[i], whether or not a pt is on day x[i]. 
lam = rep(mu,n)
ks = rep(1/mu,n)

for(i in 2:n){
  lam[i] = mu + sum(ks[1:(i-1)]*beta * exp(-beta*(t[i]-t[1:(i-1)])))
  ks[i] = c*lam[i]^(-p)
}

for(i in 1:n1) if(x[i]>t[1]) a1[i] = mu + sum(ks[t<x[i]]*beta * exp(-beta*(x[i]-t[t<x[i]])))

y = rep(0,T)
for(i in 1:T) y[i] = sum((t<i) & t>= (i-1))

plot(x,a1,type="l",col=grey(0),xlab="time (days)",ylab="cases (pts/day)",ylim=c(0,20))
lines(x,y,col=grey(.5)) 
legend("topright",lty=c(1,1),col=c(grey(.5),grey(0)),legend=c("observed","forecast")) 





## Output the parameters of the model:
date_last <- date_start + (T-1)

name <- paste(as.character(date_last), "_recursive_params.txt", sep = '')
sink(file = name)
cat("mu:", theta[1],'\n')
cat("c:", theta[2],'\n')
cat("beta:", theta[3],'\n')
cat("p:", theta[4],'\n')
#cat("seed:", seed)
sink(NULL)

#correctedfitrecursive copy.txt
#Displaying correctedfitrecursive copy.txt.








mu_ = theta[1]
c_ = theta[2]
beta_ = theta[3]
p_ = theta[4]


#nforcst is the total number of cases
##do a simulation of cases and get tnew(randomized timestamps) and nforcst(not sure how to do it)

# Sarita: added seed to match the Hawkes code
seed <- 928371
set.seed(seed)
##########
## This is the part added by Rick 
ite <- c() #this is a list that saves tkeep function
for(j in 1:50){ #here we do 50 iterations
  Trange = 21
  ceil = 900 
  n1 = rpois(1,ceil*Trange)
  tcand = T+sort(runif(n1)*Trange)
  
  ## I'm supposing you already have n points up to time T. 
  ## We have to sum over these and the new points to get lambda for the new points. 
  tkeep = z$t
  Kkeep = c_ / (mu_ ^ p_) ## K for the first point in the dataset. 
  for(i in 2:n){
    
    lam = mu_ + sum(Kkeep[1:(i-1)] * beta_ * exp(-beta_ * (tkeep[i] - tkeep[1:(i-1)])))
    Kkeep[i] = c_ / (lam ^ p_)
    
  }
  ## Now go through the new candidate points. 
  m = n
  for(i in 1:n1){
    lam = mu_ + sum(Kkeep[1:m] * beta_ * exp(-beta_ * (tcand[i]-tkeep[1:m])))
    if(lam > ceil) cat("problem. lam = ", lam, "\n")
    if(runif(1) < lam/ceil){
      m = m+1
      tkeep[m] = tcand[i]
      Kkeep[m] = c_ / lam^p_
    }
  }
  
  ite[[j]] <- tkeep[tkeep>(T+1)]
  # all of outbreak
}
##this gives you one simulation
plot(c(0,T+Trange),c(0,1),xlab="t",ylab="",type="n")
points(tkeep,runif(m),pch=".") 
seven <- c() #list that stores tkeep in 7 days
fourteen <- c() #list that stores tkeep in 14 days
twentyone <- ite #list that stores tkeep in 21 days
for(i in 1:50){
  tkeeptemp <- ite[[i]]
  seven[[i]] <- tkeeptemp[tkeeptemp<(T+8)]
  fourteen[[i]] <- tkeeptemp[tkeeptemp<(T+15)]
}
sevendaytotalbeforeaverage <- 0
fourteendaytotalbeforeaverage <- 0
twentyonedaytotalbeforeaverage <- 0
x=rep(0,50)
y=rep(0,50)
z=rep(0,50)
for(i in 1:50){
  x[i] = length(seven[[i]])
  y[i] = length(fourteen[[i]])
  z[i] = length(ite[[i]])
  sevendaytotalbeforeaverage <- sevendaytotalbeforeaverage+ length(seven[[i]])
  fourteendaytotalbeforeaverage <- fourteendaytotalbeforeaverage + length(fourteen[[i]])
  twentyonedaytotalbeforeaverage <- twentyonedaytotalbeforeaverage + length(ite[[i]])
}

sevendaytotal <- sevendaytotalbeforeaverage / 50
fourteendaytotal <- fourteendaytotalbeforeaverage /50
twentyonedaytotal <- twentyonedaytotalbeforeaverage /50
sevendaytotal
fourteendaytotal
twentyonedaytotal
# forecast only
# plot(c(254+1,254+Trange),c(0,1),xlab="t",ylab="",type="n")
# for(i in 1:100){
#   l = length(ite[[i]])
#   points(ite[[i]],runif(m)[(m-l+1):m],pch=".") 
# }
# standard deviation
name2 <- paste(as.character(date_last), "_recursive_projection.txt", sep = '')
sink(file = name2)
cat("seven-day forecast:", round(sevendaytotal,1),'\n')
cat("fourteen-day forecast:", round(fourteendaytotal,1),'\n')
cat("twentyone-day forecast:", round(twentyonedaytotal,1),'\n')
cat("standard deviation of 7 day forecast:", round(sd(x),1),'\n')
cat("standard deviation of 14 day forecast:", round(sd(y),1),'\n')
cat("standard deviation of 21 day forecast:", round(sd(z),1),'\n')
cat("seed for forecasting:", seed)
sink(NULL)

dat1<-read.csv("ppdata.csv", header=T)#jittered ppdata


#### Sarita's code: added the start date for the output 
date_start <- as.Date("2018-05-03") # define the beginning of the outbreak, May 3 2018, to be consistent with the paper
n <- length(dat1$t) # number of total cases
T <- trunc(dat1$t[n]) + 1 # number of days of the outbreak
date_last <- date_start + (T-1)
#### end of Sarita's code


tpp<-as.numeric(dat1[,1])#put it as vector
tot_time<-max(tpp)-min(tpp)
fit <- NonparamETAS(tpp = tpp, nbins = (tot_time+1), tot_time = tot_time, num_iter = 1500)

binhghts<-fit$g[1:(max(which(fit$g>1e-4))+1)] #max(...) is the furthest out binhght greater than 1e-4. 
binhghts<-c(binhghts,1e-4)#Look a day ahead to let tail decay in smooth.
midpoints<-seq(fit$delta_t/2,by=fit$delta_t,length.out = length(binhghts))
smoothed<-ksmooth(midpoints,binhghts,"normal",
                  bandwidth =4,range.x=c(0,max(midpoints)+1),n.points=200)

#making histogram bin heights from the smoothed curve for simulation
new_binwidth<-diff(smoothed$x)[1] #diff(...) is the binwidth of our 'fine' hist from ksmoothing.
new_midpoints<-(smoothed$x-new_binwidth/2)[-1] # [-1] is to discard the negative midpoint
new_binhghts<-c()
for(i in 2:length(smoothed$y)){
  new_binhghts[i-1]<-mean(c(smoothed$y[i],smoothed$y[i-1]))
}

#Forecasting
seed <- 928371 # added by Sarita for output later
set.seed(928371)
castlen<-21#number of days we are forecasting
predictions<-c() #will be matrix of predictions, each row is one prediction
for (k in 1:1000){
  
  bground<-c()
  for (j in 1:length(tpp)){#for each point in the actual data prior to the period we are doing prediction on. get 1st gen children
    #numb of 1st gen children is poisson(Kp). where they are in time is by discrete prob weights given by new_binhghts (multinomial)
    bground<-c(bground,which(rmultinom(rpois(1,fit$kp),1,new_binhghts)==1, arr.ind = T)[,1]*new_binwidth-new_binwidth/2+tpp[j]) #-new_binwidth/2 bc want to simulate on the midpoint
  }
  
  bground<-c(bground,runif(n=rpois(1,fit$mu * castlen), min=max(tpp), max=max(tpp)+castlen)) #exogenous points from Mu
  bground<-subset(bground, bground<=(max(tpp)+castlen) & bground>=max(tpp)) #take only ones which fall in period to be forecasted
  
  predpoints<-bground  #store what points we have so far
  while(length(bground)>0){ #cascade of aftershocks until there is no more aftershocks
    
    aftshocks<-c() #reset aftershocks
    for (j in 1:length(bground)){
      aftshocks<-c(aftshocks,which(rmultinom(rpois(1,fit$kp),1,new_binhghts)==1, arr.ind = T)[,1]*new_binwidth-new_binwidth/2+bground[j])
    }
    newpoints<-subset(aftshocks, aftshocks<=(max(tpp)+castlen) & aftshocks>=max(tpp)) #take only ones which fall in the prediction zone
    predpoints<-c(predpoints,newpoints) #after shocks in the interval
    bground<-newpoints
  }
  
  predpoints<-sort(predpoints) #predicted points
  
  cumsum<-c()
  i<-1
  for (j in max(tpp):(max(tpp)+castlen)){
    cumsum[i]<-length(predpoints[predpoints<=j])
    i<-i+1
  }
  predictions<-rbind(predictions, cumsum)
}#end forecast simulation

#average of 1000 simulation is my forecast
avgprediction<-apply(predictions,2, mean)
#diff(avgprediction)
prediction75<-apply(predictions,2, function(x) quantile(x,0.75))
prediction25<-apply(predictions,2, function(x) quantile(x,0.25))

act_cumsum<-c() #actual cumsum of all of tpp
k<-1
for (i in 1:(max(tpp)+1)){
  act_cumsum[k]<-length(tpp[tpp<=i])
  k<-k+1
}
####################################################################################
#####################################################################################


# #Bar Plot of triggering
# barplot(fit$g, fit$delta_t, xlim=c(0,21), ylim=c(0,0.5), space=0, xlab = "Days from Infection", ylab = "Density",
#         main="Triggering Density")
# axis(1, at=seq(0,21))

# #plot of smoothed density
# plot(smoothed, type="l", ylab="Density",xlab="Days Since Reported Infection",main="", xaxt="n" )
# axis(side=1,at=seq(0,18,by=1))

fit$mu #background rate
fit$kp #productivity

#1 Week, 2 Week and 3 Week forecast
one_2_3_wk_preds<-c( mean(predictions[,7+1]),mean(predictions[,14+1]), mean(predictions[,21+1]))
one_2_3_wk_stdevs<-c( sd(predictions[,7+1]),sd(predictions[,14+1]), sd(predictions[,21+1]))

# #forecasting into the unknown
# plot(avgprediction+length(tpp), type="l", lwd=3, col="red",
#      xlab="Days forecasted", ylab="Cumulative cases", ylim=c(length(tpp),length(tpp)+max(prediction75)), xaxt="n",
#      main="1000 simulations starting XX/XX/XXXX")
# axis(1, at=seq(1,castlen+1), labels=seq(0,castlen))
# lines(prediction25+length(tpp))
# lines(prediction75+length(tpp))

# #Plot used in paper #THIS PLOT TAKES A FEW SECONDS
# plot(c(act_cumsum, (length(tpp)+avgprediction)[-1]), 
#      xaxt='n',col="red",type="l",ylab="Cumulative count", xlab="Days Since Initial Outbreak",
#      ylim=c(0,length(tpp)+max(predictions)),
#      main="")
# for (i in 1:nrow(predictions)){
#   lines(c(act_cumsum, (length(tpp)+predictions[i,])[-1]), col=rgb(0,0,1,0.1))
# }
# lines(c(act_cumsum, (length(tpp)+avgprediction)[-1]), col="red")
# lines(c(act_cumsum), lwd=2)
# abline(v=round(max(tpp)+1))
# axis(side=1, at=seq(0, 275, by=5))
# 
# #daily forecasted
# plot(diff(avgprediction),col="red",type="l",ylab="Forecasted Cases",main="Daily new cases forecasted")


#########################################################################################################
# Added to original code for online forecasting
# Sarita and Tony
# output data

## output parameters

name1 <- paste(as.character(date_last), "_hawkes_params.txt", sep = '')
sink(file = name1)
cat("mu:", fit$mu,'\n')
cat("k:", fit$kp,'\n')
cat("bin heights:", binhghts, '\n' )
sink(NULL)

## output forecasts
# find the median instead of the mean (as used previously)
one_2_3_wk_preds_med<-c( median(predictions[,7+1]),median(predictions[,14+1]), median(predictions[,21+1]))

name2 <- paste(as.character(date_last), "_hawkes_projection.txt", sep = '')
sink(file = name2)
cat("seven-day forecast:", one_2_3_wk_preds_med[1],'\n')
cat("fourteen-day forecast:", one_2_3_wk_preds_med[2],'\n')
cat("twentyone-day forecast:", one_2_3_wk_preds_med[3],'\n')
cat("standard deviation of 7 day forecast:", round(one_2_3_wk_stdevs[1], 1),'\n')
cat("standard deviation of 14 day forecast:", round(one_2_3_wk_stdevs[2], 1),'\n')
cat("standard deviation of 21 day forecast:", round(one_2_3_wk_stdevs[3], 1),'\n')
cat("seed for forecasting:", seed)
sink(NULL)





