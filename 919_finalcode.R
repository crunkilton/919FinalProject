##919 Final
##Cody Crunkilton
##May 2017

setwd("C:/Users/Cody/Dropbox/1school2016_7/945/replications/")
library(rjags)
library(runjags)

bm <- read.csv("bm.csv")

###########################the model

zeroinflmodel1<-"
model {
for (i in 1:N) {

ystar[i]~dbern(p[i])
logit(p[i])<-inprod(x[i,],beta_count) 

y[i] ~ dnegbin(mu[i], r*ystar[i]) 
logit(mu[i]) <- inprod(x[i,],beta_casualties) 
}

#priors--logit part
beta_count[1:K_beta_count]~dmnorm(beta_count_mean, beta_count_precision)

#priors--neg binom part
beta_casualties[1:K_beta_casualties] ~ dmnorm(beta_casualties_mean, beta_casualties_precision)
r ~ dunif(0,1000)



}
"

## Variables to feed in:
x1 <- bm$twonukedyad
x2 <- bm$onenukedyad
x3 <- bm$logDistance
x4 <- bm$Contiguity
x5 <- bm$logCapabilityRatio
x6 <- bm$Ally
x7 <- bm$SmlDemocracy
x8 <- bm$NIGOs

x <- cbind(1,x1,x2,x3,x4,x5,x6,x7,x8)


prec <- diag(.0000000000001,9,9)
prec[2,2] <- .2
prec

beta_count_mean <- c(0,-2,0,0,0,0,0,0,0)
#beta_count_mean=rep(0,9) for flat model
beta_count_precision=prec #changed this for flat model as well
beta_casualties_mean <- c(0,-2,0,0,0,0,0,0,0)
#beta_casualties_mean=rep(0,9) for flat model
beta_casualties_precision=prec#changed this for flat model as well
K_beta_count=length(beta_count_mean)
K_beta_casualties=length(beta_casualties_mean)


##### Running the model:

b_priors2 <- run.jags(zeroinflmodel1, monitor=c("beta_count","beta_casualties"),
                      data=list(y=bm$deaths,
                                ystar=ifelse(bm$deaths>0, 1, NA),
                                N=length(bm$deaths),
                                x=x,
                                beta_count_mean=beta_count_mean,
                                beta_count_precision=beta_count_precision,
                                beta_casualties_mean=beta_casualties_mean,
                                beta_casualties_precision=beta_casualties_precision,
                                K_beta_count=K_beta_count,
                                K_beta_casualties=K_beta_casualties),
                      adapt=100000,burnin=10000, sample=100000, thin=20, n.chains=3)


##results
summary(b_priors2)
gelman.diag(b_priors2)
plot(b_priors2, var="beta_count[1]",plot.type="density")
plot(b_priors2, var="beta_count[2]",plot.type="density")
plot(b_priors2, var="beta_casualties[1]",plot.type="density")
plot(b_priors2, var="beta_casualties[2]",plot.type="density")


####printing out results

df2 <- as.data.frame(summary(b_priors2))
df2[c(1:5)]
library(xtable)
xtable(df2[c(1:5)])

traceplot(as.mcmc(b_priors2))
