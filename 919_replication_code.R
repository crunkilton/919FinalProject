
##Cody Crunkilton
##919 Replication Paper: Horowitz 2009 Replication
##Data available via the Journal of Conflict Resolution website.
##Spring 2017

setwd("C:/Users/Cody/Dropbox/1school2016_7/919/horowitz 2009")

library(readstata13)
library(rjags)
library(runjags)


#############################################getting the data
d<-read.dta13("JCRNuclearAgeReplication.dta")

#############################################removing the missing data

remove <- d[which(d$Sdy=="NA")]
head(remove)

d_noNA <- d[!is.na(d$Sdy),]

dim(d)
dim(d_noNA)

##############################################small model

log_model1<-"
model {for (i in 1:N) {
y[i]~dbern(p[i])
logit(p[i])<-b0+b1*x1[i]+b2*x2[i]+b3*x3[i]+b4*x4[i]
}
b0~dnorm(0,.000001)
b1 ~ dnorm(0,1.0e-4)
b2 ~ dnorm(0,1.0e-4)
b3 ~ dnorm(0,1.0e-4)
b4 ~ dnorm(0,1.0e-4)

}
"
posterior1_year<- run.jags(log_model1_year, monitor=c("b1","b2","b3","b4","b0"), 
         data=with(d,
                   list(y=reciprocation,
                        x1=sideAnuclear,
                        x2=sideBnuclear,
                        x3=sideAnuclearage,
                        x4=sideBnuclearage,
                        N=length(reciprocation))),
         burnin=1000, sample=1000, thin=1,
         n.chains=2)

summary(posterior1)

######################################full model with flat priors

log_modelfull<-"
model {for (i in 1:N) {
y[i]~dbern(p[i])
logit(p[i])<-b0+b1*x1[i]+b2*x2[i]+b3*x3[i]+b4*x4[i]
               +b5*x5[i]+b6*x6[i]+b7*x7[i]+b8*x8[i]
               +b9*x9[i]+b10*x10[i]+b11*x11[i]+b12*x12[i]
}
b0~dnorm(0,1.0e-4)
b1 ~ dnorm(0,1.0e-4)
b2 ~ dnorm(0,1.0e-4)
b3 ~ dnorm(0,1.0e-4)
b4 ~ dnorm(0,1.0e-4)
b5 ~ dnorm(0,1.0e-4)
b6 ~ dnorm(0,1.0e-4)
b7 ~ dnorm(0,1.0e-4)
b8 ~ dnorm(0,1.0e-4)
b9 ~ dnorm(0,1.0e-4)
b10 ~ dnorm(0,1.0e-4)
b11 ~ dnorm(0,1.0e-4)
b12~ dnorm(0,1.0e-4)
}
"

posterior_full<- run.jags(log_modelfull, monitor=c("b1","b2","b3","b4","b0"), 
                      data=with(d_noNA,
                                list(y=reciprocation,
                                     x1=sideAnuclear,
                                     x2=sideBnuclear,
                                     x3=sideAnuclearage,
                                     x4=sideBnuclearage,
                                     x5=jointnuke,
                                     x6=sideAbof,
                                     x7=Sdy,
                                     x8=dem1,
                                     x9=dem2,
                                     x10=territory,
                                     x11=regimegovernment,
                                     x12=policy,
                                     N=length(reciprocation))),
                      burnin=1000, sample=1000, thin=1,
                      n.chains=2)

summary(posterior_full)


#########################################Full model w/ US, Russia controls

log_modelfull_USR<-"
model {for (i in 1:N) {
y[i]~dbern(p[i])
logit(p[i])<-b0+b1*x1[i]+b2*x2[i]+b3*x3[i]+b4*x4[i]
+b5*x5[i]+b6*x6[i]+b7*x7[i]+b8*x8[i]
+b9*x9[i]+b10*x10[i]+b11*x11[i]+b12*x12[i]
+b13*x13[i]+b14*x14[i]
}
b0~dnorm(0,1.0e-4)
b1 ~ dnorm(0,1.0e-4)
b2 ~ dnorm(0,1.0e-4)
b3 ~ dnorm(0,1.0e-4)
b4 ~ dnorm(0,1.0e-4)
b5 ~ dnorm(0,1.0e-4)
b6 ~ dnorm(0,1.0e-4)
b7 ~ dnorm(0,1.0e-4)
b8 ~ dnorm(0,1.0e-4)
b9 ~ dnorm(0,1.0e-4)
b10 ~ dnorm(0,1.0e-4)
b11 ~ dnorm(0,1.0e-4)
b12~ dnorm(0,1.0e-4)
b13~ dnorm(0,1.0e-4)
b14~ dnorm(0,1.0e-4)
}
"

posterior_full_USR<- run.jags(log_modelfull_USR, monitor=c("b1","b2","b3","b4","b0"), 
                          data=with(d_noNA,
                                    list(y=reciprocation,
                                         x1=sideAnuclear,
                                         x2=sideBnuclear,
                                         x3=sideAnuclearage,
                                         x4=sideBnuclearage,
                                         x5=jointnuke,
                                         x6=sideAbof,
                                         x7=Sdy,
                                         x8=dem1,
                                         x9=dem2,
                                         x10=territory,
                                         x11=regimegovernment,
                                         x12=policy,
                                         x13=US,
                                         x14=Russia,
                                         N=length(reciprocation))),
                          burnin=1000, sample=1000, thin=1,
                          n.chains=2)

summary(posterior_full_USR)

###########################################making plots

plot(posterior1, plot.type = "histogram")
plot(posterior_full, plot.type = "histogram")
plot(posterior_full_USR, plot.type = "density") 

#############################################tables for latex
library(stargazer)

model1 <- posterior1$hpd
model2 <- posterior_full$hpd
model3 <- posterior_full_USR$hpd


stargazer(model1)
stargazer(model2)
stargazer(model3)
