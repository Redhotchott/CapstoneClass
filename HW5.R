###### HW5 #####
#Stat Practicum#
################
rm(list=ls())

#At this point we have 3 models runing around - climatological, baseline and regularized. 
# after this assignment we will have the bss for the last 2. 

setwd('/Volumes/NO NAME/Capstone/HW5')
setwd('E:/CapstoneClass/HW5')
load('predictors.Rdata')
#load('/Volumes/NO NAME/Capstone/predictors.Rdata')
library(sn)
library(fields)
library(mvtnorm)
 
#read.table(file="baseline_classification.txt", header=T)

# dates:        all dates for which data are available
# stations:     names of all stations for which data are available
# lat, lon:     longitude and latitude coordinate for each of these stations
# elev:         station elevation (m)
# ptype:        all cases where one of the four precipitation types of interest was reported
# Twb.prof:     the corresponding vertical profiles of wetbulb temperature (0m to 3000m above the surface in steps of 100m)
# station.ind:  for each case, the index of the station where it was reported
# date.ind:     for each case, the index of the date on which it was reported


#############################
##Answer From Last Assignment
#############################
cols=1:16  	#Columns (i.e. levels) of the temperature profiles to be used in the 
years=as.numeric(substr(dates,1,4))
months=as.numeric(substr(dates,5,6))
all.months=as.numeric(substr(dates[date.ind],5,6))

############################
##Computing Prior Probabilities
############################
prior.probs=array(0,dim=c(length(stations),length(unique(months)),4))

for(i in 1:length(stations)){
  
  print(i)	
  
  for(j in 1:length(unique(months))){
    mon=sort(unique(months))[j]
    #Finding the right stations	
    station.i=which(station.ind==i)
    #Finding the right months
    month.labels=which(months==mon)
    month.rows=which(date.ind%in%month.labels)
    #Getting the right stations AND months
    rows.needed=intersect(station.i,month.rows)
    
    rain.nn=length(which(ptype[rows.needed]=="RA"))
    snow.nn=length(which(ptype[rows.needed]=="SN"))
    pellet.nn=length(which(ptype[rows.needed]=="IP"))
    ice.nn=length(which(ptype[rows.needed]=="FZRA"))
    
    prior.probs[i,j,1:4]=c(rain.nn,snow.nn,pellet.nn,ice.nn)/length(rows.needed)		
  }
}


########################################################
##Find the total number of testing profiles and all of their indices
########################################################
test.nn=array()
ALL.testing.rows=NULL
for(i in 1:12){
  test.years=2000+i
  test.labels=which(years==test.years & months>8 | (years==test.years+1 & months < 6))
  test.rows=which(date.ind%in%test.labels)
  test.nn[i]=length(test.rows)
  ALL.testing.rows=c(ALL.testing.rows,test.rows)}
###################################################################################################################################################################

########################################################
########Baseline (Assignment #4) Classification Model#####
########################################################
prob.hats=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))
colnames(prob.hats)=c("prob.rain","prob.snow","prob.pellets","prob.freezing","observed")
# 
train.nn=array()
test.nn=array()

mean.train=list()
mean.train[[1]]=matrix(0,nrow=16,ncol=12)
mean.train[[2]]=matrix(0,nrow=16,ncol=12)
mean.train[[3]]=matrix(0,nrow=16,ncol=12)
mean.train[[4]]=matrix(0,nrow=16,ncol=12)
names(mean.train)=c("rain","snow","pellets","ice")

cov.train=list()
cov.train[[1]]=list()
cov.train[[2]]=list()
cov.train[[3]]=list()
cov.train[[4]]=list()
names(cov.train)=c("rain","snow","pellets","ice")

ind=0


for(i in 1:12){
  train.years=1996:2000+i-1
  test.years=2000+i

  print(i)

  train.labels=head(which((years>=train.years[1] & months >8)),1):tail(which(years<=train.years[5]+1 & months <6),1)
  test.labels=which((years==test.years & months>8) | (years==test.years+1 & months < 6))


  train.rows=which(date.ind%in%train.labels)
  test.rows=which(date.ind%in%test.labels)

  train.nn[i]=length(train.rows)
  test.nn[i]=length(test.rows)

  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows]=="RA")
  snow.rows=which(ptype[train.rows]=="SN")
  pellet.rows=which(ptype[train.rows]=="IP")
  ice.rows=which(ptype[train.rows]=="FZRA")

  mean.train[[1]][,i]=apply(Twb.prof[train.rows[rain.rows],cols],2,mean)
  mean.train[[2]][,i]=apply(Twb.prof[train.rows[snow.rows],cols],2,mean)
  mean.train[[3]][,i]=apply(Twb.prof[train.rows[pellet.rows],cols],2,mean)
  mean.train[[4]][,i]=apply(Twb.prof[train.rows[ice.rows],cols],2,mean)

  cov.train[[1]][[i]]=cov(Twb.prof[train.rows[rain.rows],cols])
  cov.train[[2]][[i]]=cov(Twb.prof[train.rows[snow.rows],cols])
  cov.train[[3]][[i]]=cov(Twb.prof[train.rows[pellet.rows],cols])
  cov.train[[4]][[i]]=cov(Twb.prof[train.rows[ice.rows],cols])

  #######################################################
  ##Computing probabilities of observations belonging to
  ##each of the 4 groups
  #######################################################

  for(j in 1:test.nn[i]){
    if(j%%1000==0){print(j)}
    ind=ind+1

    station.j=station.ind[test.rows[j]]
    mon.j=months[date.ind[test.rows[j]]]
    mon.col=which(sort(unique(months))==mon.j)

    # #baseline
    # pi.smk=prior.probs[station.j,mon.col,]
    # pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[1]][,i], cov.train[[1]][[i]])
    # pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[2]][,i], cov.train[[2]][[i]])
    # pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[3]][,i], cov.train[[3]][[i]])
    # pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[4]][,i], cov.train[[4]][[i]])
    #
    #climatological
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]
    pi.den.snow=pi.smk[2]
    pi.den.pellet=pi.smk[3]
    pi.den.freeze=pi.smk[4]

    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    prob.hats[ind,1:4]=collection/sum(collection)
    prob.hats[ind,5]=ptype[test.rows[j]]
  }
}
# 
# ##question two just use the prior probs and use the basic prio probs. in divide by n part from hw 4, don't need bc cancels it out 
# #want bss to be closer to 1. 
# # 0 means they are equivalent
# # negative means the reference is better
# #optim maximizes, so we want to take negative of BSS. 
# 
# # #write.table(prob.hats,file="baseline_classification.txt")
# # write.table(prob.hats,file="myclimatological_classification.txt")
# # #prob.hats<-read.table("baseline_classification.txt", header=T)
# # prob.hats2<-read.table("climatology_classification.txt", header=T)
# # prob.hats<-read.table("michaels_classification.txt", header=T)
# # prob.hats<-cbind(prob.hats,prob.hats2[,5])
# 
# ########################################################
# ##Forecast Evaluation---START HERE
# ########################################################
# 
# classes=c("RA","SN","IP","FZRA")
# #classes=c("RA","SN","FZRA", "IP") #for michaels
# BS=0
# for(i in 1:4){
#   
#   matches=which(prob.hats[,5]==classes[i])
#   o.ik=rep(0,length(prob.hats[,5]))
#   o.ik[matches]=rep(1,length(matches))
#   
#   p.ik=prob.hats[,i]
#   
#   BS=BS+sum((p.ik-o.ik)^2,na.rm=T)
#   
# }
# 
# BS
# #BS from HW 4 0.2200711 
# 
# #BS.climatology<- 0.2185553
# #BS.base<-0.2217187
# #BS.NOAA <-
# BS/length(prob.hats[,5])
# 
# #for her's baseline - .2217187
# #climatology -  0.2185553
# #noaa -  0.08749028
# 
# 1-0.2200711/0.2185553
# 
# #Deterministic classification assessment
# observed=prob.hats[,5]
# observed[which(prob.hats[,5]=="RA")]=rep(1,length(which(prob.hats[,5]=="RA")))
# observed[which(prob.hats[,5]=="SN")]=rep(2,length(which(prob.hats[,5]=="SN")))
# observed[which(prob.hats[,5]=="IP")]=rep(3,length(which(prob.hats[,5]=="IP")))
# observed[which(prob.hats[,5]=="FZRA")]=rep(4,length(which(prob.hats[,5]=="FZRA")))
# observed=as.numeric(observed)
# summary(observed)
# table(observed)
# table(prob.hats[,5])
# 
# prob.class=as.matrix(prob.hats[,1:4])
# hard.class=as.integer(apply(prob.class,1,which.max))
# 
# library(s20x)
# 
# CX=crosstabs(hard.class~observed)
# round(CX$whole.props,4)
# sum(diag(CX$whole.props)) #0.8638785, => 86.4% of the profiles are correctly classified
# 


# observed
# class      1      2      		3      	4
# 1 0.6297 0.0614 0.0012 0.0026
# 2 0.0360 0.2280 0.0009 0.0040
# 3 0.0023 0.0021 0.0010 0.0006
# 4 0.0153 0.0093 0.0005 0.0053

################################################################################################################################################################
######################
#######Problem #2#####
######################
BS.func<-function(prob.hats){

  #   B.score<-vector(mode='double',length=nrow(prob.hat))
#   classes<-c('RA','SN', 'FZRA', 'IP')
#   for (i in length(prob.hat)){
#     B.score.i<-vector(mode='double',length=4)
#     for (j in 1:4){
#       if(prob.hat[5,i]==classes[j]){B.score.i[j]<-(as.numeric(prob.hat[j,i])-1)^2} else {B.score.i[j]<-(as.numeric(prob.hat[j,i]))^2}
#     }
#     B.score[i]<-sum(B.score.i)
#   }
#   return(sum(B.score,na.rm=T)/length(B.score))

  classes=c("RA","SN","IP","FZRA")
  #classes=c("RA","SN","FZRA", "IP") #for michaels
  BS=0
  for(i in 1:4){
    matches=which(prob.hats[,5]==classes[i])
    o.ik=rep(0,length(prob.hats[,5]))
    o.ik[matches]=rep(1,length(matches))
    p.ik=prob.hats[,i]
    BS=BS+sum((p.ik-o.ik)^2,na.rm=T)
  }
  return(BS/length(prob.hats[,5]))
  
}


ab.all<-matrix(0,nrow=12, ncol=2)
ab.start<-c(1,10)

cov.set<-function(param,i,set){
  return(param[1]*cov.train[[i]][[set]]+param[2]*diag(1,length(cols)))
}

# in this function m refers to the training set. 
ab.BSS<-function(param,set,BS.ref){
  a=param[1]
  b=param[2]
  if(a <0 | b <0){return(0)}
  
  cov.reg<-lapply(1:4,cov.set, param=c(a,b),set=set)

  prob.hats.reg=data.frame(matrix(0,nrow=sum(train.nn),ncol=5))
  ind=0
  for(j in 1:length(train.rows)){
    ind=ind+1
    
    station.j=station.ind[train.rows[j]]
    mon.j=months[date.ind[train.rows[j]]]
    mon.col=which(sort(unique(months))==mon.j)
    
    #baseline
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[train.rows[j],cols], mean.train[[1]][,set], cov.reg[[1]])
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[train.rows[j],cols], mean.train[[2]][,set], cov.reg[[2]])
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[train.rows[j],cols], mean.train[[3]][,set], cov.reg[[3]])
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[train.rows[j],cols], mean.train[[4]][,set], cov.reg[[4]])
    
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    prob.hats.reg[ind,1:4]=collection/sum(collection)
    prob.hats.reg[ind,5]=ptype[train.rows[j]]
    
  }

  BS=BS.func(prob.hats.reg)
  BSS<-1-BS/BS.ref
  return(-BSS)
}

#initialize properly 
cov.reg<-list()
cov.reg[[1]]<-list()
cov.reg[[2]]<-list()
cov.reg[[3]]<-list()
cov.reg[[4]]<-list()
#length4, each containg 12 
prob.hat.michael=data.frame(matrix(0,nrow=sum(test.nn),ncol=5))

for(i in 1:12){
  train.years=1996:2000+i-1
  test.years=2000+i
    
  print(i)
  
  train.labels=head(which((years>=train.years[1] & months >8)),1):tail(which(years<=train.years[5]+1 & months <6),1)
  test.labels=which((years==test.years & months>8) | (years==test.years+1 & months < 6))
  
  
  train.rows=which(date.ind%in%train.labels)
  test.rows=which(date.ind%in%test.labels)
  
  train.nn[i]=length(train.rows)
  test.nn[i]=length(test.rows)
  
  #######################################################
  ##Computing means and covariances for each precip type
  #######################################################
  rain.rows=which(ptype[train.rows]=="RA")
  snow.rows=which(ptype[train.rows]=="SN")
  pellet.rows=which(ptype[train.rows]=="IP")
  ice.rows=which(ptype[train.rows]=="FZRA")
  
  mean.train[[1]][,i]=apply(Twb.prof[train.rows[rain.rows],cols],2,mean)
  mean.train[[2]][,i]=apply(Twb.prof[train.rows[snow.rows],cols],2,mean)
  mean.train[[3]][,i]=apply(Twb.prof[train.rows[pellet.rows],cols],2,mean)
  mean.train[[4]][,i]=apply(Twb.prof[train.rows[ice.rows],cols],2,mean)
  
  cov.train[[1]][[i]]=cov(Twb.prof[train.rows[rain.rows],cols])
  cov.train[[2]][[i]]=cov(Twb.prof[train.rows[snow.rows],cols])
  cov.train[[3]][[i]]=cov(Twb.prof[train.rows[pellet.rows],cols])
  cov.train[[4]][[i]]=cov(Twb.prof[train.rows[ice.rows],cols])
  
  
  #######################################################
  ##Computing the climatological BS score over the training set
  #######################################################
  prob.hats.clim=data.frame(matrix(0,nrow=train.nn[i],ncol=5))
  ind=0
  for(j in 1:length(train.rows)){
    if(j%%1000==0){print(paste('training set: ',i, 'observation: ',j))}
    ind=ind+1
    
    station.j=station.ind[train.rows[j]]
    mon.j=months[date.ind[train.rows[j]]]
    mon.col=which(sort(unique(months))==mon.j)
  
    #climatological
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]
    pi.den.snow=pi.smk[2]
    pi.den.pellet=pi.smk[3]
    pi.den.freeze=pi.smk[4]
    
    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    prob.hats.clim[ind,1:4]=collection/sum(collection)
    prob.hats.clim[ind,5]=ptype[train.rows[j]]
  }
  
  BS.ref<-BS.func(prob.hats.clim) #for first training set 0.2400417
  
  #######################################################
  ##computing the optimal a and b for use in the regularization
  #######################################################

  ab.start=optim(ab.start, ab.BSS, set=i, BS.ref=BS.ref, control=list(parscale=c(.1,3)))$par
  ab.all[i,]=ab.start 

  
  cov.reg[[1]][[i]]=ab.start[1]*cov.train[[1]][[i]]+ab.start[2]*diag(1,length(cols))
  cov.reg[[2]][[i]]=ab.start[1]*cov.train[[2]][[i]]+ab.start[2]*diag(1,length(cols))
  cov.reg[[3]][[i]]=ab.start[1]*cov.train[[3]][[i]]+ab.start[2]*diag(1,length(cols))
  cov.reg[[4]][[i]]=ab.start[1]*cov.train[[4]][[i]]+ab.start[2]*diag(1,length(cols))
  
  
  
  #######################################################
  ##Reclassify using the new sigmas 
  #######################################################
  ind=0
  for(j in 1:length(test.rows)){
    if(j%%1000==0){print(j)}
    ind=ind+1

    station.j=station.ind[test.rows[j]]
    mon.j=months[date.ind[test.rows[j]]]
    mon.col=which(sort(unique(months))==mon.j)

    #baseline
    pi.smk=prior.probs[station.j,mon.col,]
    pi.den.rain=pi.smk[1]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[1]][,i], cov.reg[[1]][[i]])
    pi.den.snow=pi.smk[2]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[2]][,i], cov.reg[[2]][[i]])
    pi.den.pellet=pi.smk[3]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[3]][,i], cov.reg[[3]][[i]])
    pi.den.freeze=pi.smk[4]*dmvnorm(Twb.prof[test.rows[j],cols], mean.train[[4]][,i], cov.reg[[4]][[i]])

    collection=c(pi.den.rain,pi.den.snow,pi.den.pellet,pi.den.freeze)
    prob.hat.michael[ind,1:4]=collection/sum(collection)
    prob.hat.michael[ind,5]=ptype[test.rows[j]]
  }
  
}




write.table(ab.all,file="aball.txt")
write.table(prob.hat.michael,file="mymichaeltest.txt")






###################
##Problem #2
###################
# reliability<-function(for.probs){
#   p.type<-c("RA", "SN", "IP", "FZRA")
#   rain.rows<-which(for.probs[,5]=="RA")
#   snow.rows<-which(for.probs[,5]=="SN")
#   pellets<-which(for.probs[,5]=="IP")
#   fzrain.rows<-which(for.probs[,5]=="FZRA")
#   title("Reliability of preciptiation type probabilities")
# }
# 
# ###################
# ##Problem #4
# ###################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###################
# ##Problem #2 FROM HER ORIGINAL STUFF. 
# ###################
# 
# 
# 
# 
# years=as.numeric(substr(dates,1,4))
# months=as.numeric(substr(dates,5,6))
# 
# train.nn=array()
# test.nn=array()
# 
# mean.train=list()
# mean.train[[1]]=matrix(0,nrow=16,ncol=12)
# mean.train[[2]]=matrix(0,nrow=16,ncol=12)
# mean.train[[3]]=matrix(0,nrow=16,ncol=12)
# mean.train[[4]]=matrix(0,nrow=16,ncol=12)
# names(mean.train)=c("rain","snow","pellets","ice")
# 
# cov.train=list()
# cov.train[[1]]=list()
# cov.train[[2]]=list()
# cov.train[[3]]=list()
# cov.train[[4]]=list()
# 
# names(cov.train)=c("rain","snow","pellets","ice")
# 
# for(i in 1:12){
#   train.years=1996:2000+i-1
#   test.years=2000+i
#   
#   print(i)
#   
#   train.labels=head(which((years>=train.years[1] & months >8)),1):tail(which(years<=train.years[5]+1 & months <6),1)
#   test.labels=which((years==test.years & months>8) | (years==test.years+1 & months < 6))
#   
#   
#   train.rows=which(date.ind%in%train.labels)
#   test.rows=which(date.ind%in%test.labels)
#   
#   train.nn[i]=length(train.rows)
#   test.nn[i]=length(test.rows)
#   
#   #######################################################
#   ##Computing means and covariances for each precip type
#   #######################################################
#  
#   
# }


