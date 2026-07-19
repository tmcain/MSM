

library(haven)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

###############
#1. Import data
###############

#ITT data-a data set containing both PFS and OS data
deff.data <-... 


########################
#2. Manipulate data sets
########################

#Create subset of ITT data, convert PFS and OS time from months to days and check no PFS times >OS times

itt<-deff.data%>%
            select(study, patid, trt,trtc,pfs, cenpfs, os, cenos, randdt,pfsdt, osdt )
itt$pfs.days<-round((itt$pfs/12)*365.25)
itt$os.days<-round((itt$os/12)*365.25)
itt[itt$pfs.days>itt$os.days,]$pfs.days<-itt[itt$pfs.days>itt$os.days,]$os.days


#####################################################
# 3.Function to create MSM states for each individual
#####################################################

MSM.state<-function(data, subjectID, treatmentID, ostime, pfstime, oscens, pfscens )
{
  Patient <- vector()
  Treatment <- vector()
  State<-vector()
  Time<-vector()

  for (i in 1:nrow(data))
  {
    Patient<-c(Patient,subjectID[i])
    Treatment<-c(Treatment, treatmentID[i])
    Time<-c(Time,0)
    State<-c(State,1)
    
    #If PFS and OS times have both been censored
    if (pfscens[i]==0&oscens[i]==0)
    {
      if (pfstime[i]==ostime[i])
      {
        Patient<-c(Patient,subjectID[i])
        Treatment<-c(Treatment, treatmentID[i])
        Time<-c(Time,ostime[i])
        State<-c(State,1)
      }
      else #Times not equal
      {
        Patient<-c(Patient,rep(subjectID[i],2))
        Treatment<-c(Treatment, rep(treatmentID[i],2))
        Time<-c(Time,pfstime[i]+1, ostime[i])
        State<-c(State,99,99)   
      }
    }
    
    #If PFS and OS times have both been observed
    if (pfscens[i]==1&oscens[i]==1) 
    {
      if (pfstime[i]==ostime[i])
      {
        Patient<-c(Patient,subjectID[i])
        Treatment<-c(Treatment, treatmentID[i])
        Time<-c(Time,ostime[i])
        State<-c(State,3)
      } 
      else #Times not equal
      {
        Patient<-c(Patient,rep(subjectID[i],2))
        Treatment<-c(Treatment, rep(treatmentID[i],2))
        Time<-c(Time,pfstime[i], ostime[i])
        State<-c(State,2,3)   
      }
    } 
    
    #If PFS is observed and OS is censored
    if (pfscens[i]==1&oscens[i]==0) 
    {
      Patient<-c(Patient,rep(subjectID[i],2))
      Treatment<-c(Treatment, rep(treatmentID[i],2))
      Time<-c(Time,pfstime[i], ostime[i])
      State<-c(State,2,2)   
    }
    
    #If PFS is censored and OS is observed
    if (pfscens[i]==0&oscens[i]==1) 
    {
      Patient<-c(Patient,rep(subjectID[i],2))
      Treatment<-c(Treatment, rep(treatmentID[i],2))
      Time<-c(Time,pfstime[i]+1, ostime[i])
      State<-c(State,99,3)  
    }
  }
  data.frame(Patient, Treatment,State, Time)
}
################################################################################
#4. Create MSM data sets for each of the three data sets using the function in 3
################################################################################

#ITT
states.itt<-MSM.state(itt, itt$patid, itt$trtc, itt$os.days,itt$pfs.days, itt$cenos, itt$cenpfs )


#####################################################################################################
#5. Fit MSM models to the three data sets and output probability matrix using PFS and OS time in days
#####################################################################################################

library(msm)
qmatrix<-matrix(c(1,0,0,1,1,0,1,1,0),3,3)

#ITT
msm.itt.arm1<-msm(State~Time,subject=Patient, data=states.itt[states.itt$Treatment=="Arm1",], qmatrix = qmatrix,gen.inits = TRUE ,deathexact = 3, obstype=2,  censor = 99, censor.states = list(c(1,2))) 
msm.itt.arm2<-msm(State~Time,subject=Patient, data=states.itt[states.itt$Treatment=="Arm2",], qmatrix = qmatrix,gen.inits = TRUE ,deathexact = 3, obstype=2,  censor = 99, censor.states = list(c(1,2))) 

prob.itt.arm1<-pmatrix.msm(msm.itt.arm1)
prob.itt.arm2<-pmatrix.msm(msm.itt.arm2)


#####################################################################################
#6. Calculate probabilities of each state over time(10000 days) and convert to months
#####################################################################################

#ITT
states.days<-seq(1,10000,1)
states.months<-(states.days/365.25)*12
start.numbers <- c(460,0,0)

states.itt.arm1 <- matrix(NA,10000,3)

states.itt.arm1[1,]<-start.numbers

for(i in 2:dim(states.itt.arm1)[1])
{
  states.itt.arm1[i,] <-states.itt.arm1[i-1,] %*% prob.itt.arm1
} 

states.itt.arm2 <- matrix(NA,10000,3)

states.itt.arm2[1,]<-start.numbers

for(i in 2:dim(states.itt.arm2)[1])
{
  states.itt.arm2[i,] <-   states.itt.arm2[i-1,] %*% prob.itt.arm2
}  




ITT.OS.arm2.days  <- data.frame(time = states.days,est_OS = 1-(states.itt.arm2[,3]/start.numbers[1]))
colnames(ITT.OS.arm2.days)<-c("Time", "OS")

ITT.PFS.arm2.days  <- data.frame(time = states.days,est_PFS = 1-( (states.itt.arm2[,2]+states.itt.arm2[,3])/start.numbers[1]))
colnames(ITT.PFS.arm2.days)<-c("Time", "PFS")

ITT.OS.arm2.months  <- data.frame(time = states.months,est_OS = 1-(states.itt.arm2[,3]/start.numbers[1]))
colnames(ITT.OS.arm2.months)<-c("Time", "OS")

ITT.PFS.arm2.months  <- data.frame(time = states.months,est_PFS = 1-( (states.itt.arm2[,2]+states.itt.arm2[,3])/start.numbers[1]))
colnames(ITT.PFS.arm2.months)<-c("Time", "PFS")

ITT.OS.arm1.days  <- data.frame(time = states.days,est_OS = 1-(states.itt.arm1[,3]/start.numbers[1]))
colnames(ITT.OS.arm1.days)<-c("Time", "OS")

ITT.PFS.arm1.days  <- data.frame(time = states.days,est_PFS = 1-( (states.itt.arm1[,2]+states.itt.arm1[,3])/start.numbers[1]))
colnames(ITT.PFS.arm1.days)<-c("Time", "PFS")

ITT.OS.arm1.months  <- data.frame(time = states.months,est_OS = 1-(states.itt.arm1[,3]/start.numbers[1]))
colnames(ITT.OS.arm1.months)<-c("Time", "OS")

ITT.PFS.arm1.months  <- data.frame(time = states.months,est_PFS = 1-( (states.itt.arm1[,2]+states.itt.arm1[,3])/start.numbers[1]))
colnames(ITT.PFS.arm1.months)<-c("Time", "PFS")



###################################################
#7. Output probabilities over months as r data sets
###################################################

setwd('')

saveRDS(ITT.OS.arm2.months,file="ITT.OS.arm2.Rda")
saveRDS(ITT.PFS.arm2.months,file="ITT.PFS.arm2.Rda")
saveRDS(ITT.OS.arm1.months,file="ITT.OS.arm1.Rda")
saveRDS(ITT.PFS.arm1.months,file="ITT.PFS.arm1.Rda")


######################################################
#8. Output csv files of msm state probability matrices
######################################################
setwd('')

write.csv(prob.itt.arm1,"Matrix.ITT.arm1.csv") 
write.csv(prob.itt.arm2, "Matrix.ITT.arm2.csv")


##############################################################
#9. Create and save bootstrap smaples of each transition matrix
##############################################################

setwd('')

# ITT arm2
msmbootlog.arm2 <- file("./arm1.log")
sink(msmbootlog.arm2, append=TRUE)
sink(msmbootlog.arm2, append=TRUE,type="message")
boot.pmatrix.ITT.arm2 <- boot.msm(msm.itt.arm2, B=50,stat=prob.itt.arm2)
boot.pmatrix.ITT.arm2
sink() 


# ITT arm1
msmbootlog.itt.arm1 <- file("./arm2.log")
sink(msmbootlog.itt.arm1, append=TRUE)
sink(msmbootlog.itt.arm1, append=TRUE,type="message")
boot.pmatrix.ITT.arm1 <- boot.msm(msm.itt.arm1, B=50,stat=prob.itt.arm1)
boot.pmatrix.ITT.arm1
sink() 


