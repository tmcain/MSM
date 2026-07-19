
########################################################################################################
## Purpose: Output KM data for PFS, OS, time to progression, time to death and post progression survival 
########################################################################################################

########################
## Required libraries ##
########################

library (survival)
library (flexsurv)
library(haven)
library(rJava)
library(XLConnectJars)
library(XLConnect)
library(survtools)#See 'Packages autoPSM' to install packages
library(dplyr)


##add this in to import data. Also in file to create MSM data set
rm(list=ls())

library(haven)

BresMed.blue = rgb(69, 185, 209, max=255)
BresMed.red = rgb(225, 55, 60, max=255)

## Reading in SAS datasets.
## Program to get around funny folder names:
# directory above

setwd("")
# list me all the files

files <- list.files()
# the first one is the one that we want

files[28]
# use the stored string to set the wd
setwd(paste("./", files[28], sep = ""))
# Now we can go inside the folder without problems
setwd("./Project/Data/21_Jan_2016 cut off/")

files <- list.files()

setwd(paste("./", files[16], sep = ""))

PFSdata <- read_sas(data_file="./pfsfilename.sas7bdat")
OSdata <- read_sas("./osfilename.sas7bdat")

#####






set.seed(12345)
###############
#1. Import data
###############


os.data<-readRDS(file = "Project/Data/OSdata.RDS")%>%
  select(PT, trt, trtf,os, cenos)
#cenpfs=1 is progressed, 0=censored
pfs.data<-readRDS(file = "Project/Data/PFSdata.RDS")%>%
  select(PT, trt, trtf, pfs, cenpfs)
all.data<-merge(os.data, pfs.data)


#################
#1. Read in data
#################

setwd("Project/Data")

#Import survival data
os.data<-readRDS(file = "OSdata.RDS")%>%
  select(PT, trt, trtf,os, cenos)

#Import progression data
#cenpfs=1 is progressed, 0=censored
pfs.data<-readRDS(file = "PFSdata.RDS")
select(PT, trt, trtf, pfs, cenpfs)

###########################################################
#2. Merge data sets and remove observations before 92 weeks
###########################################################

#Combine os and pfs data

all.data<-merge(os.data, pfs.data, by.x='PT', by.y='PT')

#convert to weeks
all.data$os.weeks<-all.data$os*4
all.data$pfs.weeks<-all.data$pfs*4

#delete observations where os and pfs are <92 weeks and only use time after 92 weeks
all.data<-all.data[!(all.data$os.weeks<92 &all.data$pfs.weeks<92),]
all.data$os.weeks<-all.data$os.weeks-92
all.data$pfs.weeks<-all.data$pfs.weeks-92

#covert os and pfs back to months(based on new set of weeks)

all.data$os<-all.data$os.weeks/4
all.data$pfs<-all.data$pfs.weeks/4

##########################################
#3. Create variables (time and event) for 
#    each outcome  TTP, TTD and PPS
#########################################

#

#Create time to progression (TTP) data set (censoring death progression event)

all.data$centtp <- ifelse(all.data$cenpfs == 1 & all.data$cenos == 1 &  (all.data$os == all.data$pfs), 0, all.data$cenpfs)
all.data$ttp<-all.data$pfs

#Create time to death(TTD)(as progression event)

all.data$centtd <- ifelse(all.data$cenpfs == 1 & all.data$cenos == 1 &  (all.data$os == all.data$pfs), 1, 0)
all.data$ttd<-all.data$ppfs

#Create post-progression survival (PPS) data

all.data$pps <- all.data$os - all.data$pfs
all.data$pps <- ifelse(all.data$pps < 0, 0, all.data$pps)

all.data$cenpps<-all.data$cenos


############################
#3. Derive Kaplan-Meier data
############################

# KM data for PFS.#
KMdata.pfs <- survfit(Surv(pfs, cenpfs) ~ trtf, data = all.data,type = "kaplan-meier")

km.curve.pfs <- getKMcurve(km = KMdata.pfs, time.col = 'pfs', event.col = 'cenpfs', data = all.data, group=all.data$trt)

km.curve.A$group <- "KW-0761"
Km.curve.A2<-km.curve.A[c("group","Time","Nrisk","Survival")]

###  medians
KMmedian.A <- summary(KMdata.A)$table["median"]
KM.LCL.A <- summary(KMdata.A)$table["0.95LCL"]
KM.UCL.A <- summary(KMdata.A)$table["0.95UCL"]

KM.summary.A <- data.frame("Treatment"=names(KMmedian.A),KMmedian.A,KM.LCL.A,KM.UCL.A)





km.pfs <- survfit(Surv(pfs, cenpfs) ~ trtc, data = pfs.data,type = "kaplan-meier")

# Km data for OS.
km.os <- survfit(Surv(os, cenos)~trtc, data = os.data,type = "kaplan-meier")

# Km data for ttp.

km.ttp <- survfit(Surv(ttp, centtp)~trtc, data = all.data,type = "kaplan-meier")

# Km data for ttd.

km.ttd <- survfit( Surv(ttd, centtd) ~ trtc, data = all.data, type = "kaplan-meier")

# KM data for pps.

km.pps <- survfit(Surv(pps, cenpps) ~ trtc, data = all.data[(all.data$cenpfs == 1) & (all.data$pps != 0),],type = "kaplan-meier")

