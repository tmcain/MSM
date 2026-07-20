
########################################################
#Author: Theresa Cain              #
# Aim: To re-create the Kaplan-Meier Curves with      #
# superimposed MSM models.                            #
# in line with R convention: 1 = event, O = censor. # 
######################################################## 

rm(list=ls(all=TRUE))

library(survival)
library(survminer)
library(ggplot2)
library(haven)
library(dplyr)

blue = rgb(69, 185, 209, max=255)
red = rgb(225, 55, 60, max=255)
yellow = rgb(238, 224, 30, max=255)
pink=rgb(211,78,147,max=255)
Dblue=rgb(0,45,92,max=255)
Dyellow = rgb(214, 200, 16, max=255)

###############
#1. Import data
###############
#ITT data

setwd('')
deff.data <- read_sas(data_file = 'deff.sas7bdat')

#RPSFT/censor at crossover data
setwd('')
PFS.data <- read.csv( " Progression Free Survival -  Treatment Group method .csv" , header = TRUE) 
OS.data<-read.csv(" Overall Survival - Treatment Group method .csv" , header = TRUE)

## Loading in the MSM data for PFS and OS.

setwd('')

ITT.OS.treatmentA <- readRDS("ITT.OS.treatmentA.Rda")
ITT.OS.Placebo <- readRDS("ITT.OS.Placebo.Rda")

ITT.PFS.treatmentA <- readRDS("ITT.PFS.treatmentA.Rda")
ITT.PFS.Placebo <- readRDS("ITT.PFS.Placebo.Rda")

## Note, you do not need to repeat the treatmentA arm here. You use the same treatmentA.

RPSFT.OS.Placebo <- readRDS("RPSFT.OS.Placebo.Rda")
RPSFT.PFS.Placebo <- readRDS("RPSFT.PFS.Placebo.Rda")

crossover.OS.Placebo <- readRDS("crossover.OS.Placebo.Rda")
crossover.PFS.Placebo <- readRDS("crossover.PFS.Placebo.Rda")

## Loading in the MSM data for Individual Transitions

setwd('MSM Output/Individual Transitions')

## Transition 1-2 ##

transition.12.itt.trtA <- readRDS("transition.12.itt.trtA.Rda")
transition.12.itt.placebo <- readRDS("transition.12.itt.placebo.Rda")

transition.12.rpsft.placebo <- readRDS("transition.12.rpsft.placebo.Rda")
transition.12.cen.at.cross.placebo <- readRDS("transition.12.cen.at.cross.placebo.Rda")

## Transition 1-3 ##

transition.13.itt.trtA <- readRDS("transition.13.itt.trtA.Rda")
transition.13.itt.placebo <- readRDS("transition.13.itt.placebo.Rda")

transition.13.rpsft.placebo <- readRDS("transition.13.rpsft.placebo.Rda")
transition.13.cen.at.cross.placebo <- readRDS("transition.13.cen.at.cross.placebo.Rda")

## Transition 2-3 ##

transition.23.itt.trtA <- readRDS("transition.23.itt.trtA.Rda")
transition.23.itt.placebo <- readRDS("transition.23.itt.placebo.Rda")

transition.23.rpsft.placebo <- readRDS("transition.23.rpsft.placebo.Rda")
transition.23.cen.at.cross.placebo <- readRDS("transition.23.cen.at.cross.placebo.Rda")

#############################################
# 2. Derive Individual Transition Endpoints
#############################################
# 2.a. Time to progression.
#############################################

## TTP: in line with R convention: 1 = event, O = censor. #

########################################################################################################################
########################################################################################################################

deff.data$centtp <- ifelse(deff.data$cenpfs == 1 & deff.data$cenos == 1 &  (deff.data$os == deff.data$pfs), 0, deff.data$cenpfs)
death.as.pfs.event <- c(deff.data[deff.data$cenpfs == 1 & deff.data$cenos == 1 &  (deff.data$os == deff.data$pfs), c("patid")]) 
# Identifies patients who have death as PFS event.

# Censor at progression: in line with R convention: 1 = event, O = censor. #

# Patients with death as pfs event get censored. All other censor flags remain the same.
PFS.data$cen_cen_at_cross_ttp <- ifelse(PFS.data$patid_c %in% death.as.pfs.event$patid,0,PFS.data$cen_cen_at_cross)

# Patients with death as pfs event get censored. All other censor flags remain the same.
PFS.data$Untreated.Counterfactual.cens_ttp <- ifelse(PFS.data$patid_c %in% death.as.pfs.event$patid,0,PFS.data$Untreated.Counterfactual.Cens)

#############################################
# 2.b. Time to Death (ttd).                 #
#############################################

## Time to Death (ttd): 1 = event, O = censor. ##

deff.data$centtd <- ifelse(deff.data$cenpfs == 1 & deff.data$cenos == 1 &  (deff.data$os == deff.data$pfs), 1, 0)

# Censored at Crossover: Patients not censored at crossover ( PFS.data$event == PFS.data$cen_cen_at_cross & round(PFS.data$tte,3) == round(PFS.data$cen_at_cross,3) ) 
# who have death as pfs event (OS.data$patid_c %in% death.as.pfs.event$patid ). Get events: all others censored.

PFS.data$cen_cen_at_cross_ttd <- ifelse( (PFS.data$event == PFS.data$cen_cen_at_cross ) & 
                                           round(PFS.data$tte,3) == round(PFS.data$cen_at_cross,3) & 
                                           (PFS.data$patid_c %in% death.as.pfs.event$patid ), 1, 0)

# RPSFTM:  Patients who have death as pfs event (OS.data$patid_c %in% death.as.pfs.event$patid ) get events: all others censored.
PFS.data$Untreated.Counterfactual.cens_ttd <- ifelse(PFS.data$patid_c %in% death.as.pfs.event$patid, 1, 0)

## Post-prOgression survival: 1 = event, O = censor.##

#############################################
# 2.c. Post-progression Survival
#############################################

## PPS: ITT  ##

deff.data$pps <- deff.data$os - deff.data$pfs

## There is one patient ("113219") who has a psf one day after OS.

deff.data$pps <- ifelse(deff.data$pps < 0, 0, deff.data$pps)

## PPS: Censored at crossover  ##

## Deriving a dataset in order to derive post-progression survival for OS and PFS for censored at crossover.  

## OS ##

OS.data$OS_cen_at_cross <- OS.data$cen_at_cross
OS.data$OS_cen_cen_at_cross <- OS.data$cen_cen_at_cross

OS.data$OS_Untreated.Counterfactual.Times  <- OS.data$Untreated.Counterfactual.Times  
OS.data$OS_Untreated.Counterfactual.Cens  <- OS.data$Untreated.Counterfactual.Cens

# PFS # 
PFS.data$PFS_cen_at_cross <- PFS.data$cen_at_cross
PFS.data$PFS_cen_cen_at_cross <- PFS.data$cen_cen_at_cross

PFS.data$PFS_Untreated.Counterfactual.Times  <- PFS.data$Untreated.Counterfactual.Times  
PFS.data$PFS_Untreated.Counterfactual.Cens  <- PFS.data$Untreated.Counterfactual.Cens
  
################################################################################################

PFS_OS.data <- merge(OS.data[,c("patid_c","OS_cen_at_cross","OS_cen_cen_at_cross","OS_Untreated.Counterfactual.Times","OS_Untreated.Counterfactual.Cens","trtc")],
      PFS.data[,c("patid_c","PFS_cen_at_cross","PFS_cen_cen_at_cross","PFS_Untreated.Counterfactual.Times","PFS_Untreated.Counterfactual.Cens")],by= "patid_c")

################################################################################################

PFS_OS.data$PPS_cen_at_cross <-  PFS_OS.data$OS_cen_at_cross - PFS_OS.data$PFS_cen_at_cross

# There are a number of zero times; This appears to be some type of rounding issue from Excel. These will be rounded up.
# Those with a difference of less than a day are rounded to zero.

PFS_OS.data$PPS_cen_at_cross <- ifelse(PFS_OS.data$PPS_cen_at_cross < 0.00001,0, PFS_OS.data$PPS_cen_at_cross)  

## PPS: RPSFTM ##

PFS_OS.data$PPS_Untreated.Counterfactual.Times <-  PFS_OS.data$OS_Untreated.Counterfactual.Times - PFS_OS.data$PFS_Untreated.Counterfactual.Times

## There is one patient ("113219") who has a PFS one day after OS.

PFS_OS.data$PPS_Untreated.Counterfactual.Times <- ifelse(PFS_OS.data$PPS_Untreated.Counterfactual.Times < 0,0,PFS_OS.data$PPS_Untreated.Counterfactual.Times)

#############################
# 3.Derive Kaplan-Meier data
#############################
# KM data for PFS.#

km.pfs.itt <- survfit(Surv(pfs, cenpfs) ~ trtc, data = deff.data,type = "kaplan-meier")
km.pfs.rpsft <- survfit(Surv(Untreated.Counterfactual.Times, Untreated.Counterfactual.Cens)~trtc, data = PFS.data, type = "kaplan-meier")
km.pfs.censoring.at.crossover <- survfit(Surv(cen_at_cross, cen_cen_at_cross)~trtc, data = PFS.data, type = "kaplan-meier")

# Km data for OS.
km.os.itt <- survfit(Surv(os, cenos)~trtc, data = deff.data,type = "kaplan-meier")
km.os.rpsft <- survfit(Surv(Untreated.Counterfactual.Times, Untreated.Counterfactual.Cens)~trtc, data = OS.data, type = "kaplan-meier")
km.OS.censoring.at.crossover <- survfit(Surv(cen_at_cross, cen_cen_at_cross)~trtc, data = OS.data, type = "kaplan-meier")

# Km data for ttp.

km.ttp.itt <- survfit(Surv(pfs, centtp)~trtc, data = deff.data,type = "kaplan-meier")
km.ttp.rpsft <- survfit(Surv(Untreated.Counterfactual.Times, Untreated.Counterfactual.cens_ttp) ~ trtc, data = PFS.data, type = "kaplan-meier")
km.ttp.censoring.at.crossover <- survfit(Surv(cen_at_cross, cen_cen_at_cross_ttp) ~ trtc, data = PFS.data, type = "kaplan-meier")

# Km data for ttd.

km.ttd.itt <- survfit( Surv(pfs, centtd) ~ trtc, data = deff.data, type = "kaplan-meier")
km.ttd.rpsft <- survfit(Surv(Untreated.Counterfactual.Times, Untreated.Counterfactual.cens_ttd) ~ trtc, data = PFS.data, type = "kaplan-meier")
km.ttd.censoring.at.crossover <- survfit(Surv(cen_at_cross, cen_cen_at_cross_ttd) ~ trtc, data = PFS.data, type = "kaplan-meier")

# KM data for pps.

km.pps.itt <- survfit(Surv(pps, cenos) ~ trtc, data = deff.data[(deff.data$cenpfs == 1) & (deff.data$pps != 0),],type = "kaplan-meier")
km.pps.rpsft <- survfit(Surv(PPS_Untreated.Counterfactual.Times,OS_Untreated.Counterfactual.Cens) ~ trtc, 
                  data = PFS_OS.data[ (PFS_OS.data$PFS_Untreated.Counterfactual.Cens == 1) & (PFS_OS.data$PPS_Untreated.Counterfactual.Times != 0), ],type = "kaplan-meier")  
km.pps.cen_at_cross <- survfit(Surv(PPS_cen_at_cross, OS_cen_cen_at_cross) ~ trtc, 
      data = PFS_OS.data[(PFS_OS.data$PPS_cen_at_cross > 0) & (PFS_OS.data$PFS_cen_cen_at_cross == 1),], type = "kaplan-meier")



#############################################################################
# 4.Plot PFS Kaplan-Meier curves of treatmentA vs. each placebo separately
#############################################################################

surv.summary.km.pfs.itt <- surv_summary(km.pfs.itt,data = deff.data)

surv.summary.km.pfs.itt.Len <- surv.summary.km.pfs.itt[surv.summary.km.pfs.itt$trtc == "treatmentA",]
surv.summary.km.pfs.itt.Plac <- surv.summary.km.pfs.itt[surv.summary.km.pfs.itt$trtc == "Placebo",]

#PFS KM curve for treatmentA and Placebo-unadjusted(ITT)
PFS.itt <- ggsurvplot(km.pfs.itt,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.pfs.itt.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.pfs.itt.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = ITT.PFS.treatmentA, aes(x = Time, y = PFS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = ITT.PFS.Placebo, aes(x = Time, y = PFS, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Progression free survival (ITT) \n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" = 1,  "Placebo-unadjusted-KM" = 1, "treatmentA-MSM" = 3, "Placebo-unadjusted-MSM" = 3 )) +
  guides(colour = guide_legend(nrow = 2))

PFS.itt.at.risk <- ggsurvplot(km.pfs.itt, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                              palette=c(blue,red),xlab = "Time (months)",break.time.by = 12,risk.table = TRUE,tables.y.text = FALSE)

# PFS KM curve for treatmentA and Placebo-adjusted RPSFTM recensoring

surv.summary.km.pfs.rpsft <- surv_summary(km.pfs.rpsft,data = PFS.data)

surv.summary.km.pfs.rpsft.Len <- surv.summary.km.pfs.rpsft[surv.summary.km.pfs.rpsft$trtc == "treatmentA",]
surv.summary.km.pfs.rpsft.Plac <- surv.summary.km.pfs.rpsft[surv.summary.km.pfs.rpsft$trtc == "Placebo",]

PFS.rpsft <- ggsurvplot(km.pfs.rpsft,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.pfs.rpsft.Len , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.pfs.rpsft.Plac, mapping = aes(x = time, y = surv,colour="Placebo-RPSFT-KM",linetype="Placebo-RPSFT-KM" ), size = 0.75) +
  geom_line(data = ITT.PFS.treatmentA, aes(x = Time, y = PFS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = RPSFT.PFS.Placebo, aes(x = Time, y = PFS, colour = "Placebo-RPSFT-MSM", linetype = "Placebo-RPSFT-MSM"),size = 0.75) +
  labs(title = " Progression free survival for RPSFT adjustment method \n - 'treatment group method'\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-RPSFT-KM" = red ,"treatmentA-MSM"= blue, "Placebo-RPSFT-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1,  "Placebo-RPSFT-KM"=1, "treatmentA-MSM"=3, "Placebo-RPSFT-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

PFS.rpsft.at.risk <- ggsurvplot(km.pfs.rpsft, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                                palette=c(blue,red),xlab = "Time (months)",break.time.by = 12,risk.table = TRUE,tables.y.text = FALSE)

# PFS KM curve for treatmentA and Placebo-censoring at XO

surv.summary.km.censoring.at.crossover <- surv_summary(km.pfs.censoring.at.crossover,data = PFS.data)

surv.summary.km.censoring.at.crossover.Len  <- surv.summary.km.censoring.at.crossover[surv.summary.km.censoring.at.crossover$trtc == "treatmentA",]
surv.summary.km.censoring.at.crossover.Plac <- surv.summary.km.censoring.at.crossover[surv.summary.km.censoring.at.crossover$trtc == "Placebo",]

PFS.crossover <- ggsurvplot(km.pfs.censoring.at.crossover,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.censoring.at.crossover.Len   , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.censoring.at.crossover.Plac , mapping = aes(x = time, y = surv,colour="Placebo-Censoring-KM",linetype="Placebo-Censoring-KM" ), size = 0.75) +
  geom_line(data = ITT.PFS.treatmentA, aes(x = Time, y = PFS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = crossover.PFS.Placebo, aes(x = Time, y = PFS, colour = "Placebo-Censoring-MSM", linetype = "Placebo-Censoring-MSM"),size = 0.75) +
  labs(title =" Progression free survival \n Censored at crossover\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-Censoring-KM" = red ,"treatmentA-MSM"= blue, "Placebo-Censoring-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1,  "Placebo-Censoring-KM"=1, "treatmentA-MSM"=3, "Placebo-Censoring-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))

PFS.crossover.at.risk <- ggsurvplot(km.pfs.censoring.at.crossover, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                                    palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

##############################################################################
# 5.Plot OS Kaplan-Meier curves of treatmentA vs. each placebo separately.
##############################################################################

# OS KM curve for treatmentA and Placebo-unadjusted
surv.summary.km.os.itt <- surv_summary(km.os.itt,data = deff.data)

surv.summary.km.os.itt.Len <- surv.summary.km.os.itt[surv.summary.km.os.itt$trtc == "treatmentA",]
surv.summary.km.os.itt.Plac <- surv.summary.km.os.itt[surv.summary.km.os.itt$trtc == "Placebo",]

OS.itt <- ggsurvplot(km.os.itt,break.time.by = 24, colour=NA, censor = FALSE, xlim = c(0, 300), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.os.itt.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.os.itt.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = ITT.OS.treatmentA, aes(x = Time, y = OS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = ITT.OS.Placebo, aes(x = Time, y = OS, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Overall survival (ITT) \n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

OS.itt.at.risk <- ggsurvplot(km.os.itt, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 300),
                             palette =c(blue,red) ,xlab = "Time (months)", break.time.by = 24, risk.table = TRUE,tables.y.text = FALSE)

# OS KM curve for treatmentA and Placebo-adjusted RPSFTM recensoring

surv.summary.km.os.rpsft <- surv_summary(km.os.rpsft,data = OS.data)

surv.summary.km.os.rpsft.Len <- surv.summary.km.os.rpsft[surv.summary.km.os.rpsft$trtc == "treatmentA",]
surv.summary.km.os.rpsft.Plac <- surv.summary.km.os.rpsft[surv.summary.km.os.rpsft$trtc == "Placebo",]

OS.rpsft <- ggsurvplot(km.os.rpsft,break.time.by = 24, colour=NA, censor = FALSE, xlim = c(0, 300), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.os.rpsft.Len , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.os.rpsft.Plac, mapping = aes(x = time, y = surv,colour="Placebo-RPSFTM-KM",linetype="Placebo-RPSFTM-KM" ), size = 0.75) +
  geom_line(data = ITT.OS.treatmentA, aes(x = Time, y = OS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = RPSFT.OS.Placebo, aes(x = Time, y = OS, colour = "Placebo-RPSFTM-MSM", linetype = "Placebo-RPSFTM-MSM"),size = 0.75) +
  labs(title = " Overall survival for RPSFT adjustment method \n - 'treatment group method'\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-RPSFTM-KM" = red ,"treatmentA-MSM"= blue, "Placebo-RPSFTM-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM"= 3, "Placebo-RPSFTM-KM"=1, "Placebo-RPSFTM-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

OS.rpsft.at.risk <- ggsurvplot(km.os.rpsft, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 300),
                               palette =c(blue,red) ,xlab = "Time (months)", break.time.by = 24, risk.table = TRUE,tables.y.text = FALSE)

# OS KM curve for treatmentA and Placebo-censoring at XO

surv.summary.km.OS.censoring.at.crossover <- surv_summary(km.OS.censoring.at.crossover,data = OS.data)

surv.summary.km.OS.censoring.at.crossover.Len <- surv.summary.km.OS.censoring.at.crossover[surv.summary.km.OS.censoring.at.crossover$trtc == "treatmentA",]
surv.summary.km.OS.censoring.at.crossover.Plac <- surv.summary.km.OS.censoring.at.crossover[surv.summary.km.OS.censoring.at.crossover$trtc == "Placebo",]

OS.censoring.at.crossover <- ggsurvplot(km.OS.censoring.at.crossover,break.time.by = 24, colour=NA, censor = FALSE, xlim = c(0, 300), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.OS.censoring.at.crossover.Len , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.OS.censoring.at.crossover.Plac , mapping = aes(x = time, y = surv,colour="Placebo-Censoring-KM",linetype="Placebo-Censoring-KM" ), size = 0.75) +
  geom_line(data = ITT.OS.treatmentA, aes(x = Time, y = OS,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = crossover.OS.Placebo, aes(x = Time, y = OS, colour = "Placebo-Censoring-MSM", linetype = "Placebo-Censoring-MSM"),size = 0.75) +
  labs(title = " Overall survival \n Censored at crossover\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-Censoring-KM" = red ,"treatmentA-MSM"= blue, "Placebo-Censoring-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM"= 3, "Placebo-Censoring-KM"=1, "Placebo-Censoring-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

OS.censoring.at.crossover.at.risk <- ggsurvplot(km.OS.censoring.at.crossover, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 300),
                                                palette =c(blue,red) ,xlab = "Time (months)", break.time.by = 24, risk.table = TRUE,tables.y.text = FALSE)

#############################################################################
# 6.Plot TTP Kaplan-Meier curves of treatmentA vs. each placebo separately
#############################################################################

# PFS KM curve for treatmentA and Placebo-unadjusted(ITT)

surv.summary.km.ttp.itt <- surv_summary(km.ttp.itt,data = deff.data)

surv.summary.km.ttp.itt.Len <- surv.summary.km.ttp.itt[surv.summary.km.ttp.itt$trtc == "treatmentA", ]
surv.summary.km.ttp.itt.Plac <- surv.summary.km.ttp.itt[surv.summary.km.ttp.itt$trtc == "Placebo", ]

ttp.itt <- ggsurvplot(km.ttp.itt,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.ttp.itt.Len   , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.ttp.itt.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = transition.12.itt.trtA , aes(x = Time, y = P12,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.12.itt.placebo , aes(x = Time, y = P12, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Time to progression (ITT) \n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))

ttp.itt.at.risk <- ggsurvplot(km.ttp.itt, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                              palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

# PFS KM curve for treatmentA and Placebo-unadjusted (ITT)

surv.summary.km.ttp.rpsft <- surv_summary(km.ttp.rpsft,data = deff.data)

surv.summary.km.ttp.rpsft.Len <- surv.summary.km.ttp.rpsft[surv.summary.km.ttp.rpsft$trtc == "treatmentA",]
surv.summary.km.ttp.rpsft.Plac <- surv.summary.km.ttp.rpsft[surv.summary.km.ttp.rpsft$trtc == "Placebo",]

# ttp KM curve for treatmentA and Placebo-adjusted RPSFTM recensoring
ttp.rpsft <- ggsurvplot(km.ttp.rpsft,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.ttp.rpsft.Len   , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.ttp.rpsft.Plac, mapping = aes(x = time, y = surv,colour="Placebo-RPSFTM-KM",linetype="Placebo-RPSFTM-KM" ), size = 0.75) +
  geom_line(data = transition.12.itt.trtA , aes(x = Time, y = P12,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.12.rpsft.placebo , aes(x = Time, y = P12 , colour = "Placebo-RPSFTM-MSM", linetype = "Placebo-RPSFTM-MSM"),size = 0.75) +
  labs(title = " Time to progression for RPSFT adjustment method \n - 'treatment group method'\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-RPSFTM-KM" = red ,"treatmentA-MSM"= blue, "Placebo-RPSFTM-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM" = 3, "Placebo-RPSFTM-KM" = 1, "Placebo-RPSFTM-MSM" = 3 ))+
  guides(colour = guide_legend(nrow = 2))

ttp.rpsft.at.risk <- ggsurvplot(km.ttp.rpsft, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                                palette =c(blue,red),xlab = "Time (months)",break.time.by = 12,risk.table = TRUE,tables.y.text = FALSE)

# ttp KM curve for treatmentA and Placebo censored at cross.

surv.summary.km.ttp.censoring.at.crossover <- surv_summary(km.ttp.censoring.at.crossover,data = deff.data)

surv.summary.km.ttp.censoring.at.crossover.Len <- surv.summary.km.ttp.censoring.at.crossover[surv.summary.km.ttp.censoring.at.crossover$trtc == "treatmentA",]
surv.summary.km.ttp.censoring.at.crossover.PLac <- surv.summary.km.ttp.censoring.at.crossover[surv.summary.km.ttp.censoring.at.crossover$trtc == "Placebo",]

#ttp KM curve for treatmentA and Placebo-censoring at XO

ttp.censoring.at.crossover <- ggsurvplot(km.ttp.censoring.at.crossover,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.ttp.censoring.at.crossover.Len , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data =  surv.summary.km.ttp.censoring.at.crossover.PLac, mapping = aes(x = time, y = surv,colour="Placebo-Censoring-KM",linetype="Placebo-Censoring-KM" ), size = 0.75) +
  geom_line(data = transition.12.itt.trtA , aes(x = Time, y = P12,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.12.cen.at.cross.placebo , aes(x = Time, y = P12, colour = "Placebo-Censoring-MSM", linetype = "Placebo-Censoring-MSM"),size = 0.75) +
  labs(title = " Time to progression \n Censored at crossover\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-Censoring-KM" = red ,"treatmentA-MSM"= blue, "Placebo-Censoring-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM"= 3, "Placebo-Censoring-KM"=1, "Placebo-Censoring-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

ttp.censoring.at.crossover.at.risk <- ggsurvplot(km.ttp.censoring.at.crossover, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                                                 palette =c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

#############################################################################
# 6.Plot TTD Kaplan-Meier curves of treatmentA vs. each placebo separately
#############################################################################

# TTD KM curve for treatmentA and Placebo-unadjusted(ITT)

surv.summary.km.ttd.itt <- surv_summary(km.ttd.itt,data = deff.data)

surv.summary.km.ttd.itt.Len <- surv.summary.km.ttd.itt[surv.summary.km.ttd.itt$trtc == "treatmentA", ]
surv.summary.km.ttd.itt.Plac <- surv.summary.km.ttd.itt[surv.summary.km.ttd.itt$trtc == "Placebo", ]

ttd.itt <- ggsurvplot(km.ttd.itt,break.time.by = 12, colour= NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.ttd.itt.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.ttd.itt.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = transition.13.itt.trtA , aes(x = Time, y = P13,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.13.itt.placebo , aes(x = Time, y = P13, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Time to death (ITT) \n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))

ttd.itt.at.risk <- ggsurvplot(km.ttd.itt, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                              palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

# TTD KM curve for treatmentA and Placebo-adjusted RPSFTM Re-censoring.

surv.summary.km.ttd.rpsft <- surv_summary(km.ttd.rpsft,data = deff.data)

surv.summary.km.ttd.rpsft.Len <- surv.summary.km.ttd.rpsft[surv.summary.km.ttd.rpsft$trtc == "treatmentA",]
surv.summary.km.ttd.rpsft.Plac <- surv.summary.km.ttd.rpsft[surv.summary.km.ttd.rpsft$trtc == "Placebo",]

ttd.rpsft <- ggsurvplot(km.ttd.rpsft,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data = surv.summary.km.ttd.rpsft.Len   , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.ttd.rpsft.Plac, mapping = aes(x = time, y = surv,colour="Placebo-RPSFTM-KM",linetype="Placebo-RPSFTM-KM" ), size = 0.75) +
  geom_line(data = transition.13.itt.trtA  , aes(x = Time, y = P13,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.13.rpsft.placebo  , aes(x = Time, y = P13 , colour = "Placebo-RPSFTM-MSM", linetype = "Placebo-RPSFTM-MSM"),size = 0.75) +
  labs(title = " Time to death  \n RPSFT adjustment method \n - 'treatment group method'\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-RPSFTM-KM" = red ,"treatmentA-MSM"= blue, "Placebo-RPSFTM-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM" = 3, "Placebo-RPSFTM-KM" = 1, "Placebo-RPSFTM-MSM" = 3 ))+
  guides(colour = guide_legend(nrow = 2))

ttd.rpsft.at.risk <- ggsurvplot(km.ttd.rpsft, combine = TRUE, legend.title = "Treatment", censor=FALSE, xlim = c(0, 180),
                                palette =c(blue,red),xlab = "Time (months)",break.time.by = 12,risk.table = TRUE,tables.y.text = FALSE)

# TTD KM curve for treatmentA and Placebo-censoring at XO

# ttp KM curve for treatmentA and Placebo censored at cross.

surv.summary.km.ttd.censoring.at.crossover <- surv_summary(km.ttd.censoring.at.crossover,data = deff.data)

surv.summary.km.ttd.censoring.at.crossover.Len <- surv.summary.km.ttd.censoring.at.crossover[surv.summary.km.ttd.censoring.at.crossover$trtc == "treatmentA",]
surv.summary.km.ttd.censoring.at.crossover.PLac <- surv.summary.km.ttd.censoring.at.crossover[surv.summary.km.ttd.censoring.at.crossover$trtc == "Placebo",]

#ttd KM curve for treatmentA and Placebo-censoring at XO

ttd.censoring.at.crossover <- ggsurvplot(km.ttd.censoring.at.crossover,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.ttd.censoring.at.crossover.Len , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data =  surv.summary.km.ttd.censoring.at.crossover.PLac, mapping = aes(x = time, y = surv,colour="Placebo-Censoring-KM",linetype="Placebo-Censoring-KM" ), size = 0.75) +
  geom_line(data = transition.13.itt.trtA , aes(x = Time, y = P13,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data =  transition.13.cen.at.cross.placebo, aes(x = Time, y = P13, colour = "Placebo-Censoring-MSM", linetype = "Placebo-Censoring-MSM"),size = 0.75) +
  labs(title = " Time to death \n Censored at crossover\n (CALGB cut-off 1st February 2016)") +
  scale_colour_manual(name= "Curve Type",values = c( "treatmentA-KM" = blue , "Placebo-Censoring-KM" = red ,"treatmentA-MSM"= blue, "Placebo-Censoring-MSM" = red))+
  scale_linetype_manual(name= "Curve Type",values = c( "treatmentA-KM" = 1,  "treatmentA-MSM"= 3, "Placebo-Censoring-KM"=1, "Placebo-Censoring-MSM"=3 ))+
  guides(colour = guide_legend(nrow = 2))

ttd.censoring.at.crossover.at.risk <- ggsurvplot(km.ttd.censoring.at.crossover, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                                                 palette =c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

#############################################################################
# 6.Plot PPS Kaplan-Meier curves of treatmentA vs. each placebo separately
#############################################################################

# PPS KM curve for treatmentA and Placebo-unadjusted (ITT)

surv.summary.km.pps.itt <- surv_summary(km.pps.itt,data = deff.data)

surv.summary.km.pps.itt.Len <- surv.summary.km.pps.itt[surv.summary.km.pps.itt$trtc == "treatmentA", ]
surv.summary.km.pps.itt.Plac <- surv.summary.km.pps.itt[surv.summary.km.pps.itt$trtc == "Placebo", ]

pps.itt <- ggsurvplot(km.pps.itt,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.pps.itt.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.pps.itt.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = transition.23.itt.trtA , aes(x = Time, y = P23,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.23.itt.placebo , aes(x = Time, y = P23, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Post progression survival (ITT) \n (CALGB cut-off 1st February 2016)",caption ="One patient has OS < PFS (1 day). \n This patient has been removed from the analysis of post-progression survival.") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))+
  theme_survminer(font.caption = c(10, "plain", "black"))

# caption = "One patient has PFS > OS, one day later. This patients PPS has been set to 0."

pps.itt.at.risk <- ggsurvplot(km.pps.itt, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                              palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)

# PPS KM curve for treatmentA and Placebo-censored at crossover. 

surv.summary.km.pps.cen_at_cross <- surv_summary(km.pps.cen_at_cross,data = PFS_OS.data)

surv.summary.km.pps.cen_at_cross.Len <- surv.summary.km.pps.cen_at_cross[surv.summary.km.pps.cen_at_cross$trtc == "treatmentA", ]
surv.summary.km.pps.cen_at_cross.Plac <- surv.summary.km.pps.cen_at_cross[surv.summary.km.pps.cen_at_cross$trtc == "Placebo", ]

pps.cen_at_cross <- ggsurvplot(km.pps.cen_at_cross,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.pps.cen_at_cross.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.pps.cen_at_cross.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = transition.23.itt.trtA  , aes(x = Time, y = P23,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.23.cen.at.cross.placebo  , aes(x = Time, y = P23, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Post progression survival  \n Censored at crossover \n (CALGB cut-off 1st February 2016)",caption ="One patient has OS < PFS (1 day). \n This patient has been removed from the analysis of post-progression survival.") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))+
  theme_survminer(font.caption = c(10, "plain", "black"))

#caption = "One patient has PFS > OS, one day later. This patients PPS has been set to 0."

pps.cen_at_cross.at.risk <- ggsurvplot(km.pps.cen_at_cross, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                              palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)


# PPS KM curve for treatmentA and Placebo - RPSFT. 

surv.summary.km.pps.rpsft <- surv_summary(km.pps.rpsft,data = PFS_OS.data)

surv.summary.km.pps.rpsft.Len <- surv.summary.km.pps.rpsft[surv.summary.km.pps.rpsft$trtc == "treatmentA", ]
surv.summary.km.pps.rpsft.Plac <- surv.summary.km.pps.rpsft[surv.summary.km.pps.rpsft$trtc == "Placebo", ]

pps.rpsft <- ggsurvplot(km.pps.rpsft,break.time.by = 12, colour=NA, censor = FALSE, xlim = c(0, 180), xlab= "Time (months)", legend.title = '', conf.int = FALSE) +
  geom_step(data =  surv.summary.km.pps.rpsft.Len  , mapping = aes(x = time, y = surv,colour = "treatmentA-KM",linetype= "treatmentA-KM"), size = 0.75) +
  geom_step(data = surv.summary.km.pps.rpsft.Plac, mapping = aes(x = time, y = surv,colour="Placebo-unadjusted-KM",linetype="Placebo-unadjusted-KM" ), size = 0.75) +
  geom_line(data = transition.23.itt.trtA  , aes(x = Time, y = P23,colour = "treatmentA-MSM",linetype = "treatmentA-MSM"),size = 0.75) +
  geom_line(data = transition.23.rpsft.placebo, aes(x = Time, y = P23, colour = "Placebo-unadjusted-MSM", linetype = "Placebo-unadjusted-MSM"),size = 0.75) +
  labs(title = " Post progression survival  \n RPSFT adjustment method \n - 'treatment group method' \n (CALGB cut-off 1st February 2016)",caption ="One patient has OS < PFS (1 day). \n This patient has been removed from the analysis of post-progression survival.") +
  scale_colour_manual(name="Curve Type",values= c( "treatmentA-KM" = blue , "Placebo-unadjusted-KM" = red ,"treatmentA-MSM"= blue, "Placebo-unadjusted-MSM" = red))+
  scale_linetype_manual(name="Curve Type",values = c( "treatmentA-KM" =1  ,  "Placebo-unadjusted-KM"=1  , "treatmentA-MSM"=3, "Placebo-unadjusted-MSM"=3 )) +
  guides(colour = guide_legend(nrow = 2))+
  theme_survminer(font.caption = c(10, "plain", "black"))

#caption = "One patient has PFS > OS, one day later. This patients PPS has been set to 0."

pps.rpsft.at.risk <- ggsurvplot(km.pps.rpsft, combine = TRUE, legend.title = "Treatment",censor=FALSE,xlim = c(0, 180),
                               palette = c(blue,red), xlab = "Time (months)", break.time.by = 12, risk.table = TRUE, tables.y.text = FALSE)



######################################
#8.Save KM graphs as pdf to directory
######################################
# ITT graphs for OS and PFS


setwd("KM plots")

pdf("./ITT - PFS and OS - MSM.pdf")

gridExtra::grid.arrange(ggplotGrob(PFS.itt$plot),PFS.itt.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(OS.itt$plot),OS.itt.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttp.itt$plot),ttp.itt.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttd.itt$plot),ttd.itt.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(pps.itt$plot),pps.itt.at.risk$table, heights = c(0.75, 0.25))


dev.off()

# RPSFT graphs for OS and PFS

setwd("..../KM plots")

pdf("./RPSFT - PFS and OS - MSM.pdf")

gridExtra::grid.arrange(ggplotGrob(PFS.rpsft$plot),PFS.rpsft.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(OS.rpsft$plot),OS.rpsft.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttp.rpsft$plot),ttp.rpsft.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttd.rpsft$plot),ttd.rpsft.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(pps.rpsft$plot),pps.rpsft.at.risk$table, heights = c(0.75, 0.25))

dev.off()

# crossover censoring graphs for OS and PFS

setwd("..../KM plots")

pdf("./Censoring at crossover - PFS and OS  - MSM.pdf")

gridExtra::grid.arrange(ggplotGrob(PFS.crossover$plot),PFS.crossover.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(OS.censoring.at.crossover$plot),OS.censoring.at.crossover.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttp.censoring.at.crossover$plot),ttp.censoring.at.crossover.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(ttd.censoring.at.crossover$plot),ttd.censoring.at.crossover.at.risk$table, heights = c(0.75, 0.25))

gridExtra::grid.arrange(ggplotGrob(pps.cen_at_cross$plot),pps.cen_at_cross.at.risk$table, heights = c(0.75, 0.25))

dev.off()

#PFS and OS combined graphs

#setwd("..../KM plots")

#pdf("./Combined - PFS and OS.pdf")

#gridExtra::grid.arrange(ggplotGrob(all.PFS$plot),all.PFS.at.risk$table, heights = c(0.75, 0.25))

#gridExtra::grid.arrange(ggplotGrob(all.OS$plot),all.OS.at.risk$table, heights = c(0.75, 0.25))

#dev.off()

