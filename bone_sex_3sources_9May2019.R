#set library#
.libPaths("C:/R/library")
#load the SIAR package#
library('siar')
#set working directory#
setwd("C:/R/data")
#read in data as matrices and run SIAR
data<-read.table('sex_bone_1may2019.txt',header=TRUE)
sources<-read.table('bone_sex_prey_3sources_9may2019.txt',header=TRUE)
tef<-read.table('tef_summer_fall_3sources_8may2019_v01.txt',header=TRUE)
conc<-read.table('conc_dep_bone_sex_3sources_9may2019.txt',header=TRUE)
sex_bone_1may2019<-siarmcmcdirichletv4(data,sources,tef,conc,500000,50000)
#plot data
siarplotdata(sex_bone_1may2019)
#plot posterior distributions with correlations between source estimates
siarmatrixplot(sex_bone_1may2019)
#graph histograms of proportion densities by group
siarhistograms(sex_bone_1may2019)
#plot source proportions as box plots by group
siarproportionbygroupplot(sex_bone_1may2019)
#plot group proportions by source as box plots
siarproportionbysourceplot(sex_bone_1may2019,
                           prn=TRUE,grp=3,probs=c(5,25,75,95),xspc=0.3)
#view raw output data,
out<-sex_bone_1may2019$output
fix(out)
#summarize model output
library('plyr')
library('reshape2')
melted <- melt(out)
ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))
