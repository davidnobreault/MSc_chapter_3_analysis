#SIAR mixing model using SI consumer data from marten bone tissue, with sex as the grouping variable, and 4 sources
##playing around on GitHub with Angus on May 15th, 2019

#setup----
rm(list = ls())

#set library
.libPaths("C:/R/library")

#load the SIAR package
library('siar')

#set working directory
setwd("C:\Documents\Google Drive\projects\chapter_3_analysis")

#read in data as matrices and run solo version of SIAR
data<-read.table('sex_bone_1may2019.txt',header=TRUE)

sources<-read.table('prey_4sources_all_seasons_13may2019.txt',header=TRUE)

tef<-read.table('tef_summer_fall_4sources_10may2019.txt',header=TRUE)

conc<-read.table('conc_dep_4sources_all_seasons_13may2019.txt',header=TRUE)

bone_sex_13may2019<-siarmcmcdirichletv4(data,sources,tef,conc,500000,50000)

#plot data
siarplotdata(bone_sex_13may2019)

#plot posterior distributions with correlations between source estimates
siarmatrixplot(bone_sex_13may2019)

#graph histograms of proportion densities by group
siarhistograms(bone_sex_13may2019)

#plot source proportions as box plots by group
siarproportionbygroupplot(bone_sex_13may2019)

#plot group proportions by source as box plots
siarproportionbysourceplot(bone_sex_13may2019,
                           prn=TRUE,grp=4,probs=c(5,25,75,95),xspc=0.3)

#view raw output data
out<-bone_sex_13may2019$output
fix(out)

#summarize model output
library('plyr')
library('reshape2')
melted <- melt(out)
ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))
