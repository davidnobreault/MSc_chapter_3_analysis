# Load data and clean
##
## This script we load the data and clean it.

## David Breault 2019-05-15

# prep ----
rm(list = ls())

setwd("C:/Users/david/Documents/Google Drive/projects/chapter_3_analysis")

library(tidyverse)
library(lubridate)

## loading data
dat.nec <- read.csv("data/marten_necropsy.csv", header = TRUE, sep = ",")
dat.bone <- read.csv("data/marten_bone_isotopes.csv", header = TRUE, sep = ",")


# cleaning data ----
summary(dat.nec)
