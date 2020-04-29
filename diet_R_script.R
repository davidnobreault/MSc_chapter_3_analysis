##These are the analyses carried out for Chapter 3: diet of marten on HG by SI analysis
#performed by D. Breault in 2019 for MSc Thesis

install.packages("extrafont")

install.packages("ggpubr")

install.packages("PMCMR")

# prep ----

#set library
.libPaths("C:/Users/david/Documents/Google Drive/UNBC/HG Marten Project/R/library")

#remove any objectives from previous analysis

rm(list = ls())

#Load packages

library(tidyverse)

library(lubridate)

library(ggpubr)

theme_set(theme_pubr())

library('siar')

library('plyr')

library('reshape2')

library(extrafont)

library('PMCMR')

font_import()
y

loadfonts(device="win")       #Register fonts for Windows bitmap output

fonts()

#set working directory

setwd(("C:/Users/david/Documents/Google Drive/UNBC/HG Marten Project/R/projects/chapter_3_analysis"))

# load marten tissue SI data

dat.tis <- read.csv("data/marten_all_info_12jun2019_v04.csv", header = TRUE, sep = ",")

#plot the hair bottom data----

ggplot(dat.tis, aes(x = bothair_c13 , y = bothair_n15)) +
  geom_point() +
  stat_smooth()

# test linear correlation between C13 and N15 for each tissue type ----
##This is to show if SI signatures are concistent with a marine diet (Pauli et al. 2012)

cor(dat.tis$muscle_n15, dat.tis$muscle_c13, use = "complete.obs")

cor(dat.tis$bone_n15, dat.tis$bone_c13, use = "complete.obs")

cor(dat.tis$bothair_n15, dat.tis$bothair_c13, use = "complete.obs")

cor(dat.tis$tophair_n15, dat.tis$tophair_c13, use = "complete.obs")

##Appears that bottom hair sections show strongest correlation (r = 0.71) (consistent with fall/salmon)

#on SI values from bottom hair, perform linear regression with x = C13 and y = N15

bot.lin <- lm(bothair_n15 ~ bothair_c13, data = dat.tis)

bot.lin

# plot the data with the linear regression line and 95% CI


ggplot(dat.tis, aes(bothair_c13, bothair_n15))  +
  geom_point() +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "??13C", y = "??15N") +
  theme(axis.title.x = element_text(family="Times New Roman", size=11),
        axis.title.y = element_text(family="Times New Roman", size=11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "11"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "11")) +
  scale_x_continuous(breaks = seq(-25, -18, 1)) +
  scale_y_continuous(breaks = seq(5, 13, 1))

#test significance of relationship----

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(bot.lin)

#calculate the 95% CI around the beta coefficients

confint(bot.lin)

#this analysis was carried out by David Breault on May 20, 2019 (and again on June 17, 2019)

#This analysis is to test for linear correlation between muscle SI values and time----
#I did this on June 2nd, 2019

# load marten muscle SI data

dat.mus <- read.csv("data/marten_muscle_time_21jun2019.csv", header = TRUE, sep = ",")

#plot the data----

ggplot(dat.mus, aes(x = days_sf, y = muscle_c13)) +
  geom_point() +
  stat_smooth()

# test linear correlation between C13 and N15 and Time

cor(dat.mus$muscle_c13, dat.mus$days_sf)

cor(dat.mus$muscle_n15, dat.mus$days_sf)

#perform linear regression with x = Time and y = C13, and y = N15

lin.c13 <- lm(muscle_c13 ~ days_sf, data = dat.mus)

lin.c13

lin.n15 <- lm(muscle_n15 ~ days_sf, data = dat.mus)

lin.n15

# plot the data with the linear regression line and 95% CI

ggplot(dat.mus, aes(days_sf, muscle_c13)) +
  geom_point() +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Days since December 3rd", y = "??13C") +
  theme(axis.title.x = element_text(family="Times New Roman", size=11),
        axis.title.y = element_text(family="Times New Roman", size=11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "11"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "11")) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  scale_y_continuous(breaks = seq(-30, -20, 1))

ggplot(dat.mus, aes(days_sf, muscle_n15)) +
  geom_point() +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Days since December 3rd", y = "??15N") +
  theme(axis.title.x = element_text(family="Times New Roman", size=11),
        axis.title.y = element_text(family="Times New Roman", size=11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "11"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "11")) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  scale_y_continuous(breaks = seq(5, 14, 1))

#test significance of relationship----

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.c13)

summary(lin.n15)

#calculate the 95% CI around the beta coefficients

confint(lin.c13)

confint(lin.n15)

#This analysis is to test for differences between seasons within each prey group----
#This was done on May 31, 2019

#compare berries across seasons

dat.berr <- dat.prey %>%
  filter(Source == "berry")

kruskal.test(c13 ~ Season, data = dat.berr)

pairwise.wilcox.test(dat.berr$c13, dat.berr$Season,
                     p.adjust.method = "bonferroni")

kruskal.test(n15 ~ Season, data = dat.berr)

pairwise.wilcox.test(dat.berr$n15, dat.berr$Season,
                     p.adjust.method = "bonferroni")

#compare t.vertebrates across seasons

dat.fauna <- dat.prey %>%
  filter(Source == "t.fauna")

kruskal.test(c13 ~ Season, data = dat.fauna)

pairwise.wilcox.test(dat.fauna$c13, dat.fauna$Season,
                     p.adjust.method = "bonferroni")

kruskal.test(n15 ~ Season, data = dat.fauna)

pairwise.wilcox.test(dat.fauna$n15, dat.fauna$Season,
                     p.adjust.method = "bonferroni")

#compare m.invertebrates across seasons

dat.inv <- dat.prey %>%
  filter(Source == "m.invertebrate")

kruskal.test(c13 ~ Season, data = dat.inv)

pairwise.wilcox.test(dat.inv$c13, dat.inv$Season,
                     p.adjust.method = "bonferroni")

kruskal.test(n15 ~ Season, data = dat.inv)

pairwise.wilcox.test(dat.inv$n15, dat.inv$Season,
                     p.adjust.method = "bonferroni")

#compare salmon across seasons

dat.sal <- dat.prey %>%
  filter(Source == "salmon")

kruskal.test(c13 ~ Season, data = dat.sal)

pairwise.wilcox.test(dat.sal$c13, dat.sal$Season,
                     p.adjust.method = "bonferroni")

kruskal.test(n15 ~ Season, data = dat.sal)

pairwise.wilcox.test(dat.sal$n15, dat.sal$Season,
                     p.adjust.method = "bonferroni")





#This is the analysis of significant differences between 4 prey groups, sampled across all seasons,
#for use in SIAR mixing models of lifetime marten diet, using marten bone tissue C and N ratios

#This is to show my rational for picking the prey groups that I did----

#Load diet item data

dat.prey <- read.csv("data/prey_all_info_c_n_17jun2019.csv", header = TRUE, sep = ",")

#check groups are correct

levels(dat.prey$Source)

##test for normality----

# Shapiro-Wilk normality test for C13 values by prey group

with(dat.prey, shapiro.test(c13[Source == "salmon"]))

with(dat.prey, shapiro.test(c13[Source == "m.invertebrate"]))

with(dat.prey, shapiro.test(c13[Source == "t.fauna"]))

with(dat.prey, shapiro.test(c13[Source == "berry"]))

#violations of normality assumption occurred, therefore move to nonparametric test of differences

#Nonparametric test of differences between groups----

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(c13 ~ Source, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

posthoc.kruskal.nemenyi.test(x = c13, g = Source, dist="Tukey")

#kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(n15 ~ Source, data = dat.prey)

require(PMCMR)

data(dat.prey)

attach(dat.prey)

posthoc.kruskal.nemenyi.test(x = n15, g = Source, dist="Tukey")

#All prey groups are significantly differed from each other in C and N values

#Mixing Models----

#Marten tissues as consumer groupings----

##Run SIAR with consumer table coded by tissue type, and 4 prey groups, all sampled across all
#seasons (only looking at output from bone, for lifetime diet)

#Run model with tissue type as grouping variable
#hair top = 1, hair bottom = 2, muscle = 3, bone = 4

#use marten tissue data 

dat.mar <- read.table('data/marten_tissues_coded_21jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

tis.mod <- siarmcmcdirichletv4(dat.mar,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(tis.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(tis.mod)

#graph histograms of proportion densities by group

siarhistograms(tis.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(tis.mod)

#plot group (tissue type) proportions by source as box plots

siarproportionbysourceplot(tis.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(tis.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(tis.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(tis.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- tis.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))


#This analysis is to determine if tissues differ significantly in SI values----
#I did this on June 3rd, 2019

#load marten tissue data

dat.mar <- read.table('data/marten_tissues_coded_12jun2019.txt',header=TRUE)

# Shapiro-Wilk normality test for C13 and N15 by tissue type

with(dat.mar, shapiro.test(d13C[Code == 1]))

with(dat.mar, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.mar)

#pairwise comparisons between group levels

pairwise.wilcox.test(dat.mar$d13C, dat.mar$Code,
                     p.adjust.method = "bonferroni")

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.mar)

#pairwise comparisons between group levels

pairwise.wilcox.test(dat.mar$d15N, dat.mar$Code,
                     p.adjust.method = "bonferroni")

#remove any objectives from previous analysis

rm(list = ls())

#Salmon density as consumer groupings (June 12th 2019)----

##Run SIAR with consumer SI values (from bone) coded by salmon density----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with salmon density as grouping variable
# =< 0.22 = 1, 0.22 < x =< 0.60 = 2, > 0.60 = 3

#use marten tissue data 

dat.sabon <- read.table('data/marten_coded_salmon_bone_12jun2019.txt',header=TRUE)

#This analysis is to determine if marten from areas with different salmon stream density
#differ significantly in SI values----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.sabon, shapiro.test(d13C[Code == 1]))

with(dat.sabon, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.sabon)

require(PMCMR)

data(dat.sabon)

attach(dat.sabon)

posthoc.kruskal.nemenyi.test(x = d13C, g = Code, dist="Tukey")

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.sabon)

require(PMCMR)

data(dat.sabon)

attach(dat.sabon)

posthoc.kruskal.nemenyi.test(x = d15N, g = Code, dist="Tukey")

##marten bone not significantly different in either isotope between salmon densities

##Run SIAR with consumer SI values (from hair_bottom) coded by salmon density----
#and 4 prey groups, all sampled across all seasons (did this on June 13th, 2019)

#Run model with salmon density as grouping variable
# =< 0.22 = 1, 0.22 < x =< 0.60 = 2, > 0.60 = 3

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.sabot <- read.table('data/marten_coded_salmon_bottom_12jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

sabot.mod <- siarmcmcdirichletv4(dat.sabot,dat.sou,dat.tef,dat.con,500000,50000)

#plot group proportions by source as box plots

siarproportionbysourceplot(sabot.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sabot.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sabot.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sabot.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)
 
#view raw output data

out <- sabot.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if marten from areas with different salmon stream density
#differ significantly in SI values----
#I did this on June 13th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.sabot, shapiro.test(d13C[Code == 1]))

with(dat.sabot, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.sabot)

require(PMCMR)

data(dat.sabot)

attach(dat.sabot)

posthoc.kruskal.nemenyi.test(x = d13C, g = Code, dist="Tukey")

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.sabot)

require(PMCMR)

data(dat.sabot)

attach(dat.sabot)

posthoc.kruskal.nemenyi.test(x = d15N, g = Code, dist="Tukey")

##Run SIAR with consumer SI values (from bone) coded by distance to marine shoreline----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with distance to shoreline as grouping variable
# =< 2577 m = 1, 2577 < x =< 10729 m = 2, > 10729 = 3

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.cobon <- read.table('data/marten_coast_bone_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

cobon.mod <- siarmcmcdirichletv4(dat.cobon,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(cobon.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(cobon.mod)

#graph histograms of proportion densities by group

siarhistograms(cobon.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(cobon.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(cobon.mod,prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cobon.mod,prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cobon.mod,prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cobon.mod,prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- cobon.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if marten (bone) from areas with different distance to marine shoreline
#differ significantly in SI values----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by distance to shoreline

with(dat.cobon, shapiro.test(d13C[Code == 1]))

with(dat.cobon, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.cobon)

require(PMCMR)

data(dat.cobon)

attach(dat.cobon)

posthoc.kruskal.nemenyi.test(x = d13C, g = Code, dist="Tukey")

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.cobon)

require(PMCMR)

data(dat.cobon)

attach(dat.cobon)

posthoc.kruskal.nemenyi.test(x = d15N, g = Code, dist="Tukey")

##Run SIAR with consumer SI values (from hair_bottom) coded by distance to marine shoreline----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with distance to shoreline as grouping variable
# =< 2577 m = 1, 2577 < x =< 10729 m = 2, > 10729 = 3

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.cobot <- read.table('data/marten_coast_bottom_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

cobot.mod <- siarmcmcdirichletv4(dat.cobot,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(cobot.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(cobot.mod)

#graph histograms of proportion densities by group

siarhistograms(cobot.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(cobot.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(cobot.mod,
                           prn=TRUE,grp=4,probs=c(5,25,75,95),xspc=0.3)

#view raw output data

out <- cobot.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if marten from areas with different distance to marine shoreline
#differ significantly in SI values----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by marine distance

with(dat.cobot, shapiro.test(d13C[Code == 1]))

with(dat.cobot, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.cobot)

#pairwise comparisons between group levels

pairwise.wilcox.test(dat.cobot$d13C, dat.cobot$Code,
                     p.adjust.method = "bonferroni")

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.cobot)

#pairwise comparisons between group levels

pairwise.wilcox.test(dat.cobot$d15N, dat.cobot$Code,
                     p.adjust.method = "bonferroni")

##Run SIAR with consumer SI values (from bone) coded by marten age class----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with age class as grouping variable
# Juvenile (age 0) = 1, Adult (age >= 1) = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.age <- read.table('data/marten_age_bone_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

age.mod <- siarmcmcdirichletv4(dat.age,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(age.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(age.mod)

#graph histograms of proportion densities by group

siarhistograms(age.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(age.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(age.mod,
                           prn=TRUE,grp=4,probs=c(5,25,75,95),xspc=0.3)

#view raw output data

out <- age.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if juvenile and adult marten differ significantly in SI values----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.age, shapiro.test(d13C[Code == 1]))

with(dat.age, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.age)

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.age)

##Run SIAR with consumer SI values (from hair_bottom) coded by marten age class----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with age class as grouping variable
# Juvenile (age 0) = 1, Adult (age >= 1) = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.abot <- read.table('data/marten_age_bottom_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

abot.mod <- siarmcmcdirichletv4(dat.abot,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(abot.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(abot.mod)

#graph histograms of proportion densities by group

siarhistograms(abot.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(abot.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(abot.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(abot.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(abot.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(abot.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- abot.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if juvenile and adult marten (hair bottom) differ significantly in SI values----
#I did this on June 13th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.abot, shapiro.test(d13C[Code == 1]))

with(dat.abot, shapiro.test(d15N[Code == 1]))

with(dat.abot, shapiro.test(d13C[Code == 2]))

with(dat.abot, shapiro.test(d15N[Code == 2]))

#check for homogeneity of variances

ftest.d13c <- var.test(d13C ~ Code, data = dat.abot)

ftest.d13c

ftest.d15n <- var.test(d15N ~ Code, data = dat.abot)

ftest.d15n

# Compute t-test

tt.abc <- t.test(d13C ~ Code, data = dat.abot, var.equal = TRUE)

tt.abc

tt.abn <- t.test(d15N ~ Code, data = dat.abot, var.equal = TRUE)

tt.abn


##Run SIAR with consumer SI values (from bone) coded by marten sex----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with sex as grouping variable
# Female = 1, Male = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.sex <- read.table('data/marten_sex_bone_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

sex.mod <- siarmcmcdirichletv4(dat.sex,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(sex.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(sex.mod)

#graph histograms of proportion densities by group

siarhistograms(sex.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(sex.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(sex.mod,
                           prn=TRUE,grp=4,probs=c(5,25,75,95),xspc=0.3)

#view raw output data

out <- sex.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if male and female marten differ significantly in SI values----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.sex, shapiro.test(d13C[Code == 1]))

with(dat.sex, shapiro.test(d15N[Code == 1]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.sex)

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.sex)

##Run SIAR with consumer SI values (from hair_bottom) coded by marten sex----
#and 4 prey groups, all sampled across all seasons (did this on June 12th, 2019)

#Run model with sex as grouping variable
# Female = 1, Male = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.sbot <- read.table('data/marten_sex_bottom_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

sbot.mod <- siarmcmcdirichletv4(dat.sbot,dat.sou,dat.tef,dat.con,500000,50000)

#plot data with TEFs applied

siarplotdata(sbot.mod)

#plot posterior distributions with correlations between source estimates

siarmatrixplot(sbot.mod)

#graph histograms of proportion densities by group

siarhistograms(sbot.mod)

#plot source proportions as box plots by group

siarproportionbygroupplot(sbot.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(sbot.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sbot.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sbot.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(sbot.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- sbot.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))


#This analysis is to determine if male and female marten differ significantly in SI values from hair_bottom----
#I did this on June 12th, 2019

# Shapiro-Wilk normality test for C13 and N15 by sex

with(dat.sbot, shapiro.test(d13C[Code == 1]))

with(dat.sbot, shapiro.test(d15N[Code == 1]))

with(dat.sbot, shapiro.test(d13C[Code == 2]))

with(dat.sbot, shapiro.test(d15N[Code == 2]))

#check for homogeneity of variances

ftest.d13c <- var.test(d13C ~ Code, data = dat.sbot)

ftest.d13c

ftest.d15n <- var.test(d15N ~ Code, data = dat.sbot)

ftest.d15n

# Compute t-test

tt.sbc <- t.test(d13C ~ Code, data = dat.sbot, var.equal = TRUE)

tt.sbc

tt.sbn <- t.test(d15N ~ Code, data = dat.sbot, var.equal = TRUE)

tt.sbn

##Run SIAR with consumer SI values (from hair_bottom) coded by marten body condition----
#and 4 prey groups, all sampled across all seasons (did this on June 13th, 2019)

#Run model with condition as grouping variable
# Poor, Fair, Good = 1, Excellent, Obese = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.cbot <- read.table('data/marten_condition_bottom_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

cbot.mod <- siarmcmcdirichletv4(dat.cbot,dat.sou,dat.tef,dat.con,500000,50000)

#graph histograms of proportion densities by group

siarhistograms(cbot.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(cbot.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbot.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbot.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbot.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- cbot.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if body conditions differ significantly in SI values----
#I did this on June 13th, 2019

# Shapiro-Wilk normality test for C13 and N15 by body condition

with(dat.cbot, shapiro.test(d13C[Code == 1]))

with(dat.cbot, shapiro.test(d15N[Code == 1]))

with(dat.cbot, shapiro.test(d13C[Code == 2]))

with(dat.cbot, shapiro.test(d15N[Code == 2]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.cbot)

##kruskal-wallis test of differences between > 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.cbot)

##Run SIAR with consumer SI values (from bone) coded by marten body condition----
#and 4 prey groups, all sampled across all seasons (did this on June 13th, 2019)

#Run model with condition as grouping variable
# Poor, Fair, Good = 1, Excellent, Obese = 2

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.cbon <- read.table('data/marten_condition_bone_13jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

cbon.mod <- siarmcmcdirichletv4(dat.cbon,dat.sou,dat.tef,dat.con,500000,50000)

#graph histograms of proportion densities by group

siarhistograms(cbon.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- cbon.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if body conditions differ significantly in SI values----
#I did this on June 13th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.cbon, shapiro.test(d13C[Code == 1]))

with(dat.cbon, shapiro.test(d15N[Code == 1]))

with(dat.cbon, shapiro.test(d13C[Code == 2]))

with(dat.cbon, shapiro.test(d15N[Code == 2]))

##kruskal-wallis test of differences between > 2 groups in C13

wilcox.test(d13C ~ Code, data = dat.cbon, alternative = c("less"))

##kruskal-wallis test of differences between > 2 groups in N15

wilcox.test(d15N ~ Code, data = dat.cbon, alternative = c("greater"))

##Run SIAR with consumer SI values (from bone) coded by marten body condition----
#and 4 prey groups, all sampled across all seasons (did this on June 24th, 2019)

#Run model with condition as grouping variable
# Poor, Fair = 1, Good = 2, Excellent, Obese = 3

#remove any objectives from previous analysis

rm(list = ls())

#use marten tissue data 

dat.cbon <- read.table('data/marten_condition_bone_24jun2019.txt',header=TRUE)

#use source table with 4 groups, sampled across all seasons

dat.sou <- read.table('data/prey_groups_mean_12jun2019.txt',header=TRUE)

#use trophic enrichment fractor table with TEF values for 4 sources in this model

dat.tef <- read.table('data/tef_summer_fall_4sources_12jun2019.txt',header=TRUE)

#use mean concentrations (%) of N and C with SDs for 4 sources in this model

dat.con <- read.table('data/conc_dep_12jun2019.txt',header=TRUE)

#create model using uninformative priors (dirichlet distribution), with 500,000 iterations and 5000 burn-ins

cbon.mod <- siarmcmcdirichletv4(dat.cbon,dat.sou,dat.tef,dat.con,500000,50000)

#graph histograms of proportion densities by group

siarhistograms(cbon.mod)

#plot group proportions by source as box plots

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=1, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=2, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=3, probs=c(5,25,75,95), xspc=0.3)

siarproportionbysourceplot(cbon.mod, prn=TRUE, grp=4, probs=c(5,25,75,95), xspc=0.3)

#view raw output data

out <- cbon.mod$output

fix(out)

#summarize model output

melted <- melt(out)

ddply(melted, c("Var2"), summarise,
      mean = mean(value), max = max(value), min = min(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)), count = length(value))

#This analysis is to determine if body conditions differ significantly in SI values----
#I did this on June 13th, 2019

# Shapiro-Wilk normality test for C13 and N15 by salmon density

with(dat.cbon, shapiro.test(d13C[Code == 1]))

with(dat.cbon, shapiro.test(d15N[Code == 1]))

with(dat.cbon, shapiro.test(d13C[Code == 2]))

with(dat.cbon, shapiro.test(d15N[Code == 2]))

##kruskal-wallis test of differences between > 2 groups in C13

kruskal.test(d13C ~ Code, data = dat.cbon)

##kruskal-wallis test of differences between 2 groups in N15

kruskal.test(d15N ~ Code, data = dat.cbon)
