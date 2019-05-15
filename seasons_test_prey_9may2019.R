#set library#
.libPaths("C:/R/library")
#set working directory#
setwd("C:/R/data")
#check differences in prey SI between seasons
#load prey group data
my_data <- read.delim('seasons_test_squirrel_9may2019.txt')
#check groups are correct
levels(my_data$Group)
#kruskal-wallis test of differences between > 2 groups in C13
kruskal.test(c13 ~ Group, data = my_data)
#pairwise comparisons between group levels
pairwise.wilcox.test(my_data$c13, my_data$Group,
                     p.adjust.method = "bonferroni")
#kruskal-wallis test of differences between > 2 groups in N15
kruskal.test(n15 ~ Group, data = my_data)
#pairwise comparisons between group levels
pairwise.wilcox.test(my_data$n15, my_data$Group,
                     p.adjust.method = "bonferroni")
