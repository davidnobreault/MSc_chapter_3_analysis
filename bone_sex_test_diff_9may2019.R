#set library#
.libPaths("C:/R/library")
#set working directory#
setwd("C:/R/data")
#ANOVA and/or two-sample unpaired t-test, to check differences between marten tissues
#load marten tissue data
my_data <- read.delim('sex_bone_test_diff_9may2019.txt')
#check groups are correct
levels(my_data$Code)
#test for normality (small p-value indicates different from normal)
# Shapiro-Wilk normality test for female C13 values
with(my_data, shapiro.test(d13C[Code == "female"]))
# Shapiro-Wilk normality test for male resource c13 values
with(my_data, shapiro.test(d13C[Code == "male"])) 
#test for normality
# Shapiro-Wilk normality test for female n15 values
with(my_data, shapiro.test(d15N[Code == "female"]))
# Shapiro-Wilk normality test for male n15 values
with(my_data, shapiro.test(d15N[Code == "male"])) 
#low p-values indicate non-normal distributions; will do nonparametric test
#unpaired two-sample wilcoxon signed-rank test - are groups different in c13 values?
res <- wilcox.test(d13C ~ Code, data = my_data,
                   exact = FALSE)
res
#unpaired two-sample wilcoxon signed-rank test - are groups different in n15 values?
res <- wilcox.test(d15N ~ Code, data = my_data,
                   exact = FALSE)
res
