#######  CRISPRv - Ultrafast quantitative genetic profiling of cellular phenotypes using CRISPR-Cas9 based gene essentiality  #######
#######  by Yuchen Cheng, Jingyuan Chen and Mikael Bjorklund  #######

library(lmtest) #needed for homoscedasticity test


# Select data which has both CERES gene essentiality scores and phenotype values linking via DepMap_ID
# Import quantitative phenotype with cell lines described with DepMap_ID in the first column and phenotype values in the next one
# Here we are using cell size as an example.
Phenotype<-read.csv(file="cell_size.csv", header=TRUE, sep=",")

# Import gene essentiality data downloaded from https://depmap.org/portal/download/
CERES<-read.csv(file="Achilles_gene_effect.csv", header=TRUE, sep=",")

# Select common cell lines in both phenotype and essentiality dataset
commonIDs<-as.vector(Phenotype$DepMap_ID)
CERES_common<-CERES[CERES$DepMap_ID %in% commonIDs, ] 
commonIDs2<-as.vector(CERES_common$DepMap_ID)
Phenotype_common<-Phenotype[Phenotype$DepMap_ID %in% commonIDs2, ]

# Select data for association analysis
# The code can run multiple phenotypic associations simultaneously by adding more columns with phenotype data
a<-as.data.frame(Phenotype_common[,5]) # Uses Log10CellVolume column as an example
b<-CERES_common[,2:ncol(CERES_common)]

# t statistic value for linear model (significance and directionality of association)
t_stats<-apply(b, 2, function(x) apply(a, 2, function(y) summary(lm(x~y))$coefficients[2,3]))

# slope for linear model (strength of association)
slope<-apply(b, 2, function(x) apply(a, 2, function(y) summary(lm(x~y))$coefficients[2,1]))

# p value for linear model (significance of association, partially redundant with t stat)
p_val<-apply(b, 2, function(x) apply(a, 2, function(y) summary(lm(x~y))$coefficients[2,4]))

# Homoscedasticity test (Breush-Pagan test)
BP_test<-apply(b, 2, function(x) apply(a, 2, function(y) lmtest::bptest(lm(x~y))$p.value))

# Normality of residuals by Shapiro-Wilk test. Note that this test is sensitive to large sample sizes
shapiro_test<-apply(b, 2, function(x) apply(a, 2, function(y) shapiro.test(resid(lm(x~y)))$p.value))

# Write summary 
summary<-cbind(t_stats, slope, p_val, BP_test, shapiro_test) 

write.table(summary, file = "/CRISPRv_tstat_slope_pval_homosced_normality_residuals.txt", sep = "\t")

