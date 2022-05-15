#######  CRISPRv - Ultrafast quantitative genetic profiling of cellular phenotypes using CRISPR-Cas9 based gene essentiality  #######
#######  by Yuchen Cheng, Jingyuan Chen and Mikael Bjorklund  #######

library(lmtest) #needed for homoscedasticity test


# Select data which has both CERES gene essentiality scores and phenotype values linking via DepMap_ID
# Import quantitative phenotype with cell lines described with DepMap_ID in the first column and phenotype values in the next one
# Here we are using cell size as an example.
Phenotype<-read.csv(file="/Users/chenghui/Movies/lab/Mikael/virtual_screening/CRISPRv/cell_size.csv", header=TRUE, sep=",")

# Import gene essentiality data downloaded from https://depmap.org/portal/download/
CERES<-read.csv(file="/Users/chenghui/Movies/lab/Mikael/2020 summer/sgRNAKO_drug_target/Achilles_gene_effect.csv", header=TRUE, sep=",")

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

# plot GPX4
colnames(ddd)[1]='CellVol'
colnames(ddd)[2]='gpx4'
model=lm(gpx4~CellVol,data=ddd)
new_mpg <- seq(min(ddd$CellVol, max(ddd$CellVol), 0.01))
pred_wt <- data.frame(predict(model, newdata = data.frame(CellVol = new_mpg),
                              interval = "confidence"), 
                      new_mpg = new_mpg)
print(head(pred_wt))

ggplot(ddd,aes(x=CellVol,y=gpx4))+
  geom_point()+
  #geom_abline(slope = -1.112,intercept = 3.015,color='red',size=2)+
  xlab('Log10 Cell Volume (pl)')+ylab('GPX4 CERES Score')+
  geom_smooth(method = 'lm')+
  geom_text_repel(label=ddd$`Phenotype_common$CellLine`)+
  theme_classic()+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=15),
        axis.title=element_text(size=16),axis.text = element_text(size=12) ,legend.text = element_text(size=12),
        legend.title = element_text(size=12))


ddd2=cbind(a,b$ATP5F1D..513.,Phenotype_common$CellLine)
colnames(ddd2)[1]='CellVol'
colnames(ddd2)[2]='atp5f1d'
model=lm(atp5f1d~CellVol,data=ddd2)
new_mpg <- seq(min(ddd$CellVol, max(ddd$CellVol), 0.01))
pred_wt <- data.frame(predict(model, newdata = data.frame(CellVol = new_mpg),
                              interval = "confidence"), 
                      new_mpg = new_mpg)
print(head(pred_wt))

ggplot(ddd2,aes(x=CellVol,y=atp5f1d))+
  geom_point()+
  #geom_abline(slope = -1.112,intercept = 3.015,color='red',size=2)+
  xlab('Log10 Cell Volume (pl)')+ylab('ATP5F1D CERES Score')+
  geom_smooth(method = 'lm')+
  geom_text_repel(label=ddd2$`Phenotype_common$CellLine`)+
  theme_classic()+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=15),
        axis.title=element_text(size=16),axis.text = element_text(size=12) ,legend.text = element_text(size=12),
        legend.title = element_text(size=12))

ddd2=cbind(a,b$MYCN..4613.,Phenotype_common$CellLine)
colnames(ddd2)[1]='CellVol'
colnames(ddd2)[2]='mycn'
model=lm(mycn~CellVol,data=ddd2)
new_mpg <- seq(min(ddd$CellVol, max(ddd$CellVol), 0.01))
pred_wt <- data.frame(predict(model, newdata = data.frame(CellVol = new_mpg),
                              interval = "confidence"), 
                      new_mpg = new_mpg)
print(head(pred_wt))

ggplot(ddd2,aes(x=CellVol,y=mycn))+
  geom_point()+
  #geom_abline(slope = -1.112,intercept = 3.015,color='red',size=2)+
  xlab('Log10 Cell Volume (pl)')+ylab('MYCN CERES Score')+
  geom_smooth(method = 'lm')+
  geom_text_repel(label=ddd2$`Phenotype_common$CellLine`)+
  theme_classic()+
  theme(axis.title.x =element_text(size=16), axis.title.y=element_text(size=15),
        axis.title=element_text(size=16),axis.text = element_text(size=12) ,legend.text = element_text(size=12),
        legend.title = element_text(size=12))
