###################################
###### 
###################################

rm(list=ls(all=TRUE))
wd <- getwd()
setwd(wd)

############# GENERATION LENGTH FOR ALL MAMMALS

GT = read.table("./Body/1_Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ SYN MUT
SynNuc = read.table("./Body/2_Derived/GcAtSkewNucContBig.csv", header = TRUE)
SynNuc$FrA = SynNuc$A / SynNuc$SitesNumber
SynNuc$FrT = SynNuc$T / SynNuc$SitesNumber
SynNuc$FrG = SynNuc$G / SynNuc$SitesNumber
SynNuc$FrC = SynNuc$C / SynNuc$SitesNumber

############ merge (work only with mammals since GT is only for mammals)
SynNucGT = merge(GT,SynNuc)
table()
length(unique(SynNucGT$Species))  # 705 species
table(SynNucGT$TAXON)

########### question 1: which nucleotides better correlate with GT: log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC) (in line with our mutational spectrum result that T->C correlates with generation time)
AGG = aggregate(list(SynNucGT$FrA,SynNucGT$FrT,SynNucGT$FrG,SynNucGT$FrC), by = list(SynNucGT$Species,SynNucGT$GenerationLength_d), FUN = mean)
names(AGG) = c('Species','GenerationLength_d','FrA','FrT','FrG','FrC')

###### start from pairwise correlations and go to multiple linear model:

cor.test(log2(AGG$GenerationLength_d),AGG$FrA) # a bit negative, but very weak
cor.test(log2(AGG$GenerationLength_d),AGG$FrT) # negative
cor.test(log2(AGG$GenerationLength_d),AGG$FrG) # a bit positive
cor.test(log2(AGG$GenerationLength_d),AGG$FrC) # positive

A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrT + AGG$FrG + AGG$FrC)
summary(A)

A <- lm(log2(AGG$GenerationLength_d) ~ scale(AGG$FrT) + scale(AGG$FrG) + scale(AGG$FrC))
summary(A)
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    11.06737    0.03694 299.586  < 2e-16 ***
#  scale(AGG$FrT) -0.28960    0.04532  -6.390 3.02e-10 ***
#  scale(AGG$FrG) -0.06617    0.04361  -1.517     0.13    
# scale(AGG$FrC)  0.33493    0.04346   7.706 4.45e-14 ***
# Residual standard error: 0.9816 on 702 degrees of freedom
# Multiple R-squared:  0.2113,	Adjusted R-squared:  0.208 
# F-statistic: 62.71 on 3 and 702 DF,  p-value: < 2.2e-16
# log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC)

# plot it:

pdf("./Body/4_Figures/APMutSpectrum.SynSites.MutagensWithinMammals.R.01.pdf", width = 50, height = 30)
par(mfrow=c(2,2),oma = c(0, 0, 2, 0),cex.main = 2, cex.lab = 2)
plot(log2(AGG$GenerationLength_d),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrA'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrT'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrG'); abline(h =0.4, lt = 2, col = 'red')
plot(log2(AGG$GenerationLength_d),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)', main = 'FrC'); abline(h =0.4, lt = 2, col = 'red')
mtext("log2(GT) = 11 - 0.29*scale(FrT) + 0.33*scale(FrC)", outer = TRUE, cex = 1.5)


########### question 2: which genes better correlate with GT (why T in ATP6,COX3 and ND4 do not correlate with GT and high absolute value - fast replication, no tRNA before them?)
## T is negatively and C is positively (T->C)
VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2') # ND1 ND2 ND6
length(VecOfGenes)
par(mfcol=c(4,13))
for (i in 1:length(VecOfGenes))
{ # i = 1
  OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
  main = VecOfGenes[i]
  plot(log2(OneGene$GenerationLength_d),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrT, col = 'blue', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrG, col = 'green', main = main, ylim = c(0,0.6), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrC, col = 'cyan', main = main, ylim = c(0,0.6), xlab = 'log2(GT)') # a bit negative
}
dev.off()

