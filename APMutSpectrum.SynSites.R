###################################
###### 
###################################

############# GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
rm(list=ls(all=TRUE))
GT = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ Syn mut
SynNuc = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GcAtSkewNucCont.csv', header = TRUE)
nrow(SynNuc)

############ merge (work only with mammals since GT is only for mammals)
SynNucGT = merge(GT,SynNuc)
length(unique(SynNucGT$Species))  # 705 species

########### why ND5 and ND2 are almost absent???? 
table(SynNucGT$Gene)
# TP6 ATP8 COX1 COX2 COX3 CytB  ND1  ND2  ND3  ND4 ND4L  ND5  ND6 
# 706  706  706  704  706  706  703    3  706  705  706    3  706 

SynNucGT = SynNucGT[SynNucGT$Gene != 'ND2',]
SynNucGT = SynNucGT[SynNucGT$Gene != 'ND5',]

########### derive fractions
SynNucGT$FrA = SynNucGT$A / SynNucGT$SitesNumber
SynNucGT$FrT = SynNucGT$T / SynNucGT$SitesNumber
SynNucGT$FrG = SynNucGT$G / SynNucGT$SitesNumber
SynNucGT$FrC = SynNucGT$C / SynNucGT$SitesNumber

########### question 1: which nucleotides better correlate with GT: GT = -T + C (in line with our mutational spectrum result that T->C correlates with generation time)
AGG = aggregate(list(SynNucGT$FrA,SynNucGT$FrT,SynNucGT$FrG,SynNucGT$FrC), by = list(SynNucGT$Species,SynNucGT$GenerationLength_d), FUN = mean)
names(AGG) = c('Species','GenerationLength_d','FrA','FrT','FrG','FrC')
par(mfrow=c(2,2))
plot(log2(AGG$GenerationLength_d),AGG$FrA, col = 'gray', ylim = c(0,1), xlab = 'log2(GT)')
plot(log2(AGG$GenerationLength_d),AGG$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)')
plot(log2(AGG$GenerationLength_d),AGG$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)')
plot(log2(AGG$GenerationLength_d),AGG$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)')

cor.test(log2(AGG$GenerationLength_d),AGG$FrA) # nothing
cor.test(log2(AGG$GenerationLength_d),AGG$FrT)
cor.test(log2(AGG$GenerationLength_d),AGG$FrG)
cor.test(log2(AGG$GenerationLength_d),AGG$FrC)

A <- lm(log2(AGG$GenerationLength_d) ~ AGG$FrT + AGG$FrG + AGG$FrC)
summary(A)

A <- lm(log2(AGG$GenerationLength_d) ~ scale(AGG$FrT) + scale(AGG$FrG) + scale(AGG$FrC))
summary(A)

########### question 2: which genes better correlate with GT (why T in ATP6,COX3 and ND4 do not correlate with GT and high absolute value - fast replication, no tRNA before them?)
## T is negatively and C is positively (T->C)
VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2') # ND1 ND2 ND6
length(VecOfGenes)
pdf('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/SynNucCont1.pdf')
par(mfrow=c(2,2))
for (i in 1:length(VecOfGenes))
{ # i = 1
  OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
  main = VecOfGenes[i]
  plot(log2(OneGene$GenerationLength_d),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrT, col = 'blue', ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrG, col = 'green', ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  plot(log2(OneGene$GenerationLength_d),OneGene$FrC, col = 'cyan', ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
}
dev.off()

VecOfGenes = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','CytB') # ND1 ND2 ND6

pdf('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/SynNucCont2A.pdf', height = 30, width = 100)
par(mfrow=c(4,9))
  for (i in 1:length(VecOfGenes))
  { # i = 1
    OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
    main = VecOfGenes[i]
    plot(log2(OneGene$GenerationLength_d),OneGene$FrA, col = 'gray', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  }
  for (i in 1:length(VecOfGenes))
  { # i = 1
    OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
    main = VecOfGenes[i]
    plot(log2(OneGene$GenerationLength_d),OneGene$FrT, col = 'blue', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  }
  for (i in 1:length(VecOfGenes))
  { # i = 1
    OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
    main = VecOfGenes[i]
    plot(log2(OneGene$GenerationLength_d),OneGene$FrG, col = 'green', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  }
  for (i in 1:length(VecOfGenes))
  { # i = 1
    OneGene = SynNucGT[SynNucGT$Gene == VecOfGenes[i],]
    main = VecOfGenes[i]
    plot(log2(OneGene$GenerationLength_d),OneGene$FrC, col = 'cyan', main = main, ylim = c(0,1), xlab = 'log2(GT)') # a bit negative
  }
dev.off()
  
  


cor.test(SynNucGT[SynNucGT$Gene == 'CytB',]$FrA,SynNucGT[SynNucGT$Gene == 'CytB',]$GenerationLength_d, method = 'spearman') # a bit negative
cor.test(SynNucGT[SynNucGT$Gene == 'CytB',]$FrT,SynNucGT[SynNucGT$Gene == 'CytB',]$GenerationLength_d, method = 'spearman') # good negative
cor.test(SynNucGT[SynNucGT$Gene == 'CytB',]$FrG,SynNucGT[SynNucGT$Gene == 'CytB',]$GenerationLength_d, method = 'spearman') # positive
cor.test(SynNucGT[SynNucGT$Gene == 'CytB',]$FrC,SynNucGT[SynNucGT$Gene == 'CytB',]$GenerationLength_d, method = 'spearman') # good positive


plot()


################################
################## COMPARE CLASSES
################################

rm(list=ls(all=TRUE))

############ Syn mut
SynNuc = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GcAtSkewNucCont.csv', header = TRUE)
nrow(SynNuc)

table(SynNuc$Gene) # why ND5 and ND2 are almost absent???? 
SynNuc = SynNuc[SynNuc$Gene != 'ND2',]
SynNuc = SynNuc[SynNuc$Gene != 'ND5',]

### derive fractions
SynNuc$FrA = SynNuc$A / SynNuc$SitesNumber
SynNuc$FrT = SynNuc$T / SynNuc$SitesNumber
SynNuc$FrG = SynNuc$G / SynNuc$SitesNumber
SynNuc$FrC = SynNuc$C / SynNuc$SitesNumber

### Classes from AnAGe
AnAge = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/anage_data.txt', quote = '', sep = '\t',header = TRUE)
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep = '_')
AnAge = AnAge[,grepl("Species|Class", names(AnAge))]
AnAge = unique(AnAge)

### merge
SynNuc = merge(SynNuc,AnAge)
table(SynNuc$Class)
# Actinopterygii           Amphibia               Aves Cephalaspidomorphi     Chondrichthyes           Mammalia           Reptilia      Sarcopterygii 
#  2408                361               2093                 58                407               5406               1097                 44 

## question: cold-blooded versus warmblooded: G>A

par(mfrow=c(2,2))
boxplot(SynNuc[SynNuc$Class == 'Actinopterygii',]$FrA,SynNuc[SynNuc$Class == 'Amphibia',]$FrA,SynNuc[SynNuc$Class == 'Reptilia',]$FrA,SynNuc[SynNuc$Class == 'Mammalia',]$FrA,SynNuc[SynNuc$Class == 'Aves',]$FrA, notch = TRUE, ylab = 'FrA', names=c('act','amph','rep','mam','birds'))
boxplot(SynNuc[SynNuc$Class == 'Actinopterygii',]$FrT,SynNuc[SynNuc$Class == 'Amphibia',]$FrT,SynNuc[SynNuc$Class == 'Reptilia',]$FrT,SynNuc[SynNuc$Class == 'Mammalia',]$FrT,SynNuc[SynNuc$Class == 'Aves',]$FrT, notch = TRUE, ylab = 'FrT', names=c('act','amph','rep','mam','birds'))
boxplot(SynNuc[SynNuc$Class == 'Actinopterygii',]$FrG,SynNuc[SynNuc$Class == 'Amphibia',]$FrG,SynNuc[SynNuc$Class == 'Reptilia',]$FrG,SynNuc[SynNuc$Class == 'Mammalia',]$FrG,SynNuc[SynNuc$Class == 'Aves',]$FrG, notch = TRUE, ylab = 'FrG', names=c('act','amph','rep','mam','birds'))
boxplot(SynNuc[SynNuc$Class == 'Actinopterygii',]$FrC,SynNuc[SynNuc$Class == 'Amphibia',]$FrC,SynNuc[SynNuc$Class == 'Reptilia',]$FrC,SynNuc[SynNuc$Class == 'Mammalia',]$FrC,SynNuc[SynNuc$Class == 'Aves',]$FrC, notch = TRUE, ylab = 'FrC', names=c('act','amph','rep','mam','birds'))






