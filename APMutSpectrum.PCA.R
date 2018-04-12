###################################
### PCA and analysis of the first and second principal components from ecological point of view (temperature, longevity)
###################################
# 04.04.2018; 12.04.2018

####################
#### PCA FOR MAMMALS 
####################

rm(list=ls(all=TRUE))

user = 'Kostya'
if (user == 'Kostya') {setwd('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/')}

###### LOAD MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'
METHOD = 'MAXLIKELIHOOD'

if (METHOD == 'MAXLIKELIHOOD') {MUT = read.table('2_DERIVED/Table_fixed.txt', header = TRUE)}
if (METHOD == 'PARSIMONY')     {MUT = read.table('2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}

##### FILTER: only mammals.  THIS IS STUPID TO TAKE MAMMALIAN LIST FROM GENERATION TIME DATA BASE. WE NEED TAXONOMY DIRECTLY FROM PIPELINE!!! ALYA, KOSTYA
GenerTime = read.table('1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerTime$Species = gsub(' ','_',GenerTime$Scientific_name)
ListOfMammals = GenerTime$Species; length(ListOfMammals) # 5426 - at least will use it to work with mammals only
MUT = MUT[MUT$Species %in% ListOfMammals,]

##### FILTER: normal substitutions
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')
nrow(MUT)
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
MUT = MUT[MUT$Gene == 'CytB',]

##### FILTER: species with many substitutions
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= quantile(AGG1$x,0.15),]$Group.1; length(ListOfSpeciesWithManySubst) # 
MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]

###### DERIVE MUTATIONAL SPECTRUM:
### NORMALIZATION of the 'NumberOfSynMutPerSpecies' by ancestral nucleotide count in the third position of four-fold synonymous substitutions:
NUC = read.table('2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less than before !!!! WHY???? for some species we don't have nucleotide count. Why???? Alya!!!!????

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

### COUNT THE TOTAL NUMBER OF NORMALIZED MUTATIONS PER SPECIES
# create a dataset with total number of mutations and with all 12 types of substitutions
AggTotalMutSpectrumPerSpecies = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum); 
names(AggTotalMutSpectrumPerSpecies) = c('Species','TotalMutRate')
MutTypes = data.frame(VecOfNormalSubstitutions); names(MutTypes)=c('MutType')
nrow(AggTotalMutSpectrumPerSpecies)
AggTotalMutSpectrumPerSpecies = merge(AggTotalMutSpectrumPerSpecies,MutTypes)
nrow(AggTotalMutSpectrumPerSpecies)

# count total number of all types of substitutions and merge with dataset above
AggTotalMutSpectrumPerSpeciesPerMutType = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species,MUT$Subs), FUN = sum); 
names(AggTotalMutSpectrumPerSpeciesPerMutType) = c('Species','MutType','MutTypeRate')
nrow(AggTotalMutSpectrumPerSpeciesPerMutType) # 
ALL = merge(AggTotalMutSpectrumPerSpecies,AggTotalMutSpectrumPerSpeciesPerMutType, by = c('Species','MutType'), all.x=TRUE)
ALL[is.na(ALL)]<-0   # some classes of substitutions are absent from some animals => they are transformed to zeroes. 
ALL$Fraction = ALL$MutTypeRate/ALL$TotalMutRate; summary(ALL$Fraction)

###### create matrix for PCA:
AT = ALL[ALL$MutType == 'A_T',]; AT = AT[c(1,5)]; names(AT) = c('Species','AT'); AT = AT[order(AT$Species),]
AG = ALL[ALL$MutType == 'A_G',]; AG = AG[c(1,5)]; names(AG) = c('Species','AG'); AG = AG[order(AG$Species),]
AC = ALL[ALL$MutType == 'A_C',]; AC = AC[c(1,5)]; names(AC) = c('Species','AC'); AC = AC[order(AC$Species),]
TA = ALL[ALL$MutType == 'T_A',]; TA = TA[c(1,5)]; names(TA) = c('Species','TA'); TA = TA[order(TA$Species),]
TG = ALL[ALL$MutType == 'T_G',]; TG = TG[c(1,5)]; names(TG) = c('Species','TG'); TG = TG[order(TG$Species),]
TC = ALL[ALL$MutType == 'T_C',]; TC = TC[c(1,5)]; names(TC) = c('Species','TC'); TC = TC[order(TC$Species),]
CA = ALL[ALL$MutType == 'C_A',]; CA = CA[c(1,5)]; names(CA) = c('Species','CA'); CA = CA[order(CA$Species),]
CG = ALL[ALL$MutType == 'C_G',]; CG = CG[c(1,5)]; names(CG) = c('Species','CG'); CG = CG[order(CG$Species),]
CT = ALL[ALL$MutType == 'C_T',]; CT = CT[c(1,5)]; names(CT) = c('Species','CT'); CT = CT[order(CT$Species),]
GA = ALL[ALL$MutType == 'G_A',]; GA = GA[c(1,5)]; names(GA) = c('Species','GA'); GA = GA[order(GA$Species),]
GC = ALL[ALL$MutType == 'G_C',]; GC = GC[c(1,5)]; names(GC) = c('Species','GC'); GC = GC[order(GC$Species),]
GT = ALL[ALL$MutType == 'G_T',]; GT = GT[c(1,5)]; names(GT) = c('Species','GT'); GT = GT[order(GT$Species),]
MATRIX = cbind(AT,AG[,2],AC[,2],TA[,2],TG[,2],TC[,2],CA[,2],CG[,2],CT[,2],GA[,2],GC[,2],GT[,2]); names(MATRIX) = c('Species','AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT')

### scale 12 substitutions. 
# If I don't scale PCA1 strongly depends on GA - probably this is just because GA is very common and thus fluctuations are very important by absolute values - not by relative...

### !!!!!!!!! DO IT LIKE BELOW OR LIKE THIS: df.pca <- prcomp(as.matrix(df.rw), center = TRUE, scale. = TRUE) 

#for (i in 2:13)
#{ # i =2
#  summary(MATRIX[,i])
#  MATRIX[,i] = as.numeric(scale(MATRIX[,i]))
#  summary(as.numeric(MATRIX[,i]))
#}

###### PCA 
row.names(MATRIX)=MATRIX$Species
head(MATRIX)
matrix = MATRIX[,c(2:13)]
PCA = prcomp(matrix, center = TRUE, scale. = TRUE) # !!! or scale above!! discuss - it seems identical things
print(PCA)  # PC1: +GA; PC2: -AG; PC3: +TC 
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]
MATRIX$Pca5 = PCA$x[,5]
MATRIX$Pca6 = PCA$x[,6]
MATRIX$Pca7 = PCA$x[,7]
MATRIX$Pca8 = PCA$x[,8]
MATRIX$Pca9 = PCA$x[,9]
MATRIX$Pca10 = PCA$x[,10]
MATRIX$Pca11 = PCA$x[,11]
MATRIX$Pca12 = PCA$x[,12]

# QUESTION: IF WE USE MAMMALIAN SPECIES TO GET PCA (TO UNDERSTAND RULES OF PCA), TO SAVE ALL THESE LINEAR COMBINATIONS AND LATER TO USE THESE COMBINATIONS TO RECONSTRUCT PCs FOR OTHER MAMMALS, CHORDATA AND EVEN CANCERS!?
# FOR ME IT IS MEANINGFUL. ANY COMMENTS?

####################
#### FIND ECOLOGICAL INTERPRETATIONS OF PCs
####################

###### EXTREMES
MATRIX[MATRIX$Pca1 < quantile(MATRIX$Pca1,0.05),]$Species 
MATRIX[MATRIX$Pca1 > quantile(MATRIX$Pca1,0.95),]$Species 
MATRIX[MATRIX$Pca2 < quantile(MATRIX$Pca2,0.05),]$Species 
MATRIX[MATRIX$Pca2 > quantile(MATRIX$Pca2,0.95),]$Species 
MATRIX[MATRIX$Pca3 < quantile(MATRIX$Pca3,0.05),]$Species 
MATRIX[MATRIX$Pca3 > quantile(MATRIX$Pca3,0.95),]$Species   

###### GENERATION LENGTH
GenerTime = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GenerTime$Species = gsub(' ','_',GenerTime$Scientific_name)
GenerTime = GenerTime[,c(11,13)]

Test = merge(MATRIX,GenerTime) 
cor.test(Test$Pca1,Test$GenerationLength_d, method = 'spearman') # no
cor.test(Test$Pca2,Test$GenerationLength_d, method = 'spearman') # super negative (-AG) -0.34
cor.test(Test$Pca3,Test$GenerationLength_d, method = 'spearman') # positive (+TC)       +0.26
cor.test(Test$Pca4,Test$GenerationLength_d, method = 'spearman') # no
cor.test(Test$Pca5,Test$GenerationLength_d, method = 'spearman') # no

####### AnAge 
AnAge = read.table('1_RAW/anage_data.txt', header = TRUE, sep = '\t')
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = AnAge[AnAge$Class == 'Mammalia',];  
names(AnAge)

### Adult.weight..g.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Adult.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Maximum.longevity..yrs.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Maximum.longevity..yrs."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # a bit positive (but unlikely survive multiple test correction)

### Female.maturity..days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Female.maturity..days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Weaning..days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Weaning..days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Birth.weight..g.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Birth.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # a bit positive (unlikely will pass multiple test correction)

### Litters.Clutches.per.year
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Litters.Clutches.per.year"))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca3,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Litter.Clutch.size
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Litter.Clutch.size"))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # positive
cor.test(Test$Pca3,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Weaning.weight..g.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Weaning.weight..g."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # negative
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Metabolic.rate..W. 
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Metabolic.rate..W."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)  
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Temperature..K.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Temperature..K."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

### Growth.rate..1.days.
AnAge1 = AnAge[,(names(AnAge) %in% c("Species","Growth.rate..1.days."))]; AnAge1 = AnAge1[!is.na(AnAge1[,2]),]
Test = merge(MATRIX,AnAge1); nrow(Test)    
cor.test(Test$Pca1,Test[,19], method = 'spearman') # no
cor.test(Test$Pca2,Test[,19], method = 'spearman') # no
cor.test(Test$Pca3,Test[,19], method = 'spearman') # no
cor.test(Test$Pca4,Test[,19], method = 'spearman') # no
cor.test(Test$Pca5,Test[,19], method = 'spearman') # no

###### HETEROTERMS
Hib = read.table('1_RAW/hibernation/tornew5.csv', header = TRUE, sep = ',')
Hib$Species = gsub(" ",'_',Hib$Species)
HibOnlySpecies = Hib[Hib$Type == 'HIB',]$Species; length(HibOnlySpecies)
DtOnlySpecies = Hib[Hib$Type == 'DT',]$Species; length(DtOnlySpecies)
HibAndDtSpecies = Hib[Hib$Type == 'HIB' || Hib$Type == 'DT',]$Species; length(HibAndDtSpecies)
boxplot(MATRIX[MATRIX$Species %in% HibOnlySpecies,]$Pca1,MATRIX[!MATRIX$Species %in% HibOnlySpecies,]$Pca1, notch = TRUE, names = c('HibOnlySpecies','Other'), ylab = 'PC1'); # !!!! HibSpecies have a bit lower Pc1
boxplot(MATRIX[MATRIX$Species %in% DtOnlySpecies,]$Pca1,MATRIX[!MATRIX$Species %in% DtOnlySpecies,]$Pca1, notch = TRUE, names = c('HibOnlySpecies','Other'), ylab = 'PC1');   # nothing
# control for body mass? they are small! # other tests - decreased BMR as compared to body mass - compare strong outliers (cold and hot)...; animals, which live more than should according to bosy mass (naked mole rat)

###### MARSUPIALS, BATS... AGAIN NEED NORMAL TAXONOMY !!!!!!!!
Marsupials = c(AnAge[AnAge$Order == 'Diprotodontia' | AnAge$Order == 'Didelphimorphia' | AnAge$Order == 'Dasyuromorphia',]$Species)
Placental = c(AnAge[AnAge$Order == 'Artiodactyla' | AnAge$Order == 'Cetacea' | AnAge$Order == 'Carnivora',]$Species)
boxplot(MATRIX[MATRIX$Species %in% Marsupials,]$Pca1,MATRIX[MATRIX$Species %in% Placental,]$Pca1, notch = TRUE); # !!!! HibSpecies have a bit lower Pc1

###### FIGURES:

MATRIX = merge(MATRIX,GenerTime)
MATRIX = MATRIX[order(MATRIX$GenerationLength_d),]
MATRIX$Col = c(rep('green',150),rep('gray',187),rep('red',150))
summary(PCA)
print(PCA)

pdf('4_FIGURES/PCA.pdf', width = 14, height = 14)
par(mfcol=c(2,3))
summary(PCA)
#plot(PCA)
plot(MATRIX$Pca1,MATRIX$Pca2, col = MATRIX$Col)
plot(MATRIX$Pca2,MATRIX$Pca3, col = MATRIX$Col)
# plot(PCA$x[,1],MATRIX$GenerationLength_d); cor.test(PCA$x[,1],MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
plot(MATRIX$Pca2,log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman')
plot(MATRIX$Pca3,log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') 
biplot(PCA, col = c('grey','black'), cex = 0.5)
biplot(PCA, choices=c(2,3), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()

####################
#### FIND DNA POLYMERAZE SIGNATURE (last PCs)!!!!!
####################
# LOGIC IS THE NEXT: first - third components are driven by temperature and ecology. 
# other components are not driven by ecology, physiology... sow we need to subtract effect of the first three PC and get naked signature of DNA polymeraze!!!
# how to do it carefully?!

pdf('4_FIGURES/DnaPolymerazeSignature.pdf', width = 14, height = 14)
biplot(PCA, choices=c(4,5), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()

PCA$x # PC's
PCA$sdev # the eigenvalues (res$sdev) giving information on the magnitude of each PC, 
PCA$rotation # and the loadings (res$rotation).

### IF I WANT TO DECREASE DIMENSIONALITY OF MY DATASET, I TRUNCATE IT, USING ONLY THE MOST IMPORTANT PCs:
pc.use <- 3 
trunc <- PCA$x[,1:pc.use] %*% t(PCA$rotation[,1:pc.use]) # trunc <- res$x[,1:pc.use] %*% t(res$rotation[,1:pc.use])

# add the center (and re-scale) back to data
PCA$scale
trunc <- scale(trunc, center = FALSE , scale=1/PCA$scale)
PCA$center
trunc <- scale(trunc, center = -1 * PCA$center, scale=FALSE)

### NOW - OPPOSITE EXERCISE: I WANT TO USE ONLY 4-12 PCs to reconstruct signature of gamma polymeraze:
start = 4
end = 12
gamma <- PCA$x[,start:end] %*% t(PCA$rotation[,start:end]) 
gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)

barplot(
c(
mean(gamma$AT),
mean(gamma$AG),
mean(gamma$AC),
mean(gamma$TA),
mean(gamma$TG),
mean(gamma$TC),
mean(gamma$CA),
mean(gamma$CG),
mean(gamma$CT),
mean(gamma$GA),
mean(gamma$GC),
mean(gamma$GT))
, names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'))

###### doesn't work!!! scaling and centering ...
### try without scaling:
par(mfrow = c(2,1))
PcaNoScale = prcomp(matrix, center = FALSE, scale. = FALSE) # !!! or scale above!! discuss - it seems identical things
start = 1; end = 3
gamma <- PcaNoScale$x[,start:end] %*% t(PcaNoScale$rotation[,start:end]) 
#gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
#gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)
barplot(c(mean(gamma$AT),mean(gamma$AG),mean(gamma$AC),mean(gamma$TA),mean(gamma$TG),mean(gamma$TC),mean(gamma$CA),mean(gamma$CG),mean(gamma$CT),mean(gamma$GA),mean(gamma$GC),mean(gamma$GT))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'), main = 'EXO')

start = 4; end = 10
gamma <- PcaNoScale$x[,start:end] %*% t(PcaNoScale$rotation[,start:end]) 
#gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
#gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)
barplot(c(mean(gamma$AT),mean(gamma$AG),mean(gamma$AC),mean(gamma$TA),mean(gamma$TG),mean(gamma$TC),mean(gamma$CA),mean(gamma$CG),mean(gamma$CT),mean(gamma$GA),mean(gamma$GC),mean(gamma$GT))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'), , main = 'ENDO')

MatrixNoScale = data.frame(matrix); MatrixNoScale$Species = row.names(MatrixNoScale)
PcaNoScaleDataFrame = data.frame(PcaNoScale$x); PcaNoScaleDataFrame$Species = row.names(PcaNoScaleDataFrame)
NoScale = merge(MatrixNoScale,PcaNoScaleDataFrame) 

#### NAIVE CHECKS BELOW - don't work!!! don't understand!!
A = lm(AT ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$ATRes = residuals(A); 
A = lm(AG ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$AGRes = residuals(A); 
A = lm(AC ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$ACRes = residuals(A); 
A = lm(TA ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$TARes = residuals(A); 
A = lm(TG ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$TGRes = residuals(A); 
A = lm(TC ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$TCRes = residuals(A); 
A = lm(CG ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$CGRes = residuals(A); 
A = lm(CT ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$CTRes = residuals(A); 
A = lm(CA ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$CARes = residuals(A); 
A = lm(GT ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$GTRes = residuals(A); 
A = lm(GC ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$GCRes = residuals(A); 
A = lm(GA ~ PC1 + PC2 + PC3, data = NoScale); summary(A); NoScale$GARes = residuals(A); 
barplot(c(mean(NoScale$ATRes),mean(NoScale$AGRes),mean(NoScale$ACRes),mean(NoScale$TARes),mean(NoScale$TGRes),mean(NoScale$TCRes),mean(NoScale$CARes),mean(NoScale$CGRes),mean(NoScale$CTRes),mean(NoScale$GARes),mean(NoScale$GCRes),mean(NoScale$GTRes))  , names = c('AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT'))
