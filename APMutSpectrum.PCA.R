###################################
### PCA and analysis of the first and second principal components from ecological point of view (temperature, longevity)
###################################
# 04.04.2018; 12.04.2018

####################
#### PCA FOR MAMMALS 
####################

rm(list = ls(all.names = TRUE))

# For PCA analysis:
library("irlba")
# library("logisticPCA")  # May be used to explore presence/absence some type of subst. I've add it in case of our zeros-completion of NA.
# For data-wrangling:
library("tidyverse")
library("broom")
# For ploting:
library("gplots")
library("factoextra")
library("corrplot")
library("cowplot")
library("ggrepel")

source("./funDimReduction.R")


# user = 'Kostya'
# if (user == 'Kostya') {setwd('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/')}

###### LOAD MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'
METHOD <- 'MAXLIKELIHOOD'

options(stringsAsFactors = FALSE)
if (METHOD == 'MAXLIKELIHOOD') {MUT <- read_table2('2_Derived/Table_fixed.txt')}
nrow(MUT)  # 506179
length(unique(MUT$Species))  # 2404
# if (METHOD == 'PARSIMONY')     {MUT = read.table('2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}  # not loaded yet

##### FILTER: only mammals.  THIS IS STUPID TO TAKE MAMMALIAN LIST FROM GENERATION TIME DATA BASE. WE NEED TAXONOMY DIRECTLY FROM PIPELINE!!! ALYA, KOSTYA
GenerTime <- read_tsv('1_Raw/GenerationLenghtforMammals.xlsx.txt', 
                      na = c("", "NA", "no information")) %>% 
    mutate(Species = str_replace(Scientific_name, " ", "_")) %>% 
    select(-starts_with("Sources"), -Genus, -Scientific_name, -Data_AFR)
ListOfMammals <- GenerTime$Species; length(ListOfMammals) # 5426 - at least will use it to work with mammals only
MUT <- MUT %>% filter(Species %in% ListOfMammals)
rm(ListOfMammals)
nrow(MUT)  # 158771
length(unique(MUT$Species))  # 641

##### FILTER: normal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')
nrow(MUT)  # 158771

MUT %>% count(Species, sort = TRUE) -> cs
hist(cs$n, breaks = 50)
MUT2 <- MUT %>% 
    mutate(NormS = Subs %in% VecOfNormalSubstitutions) %>%
    group_by(Species) %>% 
    summarise(prop_of_normal = mean(NormS))
MUT2
hist(MUT2$prop_of_normal, breaks = 100)
cs2 <- filter(MUT2, prop_of_normal < .95) %>% inner_join(cs)
hist(cs2$n)
qplot(x = cs2$prop_of_normal, y = cs2$n)
knitr::kable(arrange(cs2, prop_of_normal), "pandoc")  # to see table of 81 species

MUT2 <- filter(MUT2, prop_of_normal >= .95)
hist(MUT2$prop_of_normal)

MUT <- MUT %>% filter(Subs %in% VecOfNormalSubstitutions, Species %in% MUT2$Species)
rm(cs, cs2, MUT2)
# MUT <- MUT[MUT$Subs %in% VecOfNormalSubstitutions, ]  # TO DO: check percent of normal substitutions across species. May be bias here!
nrow(MUT)  # if we remove just strange subst = 155185, 
# if we remove also species with a large proportion of strangers = 141601

##### FILTER: Synonymous Substitutions
MUT <- MUT %>% filter(AncestralAA == DescendantAA)
nrow(MUT)  # old = 132477 new = 121203

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT <- MUT %>% filter(AncestorCodon   %in% VecOfSynFourFoldDegenerateSites,
                      DescendantCodon %in% VecOfSynFourFoldDegenerateSites)
nrow(MUT)  # old = 62635 new = 57367

##### FILTER: Gene
table(MUT$Gene)
MUT <- MUT %>% filter(Gene == "CytB")
nrow(MUT)  # 36517

##### FILTER: species with many substitutions
MUT %>% count(Species, sort = TRUE) -> cs
hist(cs$n)
summary(cs$n)
# now quantile .15 = 13  and .25 = 21  substitutions
# so  quantile .15 = 433 and .25 = 383 species
quantile(cs$n, .15)
quantile(cs$n, .25)

# ListOfSpeciesWithManySubst <- AGG1[AGG1$x >= quantile(AGG1$x, 0.15), ]$Group.1  # why? it's only 11.65 substitutions.
# length(ListOfSpeciesWithManySubst) # for .15 = 486 for .25 = 436

# I take .15 because it's already more that's in previous case
ListOfSpeciesWithManySubst <- 
    filter(cs, n >= quantile(n, .15)) %>% 
    select(Species) %>% 
    as_vector()
MUT <- MUT <- filter(Species %in% ListOfSpeciesWithManySubst)
rm(ListOfSpeciesWithManySubst)
nrow(MUT)  # 36517

###### DERIVE MUTATIONAL SPECTRUM:
### NORMALIZATION of the 'NumberOfSynMutPerSpecies' by ancestral nucleotide count in the third position of four-fold synonymous substitutions:
NUC <- 
    read_table2('2_Derived/ATGC_counts_in_SYN_codons_wit_full_gene.txt') %>% 
    separate(Species, into = c("Species", "Gene"), sep = "[.]") %>% 
    gather(`CountA_Syn`, 
           `CountC_Syn`, 
           `CountG_Syn`, 
           `CountT_Syn`, 
           key = "Nucleotide", value = "Count_Syn") %>% 
    mutate(AncestralNuc = str_sub(Nucleotide, 6, 6)) %>%
    select(Species, Gene, AncestralNuc, Count_Syn)  # in case normalisation for synonymus it's better to shape long table

# MUT <- inner_join(MUT, NUC,  by = c("Species", "Gene"))  # compare CountA.x and CountA.y  - they should be identical.
# rm(NUC)
# nrow(MUT) # A bit less than before !!!! WHY???? for some species we don't have nucleotide count. Why???? Alya!!!!????
# It's identical length after join 36517

MUT %>% mutate(AncestralNuc = str_extract(Subs, "^.")) -> MUT
glimpse(MUT)

MUT <- inner_join(MUT, NUC,  by = c("Species", "Gene", "AncestralNuc"))
nrow(MUT)  # 36517

# Same normalisation for Syn_subst
mod_fun <- function(df){
    df %>% count(Count_Syn, sort = TRUE) -> cs
    syn <- unique(df$Count_Syn)
    cs$n / syn
}

nMUT <- 
    MUT %>% 
    group_by(AncestralNuc, Species) %>% 
    nest() %>% 
    mutate(Normalise = map_dbl(data, mod_fun)) %>% 
    group_by(Species) %>% 
    mutate(TotalMutRate = sum(Normalise)) %>% 
    unnest() %>% 
    group_by(Species, Subs) %>% 
    nest() %>% 
    mutate(MutTypeRate = map_dbl(data, mod_fun)) %>% 
    unnest() %>% 
    mutate(Fraction = MutTypeRate / TotalMutRate)

cMUT <- 
    complete(nMUT, Species, Subs) %>% 
    replace_na(list(MutTypeRate = 0, Fraction = 0))  # ???

sMUT <- 
    select(cMUT, Species, Subs, Fraction) %>% 
    distinct()
p_subs <- ggplot(sMUT) + geom_histogram(aes(Fraction)) + facet_wrap("Subs", ncol = 3)
p_subs2 <- ggplot(sMUT) + geom_freqpoly(aes(Fraction, colour = Subs, alpha = .3))

###### create matrix for PCA:
TABLE <- 
    sMUT %>% 
    spread(Subs, Fraction) %>% 
    rename("AC" = "A_C", "AG" = "A_G", "AT" = "A_T", 
           "CA" = "C_A", "CG" = "C_G", "CT" = "C_T", 
           "GA" = "G_A", "GC" = "G_C", "GT" = "G_T", 
           "TA" = "T_A", "TC" = "T_C", "TG" = "T_G")

### scale 12 substitutions. 
# If I don't scale PCA1 strongly depends on GA - probably this is just because GA is very common and thus fluctuations are very important by absolute values - not by relative...

### !!!!!!!!! DO IT LIKE BELOW OR LIKE THIS: df.pca <- prcomp(as.matrix(df.rw), center = TRUE, scale. = TRUE) 

#for (i in 2:13)
#{ # i =2
#  summary(MATRIX[, i])
#  MATRIX[, i] = as.numeric(scale(MATRIX[, i]))
#  summary(as.numeric(MATRIX[, i]))
#}


TABLE %>% as.data.frame() %>% column_to_rownames("Species") -> mtx
# for (i in seq_along(mtx)) {
#     hist(mtx[, i], main = colnames(mtx)[i], breaks = 50)
# }

# we should scale and center variables:
psych::describe(mtx)

# Check correspondence (contingency tables and chi-squre test)
# it's too large for nice visualisation,
# but at least you can use vcd::assoc() on subset
# cntgTable <- as.table(as.matrix(mtx))
# balloonplot(t(cntgTable), main = "Substitutions",
#             xlab = "", ylab = "", label = FALSE, 
#             show.margins = FALSE)
# mosaicplot(mtx, shade = TRUE, las = 2, main = "Substitutions")
# vcd::assoc(head(cntgTable, 10), shade = TRUE, las = 3)
chisq <- chisq.test(mtx, simulate.p.value = TRUE, B = 1e4)
tidy(chisq)  # not significant
#   statistic p.value parameter                     method
# 1  356.7633       1      5522 Pearson's Chi-squared test


# Contribution in percentage (%) of cell to Chi-square score:
contrib <- 100 * (chisq$residuals^2 / chisq$statistic)
for (i in seq_along(contrib)) {
    hist(contrib[, i], main = colnames(contrib)[i], breaks = 50)
}  # here we see group of strangers on GA hist with large conribution

outGA <- rownames(contrib[contrib[, "GA"] > .15, ])
mosaicplot(mtx[outGA, ], shade = TRUE, las = 2, main = "Substitutions")
cntgTableGA <- as.table(as.matrix(mtx[outGA, ]))
balloonplot(t(cntgTableGA), main = "Substitutions",
            xlab = "", ylab = "", label = FALSE,
            show.margins = FALSE)  # Probably here we see effect of replacement with 0 vs NA our completed matrix. When we normalise and had ratio, we can't just add 0 to matrix's empty column - it'll kill chi-square. But chi-square may be wrong and may be we need Yates' correction.
chi2 <- 356.7633
df <- (nrow(mtx) - 1) * (ncol(mtx) - 1)
pval <- pchisq(chi2, df = df, lower.tail = FALSE)
pval # NO ASSOCIATION BETWEEN SPECIFIC SPECIES AND SUBSTITUTIONS
rm(outGA, cntgTableGA, contrib, chisq, chi2, df, pval, mtx)

# if Y = n*p => Y = T(L)*F + E | L = K*n & F = K*p 
# (L - loadings and F - factors) 
# L should be orthogonal and F - orthonormal.
# So if we trust chi-squre test...
###### PCA 
TABLE <- calculatePcaReduction(data.use = TABLE, nPcs = 10)

p_subs3 <- ggplot(as.data.frame(subs.loadings), aes(x = PC1, y = PC2, color = Subs)) +
    geom_text(aes(label = Subs)) +
    xlab("PC1") +
    ylab("PC2") 
p_subs4 <- ggplot(as.data.frame(subs.loadings), aes(x = PC3, y = PC4, color = Subs)) +
    geom_text(aes(label = Subs)) +
    xlab("PC3") +
    ylab("PC4") 
pSubs <- plot_grid(p_subs, p_subs2, p_subs3, p_subs4, 
                   align = 'v', labels = LETTERS[1:4], 
                   scale = c(c(1, 1), c(1, 1), c(1, 1), c(1, 1)), 
                   hjust = -1, ncol = 2, axis = 't', rel_heights = c(1, 1.3))
if (!dir.exists("./4_Figures")) {
    dir.create("./4_Figures")
}
save_plot("./4_Figures/Substitutions_PCA.pdf", pSubs, base_height = 11, base_aspect_ratio = 1.1)

# QUESTION: IF WE USE MAMMALIAN SPECIES TO GET PCA (TO UNDERSTAND RULES OF PCA), TO SAVE ALL THESE LINEAR COMBINATIONS AND LATER TO USE THESE COMBINATIONS TO RECONSTRUCT PCs FOR OTHER MAMMALS, CHORDATA AND EVEN CANCERS!?
# FOR ME IT IS MEANINGFUL. ANY COMMENTS?

####################
#### FIND ECOLOGICAL INTERPRETATIONS OF PCs
####################

###### EXTREMES
TABLE %>% filter(PC1 < quantile(PC1, .05)) %>% select(Species)
TABLE %>% filter(PC1 > quantile(PC1, .95)) %>% select(Species)
TABLE %>% filter(PC2 < quantile(PC2, .05)) %>% select(Species)
TABLE %>% filter(PC2 > quantile(PC2, .95)) %>% select(Species)
TABLE %>% filter(PC3 < quantile(PC3, .05)) %>% select(Species)
TABLE %>% filter(PC3 > quantile(PC3, .95)) %>% select(Species)

###### GENERATION LENGTH
AnAge = read_tsv('1_Raw/anage_data.tsv')
Test <- inner_join(TABLE, GenerTime) %>% inner_join(AnAge) %>% type_convert()

Var1 <- list(Test %>% select(PC1:PC10) %>% colnames())
Var2 <- list(Test %>% select(AdultBodyMass_g:GenerationLength_d, 
                             Female.maturity..days.:Maximum.longevity..yrs., 
                             Metabolic.rate..W.:Temperature..K.) %>% colnames())

list_Tests <- cross2(Var1[[1]], Var2[[1]])

corrtest <- function(v) {
    fun <- function(x, y) cor.test(x, y, method = 'spearman')
    vars <- select(Test, v[[1]], v[[2]]) %>% type_convert()
    res <- with(Test, fun(vars[, 1], vars[, 2]))
    result <- data.frame(Var1 = v[[1]], Var2 = v[[2]])
    result <- cbind(result, tidy(res))
}

Tests_results <- list_Tests %>% map_df(~ corrtest(.))

knitr::kable(Tests_results, "pandoc")

###### HETEROTERMS
Hib = read_csv('1_Raw/tornew5.csv')
Hib$Species = gsub(" ", '_', Hib$Species)
HibOnlySpecies = Hib[Hib$Type == 'HIB', ]$Species; length(HibOnlySpecies)
DtOnlySpecies = Hib[Hib$Type == 'DT', ]$Species; length(DtOnlySpecies)
HibAndDtSpecies = Hib[Hib$Type == 'HIB' || Hib$Type == 'DT', ]$Species; length(HibAndDtSpecies)
boxplot(Test[Test$Species %in% HibOnlySpecies, ]$PC1, Test[!Test$Species %in% HibOnlySpecies, ]$PC1,  notch = TRUE,  names = c('HibOnlySpecies', 'Other'),  ylab = 'PC1'); # !!!! HibSpecies have a bit lower Pc1
boxplot(Test[Test$Species %in% DtOnlySpecies, ]$PC1, Test[!Test$Species %in% DtOnlySpecies, ]$PC1,  notch = TRUE,  names = c('HibOnlySpecies', 'Other'), ylab = 'PC1');   # nothing
# control for body mass? they are small! # other tests - decreased BMR as compared to body mass - compare strong outliers (cold and hot)...; animals, which live more than should according to bosy mass (naked mole rat)

###### MARSUPIALS, BATS... AGAIN NEED NORMAL TAXONOMY !!!!!!!!
Marsupials = c(AnAge[AnAge$Order == 'Diprotodontia' | AnAge$Order == 'Didelphimorphia' | AnAge$Order == 'Dasyuromorphia', ]$Species)
Placental = c(AnAge[AnAge$Order == 'Artiodactyla' | AnAge$Order == 'Cetacea' | AnAge$Order == 'Carnivora', ]$Species)
boxplot(MATRIX[MATRIX$Species %in% Marsupials, ]$Pca1, MATRIX[MATRIX$Species %in% Placental, ]$Pca1,  notch = TRUE); # !!!! HibSpecies have a bit lower Pc1

###### FIGURES:

MATRIX = merge(MATRIX, GenerTime)
MATRIX = MATRIX[order(MATRIX$GenerationLength_d), ]
MATRIX$Col = c(rep('green', 150), rep('gray', 187), rep('red', 150))
summary(PCA)
print(PCA)

pdf('4_FIGURES/PCA.pdf', width = 14, height = 14)
par(mfcol=c(2, 3))
summary(PCA)
#plot(PCA)
plot(MATRIX$Pca1, MATRIX$Pca2, col = MATRIX$Col)
plot(MATRIX$Pca2, MATRIX$Pca3, col = MATRIX$Col)
# plot(PCA$x[, 1], MATRIX$GenerationLength_d); cor.test(PCA$x[, 1], MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
plot(MATRIX$Pca2, log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2, MATRIX$GenerationLength_d, method = 'spearman')
plot(MATRIX$Pca3, log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2, MATRIX$GenerationLength_d, method = 'spearman') 
biplot(PCA, col = c('grey', 'black'), cex = 0.5)
biplot(PCA, choices=c(2, 3), col = c('grey', 'black'), cex = 0.5) #  biplot(princomp(USArrests), choices=c(1, 3))
dev.off()

####################
#### FIND DNA POLYMERAZE SIGNATURE (last PCs)!!!!!
####################
# LOGIC IS THE NEXT: first - third components are driven by temperature and ecology. 
# other components are not driven by ecology, physiology... sow we need to subtract effect of the first three PC and get naked signature of DNA polymeraze!!!
# how to do it carefully?!

pdf('4_FIGURES/DnaPolymerazeSignature.pdf', width = 14, height = 14)
biplot(PCA, choices=c(4, 5), col = c('grey', 'black'), cex = 0.5) #  biplot(princomp(USArrests), choices=c(1, 3))
dev.off()

PCA$x # PC's
PCA$sdev # the eigenvalues (res$sdev) giving information on the magnitude of each PC, 
PCA$rotation # and the loadings (res$rotation).

### IF I WANT TO DECREASE DIMENSIONALITY OF MY DATASET, I TRUNCATE IT, USING ONLY THE MOST IMPORTANT PCs:
pc.use <- 3 
trunc <- PCA$x[, 1:pc.use] %*% t(PCA$rotation[, 1:pc.use]) # trunc <- res$x[, 1:pc.use] %*% t(res$rotation[, 1:pc.use])

# add the center (and re-scale) back to data
PCA$scale
trunc <- scale(trunc, center = FALSE, scale=1/PCA$scale)
PCA$center
trunc <- scale(trunc, center = -1 * PCA$center, scale = FALSE)

### NOW - OPPOSITE EXERCISE: I WANT TO USE ONLY 4-12 PCs to reconstruct signature of gamma polymeraze:
start = 4
end = 12
gamma <- PCA$x[, start:end] %*% t(PCA$rotation[, start:end]) 
gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)

barplot(
    c(
        mean(gamma$AT), mean(gamma$AG), mean(gamma$AC),
        mean(gamma$TA), mean(gamma$TG), mean(gamma$TC),
        mean(gamma$CA), mean(gamma$CG), mean(gamma$CT),
        mean(gamma$GA), mean(gamma$GC), mean(gamma$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT')
    )

###### doesn't work!!! scaling and centering ...
### try without scaling:
par(mfrow = c(2, 1))
PcaNoScale = prcomp(matrix, center = FALSE, scale. = FALSE) # !!! or scale above!! discuss - it seems identical things
start = 1
end = 3
gamma <-
    PcaNoScale$x[, start:end] %*% t(PcaNoScale$rotation[, start:end])
#gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
#gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)
barplot(
    c(
        mean(gamma$AT), mean(gamma$AG), mean(gamma$AC),
        mean(gamma$TA), mean(gamma$TG), mean(gamma$TC),
        mean(gamma$CA), mean(gamma$CG), mean(gamma$CT),
        mean(gamma$GA), mean(gamma$GC), mean(gamma$GT)),
    names = c(
        'AT', 'AG', 'AC',
        'TA', 'TG', 'TC',
        'CA', 'CG', 'CT',
        'GA', 'GC', 'GT'),
    main = 'EXO'
)

start = 4
end = 10
gamma <-
    PcaNoScale$x[, start:end] %*% t(PcaNoScale$rotation[, start:end])
#gamma <- scale(gamma, center = FALSE , scale=1/PCA$scale)
#gamma <- scale(gamma, center = -1 * PCA$center, scale=FALSE)
gamma = as.data.frame(gamma)
barplot(
    c(
        mean(gamma$AT), mean(gamma$AG), mean(gamma$AC),
        mean(gamma$TA), mean(gamma$TG), mean(gamma$TC),
        mean(gamma$CA), mean(gamma$CG), mean(gamma$CT),
        mean(gamma$GA), mean(gamma$GC), mean(gamma$GT)),
    names = c(
        'AT', 'AG', 'AC',
        'TA', 'TG', 'TC',
        'CA', 'CG', 'CT',
        'GA', 'GC', 'GT'),
    main = 'ENDO'
)

MatrixNoScale = data.frame(matrix)
MatrixNoScale$Species = row.names(MatrixNoScale)
PcaNoScaleDataFrame = data.frame(PcaNoScale$x)
PcaNoScaleDataFrame$Species = row.names(PcaNoScaleDataFrame)
NoScale = merge(MatrixNoScale, PcaNoScaleDataFrame)

#### NAIVE CHECKS BELOW - don't work!!! don't understand!!
A = lm(AT ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$ATRes = residuals(A)

A = lm(AG ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$AGRes = residuals(A)

A = lm(AC ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$ACRes = residuals(A)

A = lm(TA ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$TARes = residuals(A)

A = lm(TG ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$TGRes = residuals(A)

A = lm(TC ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$TCRes = residuals(A)

A = lm(CG ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$CGRes = residuals(A)

A = lm(CT ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$CTRes = residuals(A)

A = lm(CA ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$CARes = residuals(A)

A = lm(GT ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$GTRes = residuals(A)

A = lm(GC ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$GCRes = residuals(A)

A = lm(GA ~ PC1 + PC2 + PC3, data = NoScale)
summary(A)
NoScale$GARes = residuals(A)

barplot(
    c(
        mean(NoScale$ATRes), mean(NoScale$AGRes), mean(NoScale$ACRes),
        mean(NoScale$TARes), mean(NoScale$TGRes), mean(NoScale$TCRes),
        mean(NoScale$CARes), mean(NoScale$CGRes), mean(NoScale$CTRes),
        mean(NoScale$GARes), mean(NoScale$GCRes), mean(NoScale$GTRes)),
    names = c(
        'AT', 'AG', 'AC',
        'TA', 'TG', 'TC',
        'CA', 'CG', 'CT',
        'GA', 'GC', 'GT')
)
