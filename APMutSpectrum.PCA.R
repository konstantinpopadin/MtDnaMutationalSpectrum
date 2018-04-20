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
rm(outGA, cntgTableGA, contrib, chisq, chi2, df, pval)

# if Y = n*p => Y = T(L)*F + E | L = K*n & F = K*p 
# (L - loadings and F - factors) 
# L should be orthogonal and F - orthonormal.
# So if we trust chi-squre test...
###### PCA 
TABLE <- calculatePcaReduction(data.use = TABLE, nPcs = 10, center = TRUE, weight.by.var = FALSE, rev.pca = FALSE, seed.use = 42)

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
cl_results <- Tests_results %>% filter(p.value < .001)
knitr::kable(cl_results, "pandoc", caption = "TABLE OF SIGNIFICANT CORRELATIONS BETWEEN PCS AND ECOLOGY")
if (!dir.exists("./3_Results")) {
    dir.create("./3_Results")
}
pdf("./3_Results/Significant_PCs_Correlations.pdf")
gridExtra::grid.table(cl_results)
dev.off()

# Table: TABLE OF SIGNIFICANT CORRELATIONS BETWEEN PCS AND ECOLOGY
# 
# Var1   Var2                                 estimate   statistic     p.value  method                            alternative 
# -----  ---------------------------------  ----------  ----------  ----------  --------------------------------  ------------
# PC2    AdultBodyMass_g                     0.2822715   1960213.3   0.0000049  Spearman's rank correlation rho   two.sided   
# PC2    Max_longevity_d                     0.2486489   1709534.2   0.0001023  Spearman's rank correlation rho   two.sided   
# PC2    Rspan_d                             0.2629621    842533.9   0.0002469  Spearman's rank correlation rho   two.sided   
# PC2    AFR_d                               0.2854507    816826.3   0.0000655  Spearman's rank correlation rho   two.sided   
# PC2    Calculated_GL_d                     0.2878064    801346.0   0.0000636  Spearman's rank correlation rho   two.sided   
# PC2    GenerationLength_d                  0.2539992   2037429.0   0.0000421  Spearman's rank correlation rho   two.sided   
# PC2    Female.maturity..days.              0.2608739    831649.9   0.0002883  Spearman's rank correlation rho   two.sided   
# PC2    Gestation.Incubation..days.         0.2519047   1291713.0   0.0001708  Spearman's rank correlation rho   two.sided   
# PC2    Weaning..days.                      0.2751545    777356.9   0.0001441  Spearman's rank correlation rho   two.sided   
# PC2    Inter.litter.Interbirth.interval    0.2901497    440547.3   0.0002501  Spearman's rank correlation rho   two.sided   
# PC2    Birth.weight..g.                    0.2657801    949859.6   0.0001540  Spearman's rank correlation rho   two.sided   
# PC2    Adult.weight..g.                    0.2768955   1905742.0   0.0000085  Spearman's rank correlation rho   two.sided   
# PC2    Maximum.longevity..yrs.             0.2569373   1265447.7   0.0001294  Spearman's rank correlation rho   two.sided   

###### HETEROTERMS
Hib = read_csv('1_Raw/tornew5.csv')
Hib$Species = gsub(" ", '_', Hib$Species)
HibOnlySpecies = Hib[Hib$Type == 'HIB', ]$Species; length(HibOnlySpecies)
DtOnlySpecies = Hib[Hib$Type == 'DT', ]$Species; length(DtOnlySpecies)
Tab1 <- Hib %>% select(Species, Taxon, Type:Boutmean, lat:prec) %>% right_join(Test)
ggplot(Tab1, aes(Type, PC1)) + geom_boxplot()
# boxplot(Test[Test$Species %in% HibOnlySpecies, ]$PC1, Test[!Test$Species %in% HibOnlySpecies, ]$PC1,  notch = TRUE,  names = c('HibOnlySpecies', 'Other'),  ylab = 'PC1'); # !!!! HibSpecies have a bit lower load of PC1
# boxplot(Test[Test$Species %in% DtOnlySpecies, ]$PC1, Test[!Test$Species %in% DtOnlySpecies, ]$PC1,  notch = TRUE,  names = c('DtOnlySpecies', 'Other'), ylab = 'PC1');   # nothing
# boxplot(Test[Test$Species %in% HibOnlySpecies, ]$PC1, Test[Test$Species %in% DtOnlySpecies, ]$PC1,  notch = TRUE,  names = c('HibOnlySpecies', 'DtOnlySpecies'),  ylab = 'PC1')  # DT have more impute from PC1
# control for body mass? they are small! # other tests - decreased BMR as compared to body mass - compare strong outliers (cold and hot)...; animals, which live more than should according to bosy mass (naked mole rat)

###### MARSUPIALS, BATS... AGAIN NEED NORMAL TAXONOMY !!!!!!!!
Tax <- select(Hib, Taxon, Order) %>% distinct() %>% filter(Taxon != "Mar") %>% mutate(Order = str_to_title(Order))
Tab1 <- left_join(Test, Tax)
Tab1 %>% filter(is.na(Taxon)) %>% janitor::tabyl(Order)

Marsupials <- c("Didelphimorphia")
Placental <- c("Artiodactyla", "Cetacea", "Cingulata", "Erinaceomorpha", "Lagomorpha", "Perissodactyla", "Pilosa", "Proboscidea", "Scandentia", "Soricomorpha")

Tab1 %>% mutate(Taxon = ifelse(is.na(Taxon), ifelse(Order %in% Placental, "Plac", "Mars"), Taxon)) %>% janitor::tabyl(Order, Taxon)
Tab1 <- Tab1 %>% mutate(Taxon = ifelse(is.na(Taxon), ifelse(Order %in% Placental, "Plac", "Mars"), Taxon))
ggplot(Tab1, aes(Taxon, PC1)) + geom_boxplot()
# boxplot(Test[Test$Species %in% Marsupials, ]$PC1, Test[Test$Species %in% Placental, ]$PC1,  notch = TRUE, names = c('Marsupials', 'Placental'),  ylab = 'PC1'); # !!!! Marsupials have a bit lower load of PC1

###### FIGURES:
Tab2 <- Tab1 %>% arrange(GenerationLength_d) %>% group_by(Taxon) %>% mutate(GL_grops = ntile(GenerationLength_d, 5))

write_csv(Tab2, "./3_Results/PC2_Grouped_by_GL_taxon.csv")

p_SpEmb1 <- ggplot(data = Tab2, aes(x = PC1, y = PC2, colour = GL_grops, alpha = .4)) + geom_point()
p_SpEmb2 <- ggplot(data = Tab2, aes(x = PC2, y = PC3, colour = GL_grops, alpha = .4)) + geom_point()

# plot(PCA$x[, 1], MATRIX$GenerationLength_d); cor.test(PCA$x[, 1], MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
# ggplot(data = Tab2, aes(x = PC2, y = GenerationLength_d, colour = GL_grops, alpha = .4)) + geom_point() + geom_smooth(data = Tab2, mapping = aes(x = PC2, y = GenerationLength_d, linetype = Taxon)) + scale_y_log10()
p_SpGL1 <- ggplot(data = Tab2, aes(x = PC2, y = GenerationLength_d, colour = GL_grops, alpha = .4)) + geom_point() + scale_y_log10(); cor.test(Tab2$PC2, Tab2$GenerationLength_d, method = 'spearman')
p_SpGL2 <- ggplot(data = Tab2, aes(x = PC3, y = GenerationLength_d, colour = GL_grops, alpha = .4)) + geom_point() + scale_y_log10(); cor.test(Tab2$PC3, Tab2$GenerationLength_d, method = 'spearman') 

pEco <- plot_grid(p_SpEmb1, p_SpEmb2, p_SpGL1, p_SpGL2, p_subs3, p_subs4, ncol = 2, nrow = 3)
save_plot('4_Figures/Ecology_PCA.pdf', pEco, base_height = 14)
####################
#### FIND DNA POLYMERAZE SIGNATURE (last PCs)!!!!!
####################
# LOGIC IS THE NEXT: first - third components are driven by temperature and ecology. 
# other components are not driven by ecology, physiology... sow we need to subtract effect of the first three PC and get naked signature of DNA polymeraze!!!
# how to do it carefully?!

p_SpEmb3 <- ggplot(data = Tab2, aes(x = PC4, y = PC5, colour = GL_grops, alpha = .4)) + geom_point()
save_plot('4_Figures/DnaPolymerazeSignature_PCA.pdf', p_SpEmb3, base_height = 7)

# Original matrix
barplot(
    c(
        mean(mtx$AT), mean(mtx$AG), mean(mtx$AC),
        mean(mtx$TA), mean(mtx$TG), mean(mtx$TC),
        mean(mtx$CA), mean(mtx$CG), mean(mtx$CT),
        mean(mtx$GA), mean(mtx$GC), mean(mtx$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT'),
    main = "Original matrix"
)

# Reproduce original matrix of substitution
species.embeddings <- species.embeddings %>% as.data.frame() %>% column_to_rownames("Species") %>% as.matrix()
subs.loadings <- subs.loadings %>% as.data.frame() %>% column_to_rownames("Subs") %>% as.matrix()
repr <- t(species.embeddings %*% t(subs.loadings))
repr1 <- scale(repr, center = -1 * pcs$center, scale = FALSE)
repr1 <- as.data.frame(t(repr1))
barplot(
    c(
        mean(repr1$AT), mean(repr1$AG), mean(repr1$AC),
        mean(repr1$TA), mean(repr1$TG), mean(repr1$TC),
        mean(repr1$CA), mean(repr1$CG), mean(repr1$CT),
        mean(repr1$GA), mean(repr1$GC), mean(repr1$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT'),
    main = "Reproduce original matrix"
)

### IF I WANT TO DECREASE DIMENSIONALITY OF MY DATASET, I TRUNCATE IT, USING ONLY THE MOST IMPORTANT PCs:
pc.use <- 3
trunc <- t((as.matrix(species.embeddings[, 1:pc.use])) %*% t(as.matrix(subs.loadings[, 1:pc.use])))

# add the center (and re-scale) back to data
trunc <- scale(trunc, center = -1 * pcs$center, scale = FALSE)
trunc <- as.data.frame(t(trunc))

barplot(
    c(
        mean(trunc$AT), mean(trunc$AG), mean(trunc$AC),
        mean(trunc$TA), mean(trunc$TG), mean(trunc$TC),
        mean(trunc$CA), mean(trunc$CG), mean(trunc$CT),
        mean(trunc$GA), mean(trunc$GC), mean(trunc$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT'),
    main = "Matrix produced by PC1:PC3"
)
trunc <- rownames_to_column(trunc, "Species")
write_csv(trunc, "./3_Results/PCs1-3_Truncated_tab.csv")


### NOW - OPPOSITE EXERCISE: I WANT TO USE ONLY 4-10 PCs to reconstruct signature of gamma polymeraze:
start = 4
end = 10
gamma <- t(species.embeddings[, start:end] %*% t(subs.loadings[, start:end]))
gamma <- scale(gamma, center = -1 * pcs$center, scale = FALSE)
gamma <- as.data.frame(t(gamma))

barplot(
    c(
        mean(gamma$AT), mean(gamma$AG), mean(gamma$AC),
        mean(gamma$TA), mean(gamma$TG), mean(gamma$TC),
        mean(gamma$CA), mean(gamma$CG), mean(gamma$CT),
        mean(gamma$GA), mean(gamma$GC), mean(gamma$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT'),
    main = "Matrix produced by PC4:PC10"
    )
gamma <- rownames_to_column(gamma, "Species")
write_csv(gamma, "./3_Results/DnaPolymerazeSignature_tab.csv")

#### NAIVE CHECKS BELOW - work well
#### Here you shoud use predict.
NaivChk <- 
    TABLE %>% 
    select(AC:TG) %>% 
    map(~ lm(.x ~ PC1 + PC2 + PC3, data = TABLE)) %>% 
    map_dfc(predict)

barplot(
    c(
        mean(NaivChk$AT), mean(NaivChk$AG), mean(NaivChk$AC),
        mean(NaivChk$TA), mean(NaivChk$TG), mean(NaivChk$TC),
        mean(NaivChk$CA), mean(NaivChk$CG), mean(NaivChk$CT),
        mean(NaivChk$GA), mean(NaivChk$GC), mean(NaivChk$GT)),
    names = c('AT', 'AG', 'AC',
              'TA', 'TG', 'TC',
              'CA', 'CG', 'CT',
              'GA', 'GC', 'GT'),
    main = "NAIVE CHECKS"
)
