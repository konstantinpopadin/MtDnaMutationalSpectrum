########################
################ ATGC FRACTIONS ALONG GENOMES
########################

rm(list=ls(all=TRUE))

### set working directory as a directory where R script was open
wd <- getwd()
setwd(wd)

############ Syn mut
SynNuc = read.table('Body/2_Derived/GcAtSkewNucContBig.csv', header = TRUE)
SynNuc$FrA = SynNuc$A / SynNuc$SitesNumber
SynNuc$FrT = SynNuc$T / SynNuc$SitesNumber
SynNuc$FrG = SynNuc$G / SynNuc$SitesNumber
SynNuc$FrC = SynNuc$C / SynNuc$SitesNumber

VecOfTaxa = unique(SynNuc$TAXON)
SynNucAll = SynNuc

pdf("./Body/4_Figures/APMutSpectrum.SynSites.AtgcAlongGenomes.R.01.pdf", height = 10, width = 15)
par(mfrow = c(2,1))

for (taxa in 1:length(VecOfTaxa))
{ # taxa = 1
  TAX = as.character(VecOfTaxa[taxa])
  SynNuc = SynNucAll
    
  Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2') # ATP6 and ND4 
  Timing = seq(1:13)
  NewData = data.frame(Gene,Timing)
  SynNuc = merge(SynNuc,NewData)
  
  VecOfSpecies  = as.character(unique(SynNuc[SynNuc$TAXON == TAX,]$Species))
  
  plot(NA, xlim=c(1,13), ylim=c(0,0.7), xlab='', ylab="Nucleotide Fractions", main = TAX, xaxt="n")
  axis(side = 1, at=c(1:13), labels=c(Gene), las = 2) 
    for (i in 1:length(VecOfSpecies))
    { # i = 1
      Temp = SynNuc[SynNuc$Species == VecOfSpecies[i],]
      Temp = Temp[order(Temp$Timing),]
      if (nrow(Temp) == 13) # 10
      {
        for (count in 1:(nrow(Temp)-1))
        {
          segments(count, Temp$FrA[count], count+1, Temp$FrA[count+1], col = rgb(1,0.1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrT[count], count+1, Temp$FrT[count+1], col = rgb(0.1,1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrG[count], count+1, Temp$FrG[count+1], col = rgb(0.1,1,1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrC[count], count+1, Temp$FrC[count+1], col = rgb(0.1,0.1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
        }
      }
    }
      legend(13,0.6,legend=c('A','T','G','C'), col = c(rgb(1,0.1,0.1,1),rgb(0.1,1,0.1,1),rgb(0.1,1,1,1),rgb(0.1,0.1,0.1,1)), pch = 16, horiz = FALSE)
  
  #### the same but without problematic genes:
  Gene = c('COX1','COX2','ATP8','ND3','ND4L','ND5','CytB', 'ND1','ND2') # ATP6 and ND4 
  Timing = seq(1:9)
  NewData = data.frame(Gene,Timing)
  SynNuc = SynNuc[,-15]
  SynNuc = merge(SynNuc,NewData)
  
  VecOfSpecies  = as.character(unique(SynNuc[SynNuc$TAXON == TAX,]$Species))
  
  plot(NA, xlim=c(1,9), ylim=c(0,0.7), xlab='', ylab="Nucleotide Fractions", main = TAX, xaxt="n")
  axis(side = 1, at=c(1:9), labels=c(Gene), las = 2) 
  for (i in 1:length(VecOfSpecies))
    { # i = 1
      Temp = SynNuc[SynNuc$Species == VecOfSpecies[i],]
      Temp = Temp[order(Temp$Timing),]
      if (nrow(Temp) == 9) # 10
      {
        for (count in 1:(nrow(Temp)-1))
        {
          segments(count, Temp$FrA[count], count+1, Temp$FrA[count+1], col = rgb(1,0.1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrT[count], count+1, Temp$FrT[count+1], col = rgb(0.1,1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrG[count], count+1, Temp$FrG[count+1], col = rgb(0.1,1,1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
          segments(count, Temp$FrC[count], count+1, Temp$FrC[count+1], col = rgb(0.1,0.1,0.1,0.1), lwd = 1) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
        }
      }
    }
  legend(9,0.6,legend=c('A','T','G','C'), col = c(rgb(1,0.1,0.1,1),rgb(0.1,1,0.1,1),rgb(0.1,1,1,1),rgb(0.1,0.1,0.1,1)), pch = 16, horiz = FALSE)
}

dev.off()

#### 
