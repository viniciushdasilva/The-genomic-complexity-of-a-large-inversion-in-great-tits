################## 
# Author: Vinicius Henrique da Silva
# Script description: PCA and Fst analysis to identify a large inversion in chromosome 1A in great tits
# Date: January - 2019
##################


################### PCA ####################
path.gen <- "path/to/gds"

### Load the gds file
(genofile <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds"), readonly=FALSE))

chrA <- as.data.frame(cbind(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
                            gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
                            gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))
))

chrA$V3 <- as.integer(as.character(chrA$V3))
colnames(chrA) <- c("probeset_id", "Chromosome", "Position")

### Select chr 1A
chrx <- subset(chrA, Chromosome == "16") ## Select chromosome 1A - numeric index==16

### Run PCA with all samples
pca <- SNPRelate::snpgdsPCA(genofile, snp.id=chrx$probeset_id, autosome.only=F)

# variance proportion (%)
pc.percent <- pca$varprop*100
pc.percent <- (round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

## Plot PCA
pdf("PCA_analysis_chr1A.pdf")
p <- ggplot2::ggplot(tab, ggplot2::aes(EV1, EV2))
p + ggplot2::geom_point(ggplot2::aes(size = 4, alpha = 0.5, shape="o"))
dev.off()
#################### END ########################


################### Fst #########################
# Select birds with inversion 
invsam <- subset(tab, EV1 >= 0.04) ####### Threshold - Based on the PCA figure 

###################### Formatation of the genotypes

snps <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "genotype"))
samplesA <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
ProNamA <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))
colnames(snps) <- samplesA 
rownames(snps) <- ProNamA
snps <- as.data.frame(snps)

snps_inv <- snps[,colnames(snps) %in% invsam$sample.id, drop=FALSE] ## Select samples with inversion
snps_inv <- snps_inv[rownames(snps_inv) %in% chrx$probeset_id,,drop=FALSE]

snps_noinv <- snps[,Hmisc::`%nin%`(colnames(snps), invsam$sample.id), drop=FALSE] ## Select samples without inversion
snps_noinv <- snps_noinv[rownames(snps_noinv) %in% chrx$probeset_id,,drop=FALSE] 

############# Transform a normal matrix in a 'SnpMatrix' to use it in 'snpStats' package
Invsamp <- as.data.frame(colnames(snps_inv))
colnames(Invsamp) <- "sam"
Invsamp$stratum <- "Inversion"
NoInvsamp <- as.data.frame(colnames(snps_noinv))
colnames(NoInvsamp) <- "sam"
NoInvsamp$stratum <- "Normal"

allSam <- rbind(Invsamp, NoInvsamp)
allSam$stratum
rownames(allSam) <- allSam$sam

save(allSam, file=file.path(path.gen, "NormInveSamples.rda"))

## Attach the classification to the GDS
gdsfmt::add.gdsn(genofile, "inversion.classification", val=allSam, replace=TRUE)
gdsfmt::closefn.gds(genofile)

### Prepare the SNP file to infer FST values
SnpsAll <- cbind(snps_inv, snps_noinv)
SnpsAll <- SnpsAll+1
SnpsAll1 <- SnpsAll
library(snpStats)
SnpsAll <- new("SnpMatrix", t(SnpsAll))

### Infer FST for each SNP marker 
f1 <- Fst(SnpsAll, as.factor(allSam$stratum), pairwise=FALSE)

##### Transform negative values into zero
f1$Fst[f1$Fst<0] <- 0

FstRes <- cbind(as.data.frame(rownames(SnpsAll1)), as.data.frame(f1$Fst), as.data.frame(f1$weight))
colnames(FstRes) <- c("probeset_id", "Fst", "weight")
chrx <- chrx[with(chrx, order(Position)), ]

FstRes <- merge(chrx, FstRes, by="probeset_id", sort=FALSE)

#### Plot Fst
pdf("FST_chr1A_inverted_samplesXnoninverted.pdf")
p <- ggplot2::ggplot(FstRes, ggplot2::aes(Position, Fst))
p + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::geom_smooth(method = "lm")
dev.off()


### LD Analysis
#### Non-phased chromosomes for samples WITH the inversion ########################################################

#### Parameters
MAFv <- 0.40
LDT <- 0.05
SK <- "InvSamples"
TI <- "inv-norm"
#MET <- "dprime" ## Rerun with dprime if needed
MET <- "r"
#LDPlotMethod <- "D'"  ## Rerun with dprime if needed
LDPlotMethod <- "r"
FigNam <- NULL
###### end parameters

allSamInv <- subset(allSam, stratum=="Inversion")
samplesX <- as.character(allSamInv$sam)

(gdsAN <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds"), readonly=FALSE))
map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.id")), 
                           chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.chromosome")), 
                           position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.position")), 
                           snp.rs.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.rs.id"))))
mapchrX <- subset(map, chr == 16) ## Select chromosome 1A

gdsfmt::closefn.gds(gdsAN)

GWASTools::gdsSubset(file.path(path.gen, "GT201420152016ChrNumeric.gds"), 
          file.path(path.gen,"LDInversionSamples.gds"), 
          snp.include=mapchrX$snp.id, sample.include=samplesX)
(genofile <- SNPRelate::snpgdsOpen(file.path(path.gen, "LDInversionSamples.gds")))

###################### LD analysis

set.seed(1000)
snpset <- SNPRelate::snpgdsLDpruning(genofile, method=MET, remove.monosnp = FALSE,
                          ld.threshold = LDT,
                          maf=MAFv) 

print(paste0("Number of SNP markers == ", length(snpset[[1]])))

LD <- SNPRelate::snpgdsLDMat(genofile, 
                  snp.id=unlist(snpset), 
                  method=c(MET), 
                  slide=0L)

mapchrXSe <- mapchrX[mapchrX$snp.id %in% snpset[[1]],]

pdf(file.path(path.gen, paste0("LD",SK,"Chr1A.maf",MAFv, ".ld", LDT, MET, ".pdf")))
LDheatmap::LDheatmap(LD$LD, as.numeric(as.character(mapchrXSe$position)), color = (heat.colors(20)),
          title=paste0("Pairwise LD ",TI), LDmeasure=LDPlotMethod)
graphics.off()
gdsfmt::closefn.gds(genofile)

######################################################## END ####################################################

#### Non-phased chromosomes for ALL samples - New SNPs ########################################################

#### Parameters 
SK <- "AllSamples"
TI <- "inv-norm + norm-norm"
###### end parameters

(gdsAN <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds"), readonly=FALSE))
map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.id")), 
                           chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.chromosome")), 
                           position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.position")),
                           snp.rs.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.rs.id"))))
mapchrX <- subset(map, chr == 16) ## Select chromosome 1A

gdsfmt::closefn.gds(gdsAN)

GWASTools::gdsSubset(file.path(path.gen, "GT201420152016ChrNumeric.gds"),
                     file.path(path.gen, "LDInversionSamples.gds"),
                     snp.include=mapchrX$snp.id, sample.include=samplesX)
(genofile <- snpgdsOpen(file.path(path.gen,"LDInversionSamples.gds")))

###################### LD analysis 

set.seed(1000)
snpset <- SNPRelate::snpgdsLDpruning(genofile, method=MET, remove.monosnp = FALSE,
                          ld.threshold = LDT,
                          maf=MAFv) 

LD <- snpgdsLDMat(genofile, 
                  snp.id=unlist(snpset), 
                  method=c(MET), 
                  slide=0L) 

mapchrXSe <- mapchrX[mapchrX$snp.id %in% snpset[[1]],]

pdf(file.path(path.gen, paste0("LD",SK,"Chr1A.maf",MAFv, ".ld", LDT, ".pdf")))
LDheatmap::LDheatmap(LD$LD, as.numeric(as.character(mapchrXSe$position)),
          color = (heat.colors(20)), title=paste0("Pairwise LD ",TI),
          LDmeasure=LDPlotMethod)
graphics.off()
gdsfmt::closefn.gds(genofile)
######################################################## END ####################################################

#### Non-phased chromosomes for samples WITHOUT the inversion ###################################################

#### Parameters 
SK <- "NormSamples"
TI <- "norm-norm"
#### end Parameters

allSamInv <- subset(allSam, stratum=="Normal")
samplesX <- as.character(allSamInv$sam)

(gdsAN <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds"), readonly=FALSE))
map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.id")),
                           chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.chromosome")), 
                           position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.position")),
                           snp.rs.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gdsAN, "snp.rs.id"))))
mapchrX <- subset(map, chr == 16) ## Select chromosome 1A

gdsfmt::closefn.gds(gdsAN)

GWASTools::gdsSubset(file.path(path.gen, "GT201420152016ChrNumeric.gds"),
                     file.path(path.gen, "LDInversionSamples.gds"), 
                     snp.include=mapchrX$snp.id, sample.include=samplesX)
(genofile <- snpgdsOpen(file.path(path.gen, "LDInversionSamples.gds")))

###################### LD analysis

set.seed(1000)
LD <- SNPRelate::snpgdsLDMat(genofile, 
                  snp.id=unlist(snpset), 
                  method=c("r"),
                  slide=0L)

pdf(file.path(path.gen, paste0("LD",SK,"Chr1A.maf",MAFv, ".ld", LDT, ".pdf")))
LDheatmap::LDheatmap(LD$LD, as.numeric(as.character(mapchrXSe$position)),
          color = (heat.colors(20)), title=paste0("Pairwise LD ",TI),
          LDmeasure=LDPlotMethod)
graphics.off()
gdsfmt::closefn.gds(genofile)

######################################################## END ####################################################

                          
