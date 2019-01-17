################## 
# Author: Vinicius Henrique da Silva
# Script description: PCA and Fst analysis to identify a large inversion in chromosome 1A in great tits
# Date: January - 2019
##################


################### PCA ####################
### Load the gds file
(genofile <- SNPRelate::snpgdsOpen("GT201420152016ChrNumeric.gds"))

chrA <- as.data.frame(cbind(read.gdsn(index.gdsn(genofile, "snp.id")),
                            read.gdsn(index.gdsn(genofile, "snp.chromosome")),
                            read.gdsn(index.gdsn(genofile, "snp.position"))
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
p <- ggplot2::ggplot(tab, aes(EV1, EV2))
p + ggplot2::geom_point(aes(size = 4, alpha = 0.5, shape="o"))
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

snps_noinv <- snps[,Hmisc::`%nin%`(colnames(snps), invsam$sample.id), drop=FALSE] ## Select samples with inversion
snps_noinv <- snps_noinv[rownames(snps_noinv) %in% chrx$probeset_id,,drop=FALSE] ## Select samples with inversion

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

### Prepare the SNP file to infer FST values
SnpsAll <- cbind(snps_inv, snps_noinv)
SnpsAll <- SnpsAll+1
SnpsAll1 <- SnpsAll
SnpsAll <- new("SnpMatrix", t(SnpsAll))

### Infer FST for each SNP marker 
f1 <- Fst(SnpsAll, as.factor(allSam$stratum), pairwise=FALSE)

#####Transform negative values into zero
f1$Fst[f1$Fst<0] <- 0

FstRes <- cbind(as.data.frame(rownames(SnpsAll1)), as.data.frame(f1$Fst), as.data.frame(f1$weight))
colnames(FstRes) <- c("probeset_id", "Fst", "weight")
chrx <- chrx[with(chrx, order(Position)), ]

FstRes <- merge(chrx, FstRes, by="probeset_id", sort=FALSE)

#### Plot Fst
pdf("FST_chr1A_inverted_samplesXnoninverted.pdf")
p <- ggplot(FstRes, aes(Position, Fst))
p + geom_point() + theme_bw() + geom_smooth(method = "lm")
dev.off()
