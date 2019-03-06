################## 
# Author: Vinicius Henrique da Silva
# Script description: Calculate inversion heterozigosity to define informative SNPs. 
# Create a heat map with informative SNP genotypes at the center of the inversion.
# Date: March - 2019
##################

## !WORKS IN UNIX ONLY! ###
################# Parameters
path.gen <- "path/to/gds"
Ncor <- 1 ## Number of cores to be used

library(Hmisc)
library(doParallel)

# Open the GDS file
(genofile <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds")))

load(file=file.path("NormInveSamples.rda"))
allSamNorm <- subset(allSam, stratum=="Normal")
allSamInv <- subset(allSam, stratum=="Inversion")

map <- as.data.frame(cbind(snp.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
                           chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")), 
                           position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")), 
                           snp.rs.id = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id"))))

mapX <- subset(map, chr == 16) ## 16 is the numeric code of chr1A

g <- SNPRelate::snpgdsGetGeno(genofile, sample.id=allSamNorm$sam, 
                              snp.id = mapX$snp.id) ### Subset normal samples in chr1A

doMC::registerDoMC(Ncor) # Define number of cores to run in parallel 

Hete.All  <- foreach::foreach(lo=1:nrow(g), .errorhandling=c('pass')) %dopar% {
  print(lo)
  if(1 %in% names(table(g[lo,])[names(table(g[lo,]))])){
    t<-table(g[lo,])[names(table(g[lo,]))==1][[1]]
  }
  if(1 %nin% names(table(g[lo,])[names(table(g[lo,]))])){
    t <- 0  
  }
  t
}

HetNorm <- as.numeric(unlist(Hete.All))/length(allSamNorm$sam)

g <- SNPRelate::snpgdsGetGeno(genofile, sample.id=allSamInv$sam, snp.id = mapX$snp.id) ### Subset inversion samples in chr1A

Hete.All  <- foreach::foreach(lo=1:nrow(g)) %dopar% {
  print(lo)
  if(1%in%names(table(g[lo,])[names(table(g[lo,]))])){
    t<-table(g[lo,])[names(table(g[lo,]))==1][[1]]
  }
  if(1%nin%names(table(g[lo,])[names(table(g[lo,]))])){
    t <- 0  
  }
  t
}

HetInv <- as.numeric(unlist(Hete.All))/length(allSamInv$sam)

HetNorm <- cbind(mapX,HetNorm)
HetNorm$Pop <- "Normal"
HetNorm$position <- as.numeric(as.character(HetNorm$position))/1e+6
HetInv <- cbind(mapX,HetInv)
HetInv$Pop <- "Inversion"
HetInv$position <- as.numeric(as.character(HetInv$position))/1e+6
colnames(HetInv) <- colnames(HetNorm)
HetAll <- rbind(HetNorm, HetInv)

## Number of informative SNPs heterozygosity > 60%
informativeSNPs  <- subset(HetInv, HetNorm>0.6)
nrow(informativeSNPs) ## Total of 6624
save(HetInv, HetNorm, informativeSNPs, file=file.path(path.gen, "heterozygosity_inversion.rda"))
gdsfmt::closefn.gds(genofile)

## Count ''AA'' genotypes in the center of the inversion
load(file=file.path(path.gen, "heterozygosity_inversion.rda"))

library(GWASTools)
library(LDheatmap)
library(SNPRelate)
library(doMC)
library(data.table)
library(gridExtra)

(gdsAN <- snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds"), allow.fork=TRUE))
map <- as.data.frame(cbind(snp.id = read.gdsn(index.gdsn(gdsAN, "snp.id")), chr = read.gdsn(index.gdsn(gdsAN, "snp.chromosome")), 
                           position = read.gdsn(index.gdsn(gdsAN, "snp.position")), snp.rs.id = read.gdsn(index.gdsn(gdsAN, "snp.rs.id"))))

########################## AA-BB Informative within the core ######
mapchrX <- subset(map, chr == 16) ## Select chromosome 16

mapchrX$position <- as.numeric(as.character(mapchrX$position))
mapchrX <- subset(mapchrX, position<= 60000000 & position>= 20000000)

#### Select only SNPs above certain Heterozygosity in inv-norm population
HetDf <- informativeSNPs

mapchrX <- mapchrX[mapchrX$snp.rs.id %in% HetDf$snp.rs.id,]

load(file=file.path(path.gen, "NormInveSamples.rda"))
allSamInv <- subset(allSam, stratum=="Inversion")
Samp <- as.character(allSamInv$sam)

registerDoMC(Ncor)

AAPerc  <- foreach(lo=1:length(Samp), .verbose=FALSE, .errorhandling=c('stop')) %dopar% {
  g <- snpgdsGetGeno(gdsAN, sample.id=Samp[[lo]], snp.id=mapchrX$snp.id)
  as.data.frame(cbind(Samp[[lo]], as.numeric(table(g)["0"])))
  
}

AAPerc <- rbindlist(AAPerc)
head(AAPerc)

pdf(file.path(path.gen,"AAInformativeWithinCore.pdf"))
AAWi <- hist(as.numeric(as.character(AAPerc$V2)), ylab ="Number of birds", xlab
             ="Number of AA genotypes")
plot(AAWi)

dev.off()


registerDoMC(Ncor)

BBPerc  <- foreach(lo=1:length(Samp), .verbose=FALSE, .errorhandling=c('stop')) %dopar% {
  g <- snpgdsGetGeno(gdsAN, sample.id=Samp[[lo]], snp.id=mapchrX$snp.id)
  as.data.frame(cbind(Samp[[lo]], as.numeric(table(g)["2"])))
  
}

BBPerc <- rbindlist(BBPerc)
head(AAPerc)

pdf(file.path(path.gen, "BBInformativeWithinCore.pdf"))
BBWi <- hist(as.numeric(as.character(AAPerc$V2)), ylab ="Number of birds", xlab
             ="Number of BB genotypes")
print(BBWi)
dev.off()

AAPercIn <- AAPerc
BBPercIn <- BBPerc

save(HetInv, HetNorm, informativeSNPs, 
     AAPercIn, BBPercIn, file=file.path(path.gen, "heterozygosity_inversion.rda"))
gdsfmt::closefn.gds(gdsAN)

#####  Select samples to generate the heatplot
# Open the GDS file
(gdsAN <- SNPRelate::snpgdsOpen(file.path(path.gen, "GT201420152016ChrNumeric.gds")))
load(file=file.path(path.gen, "heterozygosity_inversion.rda"))

#### Select 10 birds with distinct genotypes
AAPercIn$V2 <- as.numeric(as.character(AAPercIn$V2))
ten.haplo <- subset(AAPercIn, V2>50)
rest.inv <- subset(AAPercIn, V2<50)

g.inv.ten <- SNPRelate::snpgdsGetGeno(gdsAN, sample.id=ten.haplo$V1, 
                                       snp.id=informativeSNPs$snp.id)

set.seed(1000)
g.inv.rest <- SNPRelate::snpgdsGetGeno(gdsAN, sample.id=sample(rest.inv$V1, 10), 
                   snp.id=informativeSNPs$snp.id)

set.seed(1000)
g.norm <- SNPRelate::snpgdsGetGeno(gdsAN, sample.id=sample(subset(allSam, stratum=="Normal")$sam, 20), 
                                  snp.id=informativeSNPs$snp.id)

pdf(file.path(path.gen, "g.inv.ten.pdf"))
gplots::heatmap.2(g.inv.ten, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
                  density.info='none',  key=FALSE)
dev.off()

pdf(file.path(path.gen, "g.inv.rest.pdf"))
gplots::heatmap.2(g.inv.rest, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
                  density.info='none',  key=FALSE)
dev.off()

pdf(file.path(path.gen, "g.norm.pdf"))
gplots::heatmap.2(g.norm, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
                  density.info='none',  key=FALSE)
dev.off()

### Polishing 
g.inv.rest <- magick::image_read(file.path(path.gen, 'g.inv.rest.png'), density="300")
g.inv.ten <- magick::image_read(file.path(path.gen, 'g.inv.ten.png'), density="300")
g.norm <- magick::image_read(file.path(path.gen, 'g.norm.png'), density="300")

g.inv.rest <- magick::image_trim(g.inv.rest)
g.inv.rest <- magick::image_border(g.inv.rest,  "#000080", "10x10")
g.inv.ten <- magick::image_trim(g.inv.ten)
g.inv.ten <- magick::image_border(g.inv.ten,  "#000080", "10x10")
g.norm <- magick::image_trim(g.norm)
g.norm <- magick::image_border(g.norm,  "#000080", "10x10")

heat.all <- magick::image_append(c(g.inv.rest, g.inv.ten, g.norm))
heat.all  <- magick::image_border(magick::image_background(heat.all ,"hotpink"), "white", "0x50")
heat.all <- magick::image_annotate(heat.all, "A) inv-norm (+AB chr1A)", size = 25, location = "+80+0")
heat.all <- magick::image_annotate(heat.all, "B) inv-norm (+AA chr1A center only)", 
                                   size = 25, location = "+470+0")
heat.all <- magick::image_annotate(heat.all, "C) norm-norm (+AA chr1A)", size = 25, 
                                   location = "+960+0")

magick::image_write(heat.all, file.path(path.gen, "heatmap.png"), density=300)
