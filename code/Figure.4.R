#Methylation analysis -----
library(devtools)
library(ChAMP)
library(impute)
library("ChAMP")
setwd("D:/Urinary system tumors/work/2_Methylation_new")
library(data.table)
##ccRCC
a <- fread("D:/Urinary system tumors/work/2_Methylation_new/ccRCC/TCGA-KIRC.methylation450.tsv", data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta) 
betaData=beta$data
betaData=betaData+0.00001
sum(is.na(betaData))
pd <- fread("D:/Urinary system tumors/work/2_Methylation_new/ccRCC/TCGA-KIRC.clinical.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(betaData),]
pd$sample <- "Tumor"
pd$rowname <- row.names(pd)
pd$sample[grep("-11",pd$rowname)] <- "Normal"
pd <- pd[which(pd$gender.demographic =="male"),]
temp <- pd[,c("rowname","sample")]
betaData <- betaData[,row.names(pd)]
library(ChAMP)
myLoad=champ.filter(beta = betaData ,pd = pd)
dim(myLoad$beta)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K")
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na) 
table(num.na) 
myNorm <- myNorm[,which(num.na < 250000)] 
pdf("ccRCC_Methy_number_NA.pdf", width = 4, height = 4)
hist(num.na)
dev.off()

library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
pd <- myLoad$pd[colnames(myNorm),] 

group_list=pd$sample
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
pdf("ccRCC_Methy_PCA.pdf", width = 4, height = 4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

pheatmap(myNorm[cg,row.names(ac)[order(ac$group)]],show_colnames =F,show_rownames = F,cluster_cols = F,
         annotation_col=ac,filename = 'ccRCC_heatmap_top1000_sd.png')
dev.off()

pdf("ccRCC_Methy_Cor.pdf", width = 4, height = 4)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,filename = 'ccRCC_Methy_Cor.png')
dev.off()

group_list=pd$sample
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)

myDMR <- champ.DMR(beta=myNorm,pheno=group_list,method="Bumphunter")
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=group_list)
write.table(myDMR[["BumphunterDMR"]],"ccRCC_DMR.txt",sep = "\t",quote = F)

meth_gene <- myDMP$Tumor_to_Normal
meth_gene$Chromosome <- NA
meth_gene$ChromStart <- NA
meth_gene$ChromEnd <- NA

library(stringr)
for(i in 1:dim(meth_gene)[1]){
  meth_gene$Chromosome[i] <- strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][1]
  meth_gene$ChromStart[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][1]
  meth_gene$ChromEnd[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][2]
}


pheatmapdata <- meth_gene[,c(24:26,7:8,5,14,1)]
pheatmapdata <- na.omit(pheatmapdata)
pheatmapdata <- pheatmapdata[c(order(pheatmapdata$logFC,decreasing = T)[1:50],
                               order(pheatmapdata$logFC,decreasing = F)[1:50]),]
genedata <- meth_gene[,c(14,24:26)]
genedata <- pheatmapdata[,c(1:3)]
genedata$cg <- row.names(genedata)
genedata$ChromStart <- as.numeric(genedata$ChromStart)
genedata$ChromEnd <- as.numeric(genedata$ChromEnd)

pheatmapdata <- pheatmapdata[,c(1:5)]
colnames(pheatmapdata)[1:3] <- c('chr', 'start', 'end')
bed_long <- bedMatTolong(pheatmapdata)


library(ggcirclize)
library(ggnewscale)
col <- c("#B88586","#B7AD94","#BDCB95","#D4C091","#87A9C2","#68879D","#B6CDBB","#a3c6a3","#D0CFE2","#9392A5",
         "#d1d2d2","#fbd3b9","#7B907D","#F4B0B0","#B40CD9FF", "#EC9D5DFF", "#1b7837", "#0C4388FF", "#E2D60FFF")
#Circos Plot

pdf("ccRCC_Methylation.pdf", width = 8, height = 8)
ggcirclize(mapping = aes(start = 0, end = 360, genome = "hg38")) + 
  geom_trackgenomicrect(data = genedata,
                        aes(r0 = 0.9, r1 = 0.95, chr = Chromosome, gstart = ChromStart, gend = ChromEnd, fill = Chromosome),
                        color = NA) + 
  scale_fill_manual(values = c("Other" = "#f7f7f7", "WRKY" = "#000000")) + 
  geom_trackgenomiclabel(data = genedata,
                         aes(r0 = 0.8, r1 = 0.8, chr = Chromosome, gstart = ChromStart, gend = ChromEnd,label = cg),
                         strip.label = F,
                         label.size = 2.5,
                         link_col = col,
                         label.col = col) + 
  new_scale_color() + 
  new_scale_fill() + 
  geom_trackgenomictile(data = bed_long,
                        aes(r0 = 0.35, r1 = 0.6,
                            chr = chr,gstart = start,gend = end,
                            x = x, y = y, fill = value),
                        strip.label = F) + 
  scale_fill_gradient2(low = "#cbf3f0",mid = "#f7f7f7",high = "#ff9f1c",midpoint = 0.5) + 
  new_scale_color() + 
  new_scale_fill()
dev.off()
save.image("D:/Urinary system tumors/work/2_Methylation_new/ccRCC_Methylation.RData")




myGSEA <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = myDMR, arraytype = "450K",
                     adjPval = 0.05, method = "fisher")

x3 <- myGSEA[["DMR"]][1:20,]
x3$adjPval <- -log10(x3$adjPval)

x3$Gene_List <- factor(x3$Gene_List,levels = rev(x3$Gene_List) )
ggplot(x3,aes(Gene_List,adjPval))+
  geom_col(aes(fill=adjPval))+
  coord_flip()+
  coord_flip()+scale_fill_gradient( low = "#4F091C",
                                    high = "#CE5759")

ggsave("ccRCC_myGSEA.pdf",height = 5,width = 10)
write.table(myGSEA[["DMR"]],"ccRCC_myGSEA.txt",sep = "\t",quote = F,row.names = F)
save.image("ccRCC_Methylation.RData")
##BC
a <- fread("D:/Urinary system tumors/work/2_Methylation_new/BC/TCGA-BLCA.methylation450.tsv", data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta) 
betaData=beta$data
betaData=betaData+0.00001
sum(is.na(betaData))

pd <- fread("D:/Urinary system tumors/work/2_Methylation_new/BC/TCGA-BLCA.clinical.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(betaData),]
pd$sample <- "Tumor"
pd$rowname <- row.names(pd)
pd$sample[grep("-11",pd$rowname)] <- "Normal"
pd <- pd[which(pd$gender.demographic =="male"),]
temp <- pd[,c("rowname","sample")]
betaData <- betaData[,row.names(pd)]
library(ChAMP)
myLoad=champ.filter(beta = betaData ,pd = pd) 
dim(myLoad$beta)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K")
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na) 
table(num.na) 
myNorm <- myNorm[,which(num.na < 250000)]
pdf("BC_Methy_number_NA.pdf", width = 4, height = 4)
hist(num.na)
dev.off()

library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
pd <- myLoad$pd[colnames(myNorm),]

group_list=pd$sample
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
pdf("BC_Methy_PCA.pdf", width = 4, height = 4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

pheatmap(myNorm[cg,row.names(ac)[order(ac$group)]],show_colnames =F,show_rownames = F,cluster_cols = F,
         annotation_col=ac,filename = 'BC_heatmap_top1000_sd.png')
dev.off()

pdf("BC_Methy_Cor.pdf", width = 4, height = 4)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,filename = 'BC_Methy_Cor.png')
dev.off()
myDMR <- champ.DMR(beta=myNorm,pheno=group_list,method="Bumphunter")
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=group_list)
write.table(myDMR[["BumphunterDMR"]],"BC_DMR.txt",sep = "\t",quote = F)

myDMP <- champ.DMP(beta = myNorm,pheno=group_list)


meth_gene <- myDMP$Tumor_to_Normal
meth_gene$Chromosome <- NA
meth_gene$ChromStart <- NA
meth_gene$ChromEnd <- NA

library(stringr)
for(i in 1:dim(meth_gene)[1]){
  meth_gene$Chromosome[i] <- strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][1]
  meth_gene$ChromStart[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][1]
  meth_gene$ChromEnd[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][2]
}
pheatmapdata <- meth_gene[,c(24:26,7:8,5,14,1)]
pheatmapdata <- na.omit(pheatmapdata)
pheatmapdata <- pheatmapdata[c(order(pheatmapdata$logFC,decreasing = T)[1:50],
                               order(pheatmapdata$logFC,decreasing = F)[1:50]),]
genedata <- meth_gene[,c(14,24:26)]
genedata <- pheatmapdata[,c(1:3)]
genedata$cg <- row.names(genedata)
genedata$ChromStart <- as.numeric(genedata$ChromStart)
genedata$ChromEnd <- as.numeric(genedata$ChromEnd)

pheatmapdata <- pheatmapdata[,c(1:5)]
colnames(pheatmapdata)[1:3] <- c('chr', 'start', 'end')



library(ggcirclize)
library(ggnewscale)
bed_long <- bedMatTolong(pheatmapdata)

col <- c("#B88586","#B7AD94","#BDCB95","#D4C091","#87A9C2","#68879D","#B6CDBB","#a3c6a3","#D0CFE2","#9392A5",
         "#d1d2d2","#fbd3b9","#7B907D","#F4B0B0","#B40CD9FF", "#EC9D5DFF", "#1b7837", "#0C4388FF", "#E2D60FFF")
###Circos Plot

pdf("BC_Methylation.pdf", width = 8, height = 8)
ggcirclize(mapping = aes(start = 0, end = 360, genome = "hg38")) + 
  geom_trackgenomicrect(data = genedata,
                        aes(r0 = 0.9, r1 = 0.95, chr = Chromosome, gstart = ChromStart, gend = ChromEnd, fill = Chromosome),
                        color = NA) + 
  scale_fill_manual(values = c("Other" = "#f7f7f7", "WRKY" = "#000000")) + 
  geom_trackgenomiclabel(data = genedata,
                         aes(r0 = 0.8, r1 = 0.8, chr = Chromosome, gstart = ChromStart, gend = ChromEnd,label = cg),
                         strip.label = F,
                         label.size = 2.5,
                         link_col = col,
                         label.col = col) + 
  new_scale_color() + 
  new_scale_fill() + 
  geom_trackgenomictile(data = bed_long,
                        aes(r0 = 0.35, r1 = 0.6,
                            chr = chr,gstart = start,gend = end,
                            x = x, y = y, fill = value),
                        strip.label = F) + 
  scale_fill_gradient2(low = "#cbf3f0",mid = "#f7f7f7",high = "#ff9f1c",midpoint = 0.5) + 
  new_scale_color() + 
  new_scale_fill()
dev.off()


myGSEA <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = myDMR, arraytype = "450K",
                     adjPval = 0.05, method = "fisher")

x3 <- myGSEA[["DMR"]][1:20,]
x3$adjPval <- -log10(x3$adjPval)

x3$Gene_List <- factor(x3$Gene_List,levels = rev(x3$Gene_List) )
ggplot(x3,aes(Gene_List,adjPval))+
  geom_col(aes(fill=adjPval))

ggsave("BC_myGSEA.pdf",height = 5,width = 10)
write.table(myGSEA[["DMR"]],"BC_myGSEA.txt",sep = "\t",quote = F,row.names = F)
##PCa
library(devtools)
library(ChAMP)
library(impute)
library("ChAMP")
##
setwd("D:/Urinary system tumors/work/2_Methylation")
library(data.table)
a <- fread("D:/Urinary system tumors/work/2_Methylation/PCa/TCGA-PRAD.methylation450.tsv", data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
#
beta=impute.knn(beta) 
betaData=beta$data
betaData=betaData+0.00001
sum(is.na(betaData))

#
pd <- fread("D:/Urinary system tumors/work/2_Methylation/PCa/TCGA-PRAD.GDC_phenotype.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(betaData),]
table(pd$gleason_score) 
pd$sample <- "Tumor"
pd$sample[grep("-11",pd$submitter_id.samples)] <- "Normal"
library(ChAMP)
myLoad=champ.filter(beta = betaData ,pd = pd) 
dim(myLoad$beta)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K")
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na) 
table(num.na) 
myNorm <- myNorm[,which(num.na < 250000)] 
pdf("PCa_Methy_number_NA.pdf", width = 4, height = 4)
hist(num.na)
dev.off()


library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
pd <- myLoad$pd[colnames(myNorm),] 

group_list=pd$sample
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
pdf("PCa_Methy_PCA.pdf", width = 4, height = 4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

pheatmap(myNorm[cg,row.names(ac)[order(ac$group)]],show_colnames =F,show_rownames = F,cluster_cols = F,
         annotation_col=ac,filename = 'PCa_heatmap_top1000_sd.png')
dev.off()

pdf("PCa_Methy_Cor.pdf", width = 4, height = 4)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,filename = 'PCa_Methy_Cor.png')
dev.off()

group_list=pd$group
group_list=pd$sample
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
write.table(myDMP$Normal_to_Tumor,"PCa_NvsT.txt",quote = F,sep = "\t")

PCa_myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
DMR.GUI(DMR=PCa_myDMR,beta=myNorm,pheno=group_list)

meth_gene <- myDMP$Normal_to_Tumor
meth_gene$Chromosome <- NA
meth_gene$ChromStart <- NA
meth_gene$ChromEnd <- NA

library(stringr)
for(i in 1:dim(meth_gene)[1]){
  meth_gene$Chromosome[i] <- strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][1]
  meth_gene$ChromStart[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][1]
  meth_gene$ChromEnd[i] <- strsplit(strsplit(as.character(meth_gene$UCSC_CpG_Islands_Name[i]),":")[[1]][2],"-")[[1]][2]
}
#
save.image("D:/Urinary system tumors/work/2_Methylation/PCa_Methylation.RData")
write.csv(meth_gene,"PCa_Methylation_Normal_to_Tumor.csv",quote = F)

pheatmapdata <- meth_gene[,c(24:26,7:8,5,14,1)]
pheatmapdata <- na.omit(pheatmapdata)
pheatmapdata <- pheatmapdata[c(order(pheatmapdata$logFC,decreasing = T)[1:50],
                               order(pheatmapdata$logFC,decreasing = F)[1:50]),]
genedata <- meth_gene[,c(14,24:26)]
genedata <- pheatmapdata[,c(1:3)]
genedata$cg <- row.names(genedata)
genedata$ChromStart <- as.numeric(genedata$ChromStart)
genedata$ChromEnd <- as.numeric(genedata$ChromEnd)

pheatmapdata <- pheatmapdata[,c(1:5)]
colnames(pheatmapdata)[1:3] <- c('chr', 'start', 'end')
bed_long <- bedMatTolong(pheatmapdata)


write.table(pheatmapdata$gene,"PCa_plot_gene.txt",quote = F,row.names = F,sep = "\t")

group_gene_PPI_data <- read.table("PCa_plot_genestring_interactions_short.tsv",header = F)
group_gene_PPI_data <- group_gene_PPI_data[,c(1,2)]
colnames(group_gene_PPI_data) <- c("node1","node2")
gene_pos <- pheatmapdata[,c(7,1:3)]
library(dplyr)
gene_pos <- distinct(gene_pos)

group_gene_PPI_data <- group_gene_PPI_data[group_gene_PPI_data$node1 %in% gene_pos$gene,]
group_gene_PPI_data <- group_gene_PPI_data[group_gene_PPI_data$node2 %in% gene_pos$gene,]

group_gene_PPI_data$chr1 <- NA
group_gene_PPI_data$start1 <- NA
group_gene_PPI_data$end1 <- NA

group_gene_PPI_data$chr2 <- NA
group_gene_PPI_data$start2 <- NA
group_gene_PPI_data$end2 <- NA

for(i in 1:length(group_gene_PPI_data$node1)){
  group_gene_PPI_data[i,c(3:5)] <- gene_pos[which(gene_pos$gene == group_gene_PPI_data$node1[i]),c(2:4)]
  group_gene_PPI_data[i,c(6:8)] <- gene_pos[which(gene_pos$gene == group_gene_PPI_data$node2[i]),c(2:4)]
}
library(ggcirclize)
library(ggnewscale)

col <- c("#B88586","#B7AD94","#BDCB95","#D4C091","#87A9C2","#68879D","#B6CDBB","#a3c6a3","#D0CFE2","#9392A5",
         "#d1d2d2","#fbd3b9","#7B907D","#F4B0B0","#B40CD9FF", "#EC9D5DFF", "#1b7837", "#0C4388FF", "#E2D60FFF")
####Circos Plot

pdf("PCa_Methylation.pdf", width = 8, height = 8)
ggcirclize(mapping = aes(start = 0, end = 360, genome = "hg38")) + 
  geom_trackgenomicrect(data = genedata,
                        aes(r0 = 0.9, r1 = 0.95, chr = Chromosome, gstart = ChromStart, gend = ChromEnd, fill = Chromosome),
                        color = NA) + 
  scale_fill_manual(values = c("Other" = "#f7f7f7", "WRKY" = "#000000")) + 
  geom_trackgenomiclabel(data = genedata,
                         aes(r0 = 0.8, r1 = 0.8, chr = Chromosome, gstart = ChromStart, gend = ChromEnd,label = cg),
                         strip.label = F,
                         label.size = 2.5,
                         link_col = col,
                         label.col = col) + 
  new_scale_color() + 
  new_scale_fill() + 
  geom_trackgenomictile(data = bed_long,
                        aes(r0 = 0.35, r1 = 0.6,
                            chr = chr,gstart = start,gend = end,
                            x = x, y = y, fill = value),
                        strip.label = F) + 
  scale_fill_gradient2(low = "#cbf3f0",mid = "#f7f7f7",high = "#ff9f1c",midpoint = 0.5) + 
  new_scale_color() + 
  new_scale_fill()
dev.off()


myGSEA <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = PCa_myDMR, arraytype = "450K",
                     adjPval = 0.05, method = "fisher")

x3 <- myGSEA[["DMR"]][1:20,]
x3$adjPval <- -log10(x3$adjPval)

x3$Gene_List <- factor(x3$Gene_List,levels = rev(x3$Gene_List) )
ggplot(x3,aes(Gene_List,adjPval))+
  geom_col(aes(fill=adjPval))+
  coord_flip()+
  coord_flip()+scale_fill_gradient( low = "#08605B",
                                    high = "#8DE2DE")

ggsave("PCa_myGSEA.pdf",height = 5,width = 10)
write.table(myGSEA[["DMR"]],"PCa_myGSEA.txt",sep = "\t",quote = F,row.names = F)

#SNP analysis -----
setwd("D:/Urinary system tumors/work/2_SNP")
library(maftools)
load(file = "./GDCdata/TCGA-PRAD_SNP.Rdata")
maf.PRAD <- data
maf.PCa <- read.maf(maf.PRAD)
rm("data")
load(file = "./GDCdata/TCGA-KIRC_SNP.Rdata")
load("D:/Urinary system tumors/work/2_SNP/ccRCC_TCGA_clinical.RData")
maf.KIRC <- data
maf.KIRC <- data[which(substr(maf.KIRC$Tumor_Sample_Barcode, 1, 12) %in% clinical[which(clinical$gender == "MALE"),]$bcr_patient_barcode),]
maf.ccRCC <- read.maf(maf.KIRC)
rm("data")
rm("clinical")
load(file = "./GDCdata/TCGA-BLCA_SNP.Rdata")
load("D:/Urinary system tumors/work/2_SNP/BC_TCGA_clinical.RData")
maf.BLCA <- data
maf.BLCA <- data[which(substr(maf.BLCA$Tumor_Sample_Barcode, 1, 12) %in% clinical[which(clinical$gender == "MALE"),]$bcr_patient_barcode),]
maf.BC <- read.maf(maf.BLCA)
rm("data")
rm("clinical")

plotmafSummary(maf = maf.BC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
plotmafSummary(maf = maf.ccRCC, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
plotmafSummary(maf = maf.PCa, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

BC_SNP_data <- maf.BC@data
BC_SNP_data <- BC_SNP_data[which(BC_SNP_data$Variant_Classification != "Nonsense_Mutation"),]

ccRCC_SNP_data <- maf.ccRCC@data
ccRCC_SNP_data <- ccRCC_SNP_data[which(ccRCC_SNP_data$Variant_Classification != "Nonsense_Mutation"),]

PCa_SNP_data <- maf.PCa@data
PCa_SNP_data <- PCa_SNP_data[which(PCa_SNP_data$Variant_Classification != "Nonsense_Mutation"),]

save(BC_SNP_data, ccRCC_SNP_data, PCa_SNP_data,
     file = "SNP_data.Rdata")

#CNV-peak-DMR -----
library(ArchR)
library(dplyr)
DMR <- read.table("D:/Urinary system tumors/work/2_Methylation_new/ccRCC_DMR.txt",header = T,sep = "\t")
DMR$methy_type <- "Down"
DMR[which(DMR$value > 0),]$methy_type <- "Up"
library(GenomicRanges)
DMR_GRanges <- GRanges(seqnames = DMR$seqnames,
                       ranges = IRanges(start = DMR$start, end = DMR$end))

peak <- read.csv("D:/ccRCC_scATAC/ccRCC_motif_peak_gene.csv",
                 header=T, sep = ",")
peak <- peak[which(peak$motif %in% c("FOSL1_142",
                                     "RUNX1_733",
                                     "NFIB_741",
                                     "NFIC_740",
                                     "CEBPB_140",
                                     "HNF4A_662",
                                     "ID4_75")),]
x <- as.data.frame(strsplit(peak$peakName,":"))
y <- as.data.frame(strsplit(as.character(x[2,]),"-"))
colnames(y) <- colnames(x)
z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
peak_GRanges <- GRanges(seqnames = z$seqnames,
                        ranges = IRanges(start = z$start, end = z$end))
peak_methy <- findOverlaps(peak_GRanges, DMR_GRanges)
peak_methy <- as.data.frame(peak_methy)

z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "_")
peak_methy$peak <- z[peak_methy$queryHits]

z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
zz <- z
z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "-")
peak_methy$peak2 <- z[peak_methy$queryHits]



z <- paste(paste(DMR$seqnames,DMR$start,sep = ":"),DMR$end,sep = "-")
peak_methy$methy <- z[peak_methy$subjectHits]

DMR2 <- DMR[peak_methy$subjectHits,]
peak2 <- peak[which(peak$peakName %in% peak_methy$peak2),]
CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
                  header = T, sep = "\t")
CNV<- CNV %>% dplyr::mutate(CNV=recode(state,
                                       "1"="Complete loss",
                                       "2"="One copy_loss",
                                       "3"="Neutral",
                                       "4"="One copy_gain",
                                       "5"="Two copies_gain",
                                       "6"=">Two copies_gain"))
CNV <- CNV[which(CNV$CNV != "Neutral"),]
CNV <- CNV[grepl("Tumor", CNV$cell_group_name),4:7]
CNV <- distinct(CNV)
CNV$name <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV_GRanges <- GRanges(seqnames = CNV$chr,
                       ranges = IRanges(start = CNV$start, end = CNV$end))
peak_CNV <- findOverlaps(peak_GRanges, CNV_GRanges)
peak_CNV <- as.data.frame(peak_CNV)

peak_CNV$peak <- paste(paste(zz$seqnames,zz$start,sep = ":"),zz$end,sep = "_")[peak_CNV$queryHits]
peak_CNV$CNV <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")[peak_CNV$subjectHits]
peak_CNV2 <- peak_CNV[,c(3,4)]
peak_CNV2 <- distinct(peak_CNV2)
peak_CNV2 <- peak_CNV2[which(peak_CNV2$peak %in% peak_methy$peak),]
CNV$CNV2 <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV2 <- CNV[CNV$CNV2 %in% peak_CNV2$CNV,]
write.table(CNV2,"Pca_peak_cnv.txt",quote = F,row.names = F,sep = "\t")

library(ArchR)
library(dplyr)
DMR <- read.table("D:/Urinary system tumors/work/2_Methylation_new/PCa_DMR.txt",header = T,sep = "\t")
DMR$methy_type <- "Down"
DMR[which(DMR$value > 0),]$methy_type <- "Up"
library(GenomicRanges)
DMR_GRanges <- GRanges(seqnames = DMR$seqnames,
                       ranges = IRanges(start = DMR$start, end = DMR$end))

peak <- read.csv("D:/PCa_ATAC/PCa_motif_peak_gene.csv",
                 header=T, sep = ",")
x <- as.data.frame(strsplit(peak$peakName,":"))
y <- as.data.frame(strsplit(as.character(x[2,]),"-"))
colnames(y) <- colnames(x)
z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
zz <- z
peak_GRanges <- GRanges(seqnames = z$seqnames,
                        ranges = IRanges(start = z$start, end = z$end))
peak_methy <- findOverlaps(peak_GRanges, DMR_GRanges)
peak_methy <- as.data.frame(peak_methy)

z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "_")
peak_methy$peak <- z[peak_methy$queryHits]
z <- paste(paste(DMR$seqnames,DMR$start,sep = ":"),DMR$end,sep = "-")
peak_methy$methy <- z[peak_methy$subjectHits]

DMR2 <- DMR[peak_methy$subjectHits,]


z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "-")
peak_methy$peak2 <- z[peak_methy$queryHits]

peak2 <- peak[which(peak$peakName %in% peak_methy$peak2),]
CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
                  header = T, sep = "\t")
CNV<- CNV %>% dplyr::mutate(CNV=recode(state,
                                       "1"="Complete loss",
                                       "2"="One copy_loss",
                                       "3"="Neutral",
                                       "4"="One copy_gain",
                                       "5"="Two copies_gain",
                                       "6"=">Two copies_gain"))
CNV <- CNV[which(CNV$CNV != "Neutral"),]
CNV <- CNV[grepl("Tumor", CNV$cell_group_name),4:7]
CNV <- distinct(CNV)
CNV$name <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV_GRanges <- GRanges(seqnames = CNV$chr,
                       ranges = IRanges(start = CNV$start, end = CNV$end))
peak_CNV <- findOverlaps(peak_GRanges, CNV_GRanges)
peak_CNV <- as.data.frame(peak_CNV)

peak_CNV$peak <- paste(paste(zz$seqnames,zz$start,sep = ":"),zz$end,sep = "_")[peak_CNV$queryHits]
peak_CNV$CNV <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")[peak_CNV$subjectHits]
peak_CNV2 <- peak_CNV[,c(3,4)]
peak_CNV2 <- distinct(peak_CNV2)
peak_CNV2 <- peak_CNV2[which(peak_CNV2$peak %in% peak_methy$peak),]
CNV$CNV2 <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV2 <- CNV[CNV$CNV2 %in% peak_CNV2$CNV,]
write.table(CNV2,"Pca_peak_cnv.txt",quote = F,row.names = F,sep = "\t")

met_cnv<- findOverlaps(DMR_GRanges, CNV_GRanges)
met_cnv <- as.data.frame(met_cnv)
met_cnv$met <- paste(paste(DMR$seqnames,DMR$start,sep = ":"),DMR$end,sep = "-")[met_cnv$queryHits]
met_cnv$cnv <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")[met_cnv$subjectHits]
length(unique(met_cnv$met))
length(unique(met_cnv$cnv))

write.table(peak_methy,"peak_methy.txt",quote = F,row.names = F,sep = "\t")
write.table(peak2,"peak_methy_peak.txt",quote = F,row.names = F,sep = "\t")
write.table(DMR2,"peak_methy_DMR.txt",quote = F,row.names = F,sep = "\t")
#chrX and chrY Methylation------
load("D:/Urinary system tumors/work/2_Methylation_new/ccRCC_Methylation.RData")
library(ArchR)
library(dplyr)
library(ChAMP)
myLoad=champ.filter(beta = betaData ,pd = pd,filterXY = F)
dim(myLoad$beta)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K") 
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na) 
table(num.na) 
myNorm <- myNorm[,which(num.na < 250000)] 
pdf("ccRCC_Methy_number_NA.pdf", width = 4, height = 4)
hist(num.na)
dev.off()

library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
pd <- myLoad$pd[colnames(myNorm),] 

group_list=pd$sample
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
pdf("ccRCC_Methy_PCA.pdf", width = 4, height = 4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

pheatmap(myNorm[cg,row.names(ac)[order(ac$group)]],show_colnames =F,show_rownames = F,cluster_cols = F,
         annotation_col=ac,filename = 'ccRCC_heatmap_top1000_sd.png')
dev.off()

pdf("ccRCC_Methy_Cor.pdf", width = 4, height = 4)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,filename = 'ccRCC_Methy_Cor.png')
dev.off()

group_list=pd$sample
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)

myDMR <- champ.DMR(beta=myNorm,pheno=group_list,method="Bumphunter")
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=group_list)
write.table(myDMR[["BumphunterDMR"]],"ccRCC_DMR.txt",sep = "\t",quote = F)
DMR <- read.table("ccRCC_DMR.txt",header = T,sep = "\t")
DMR$methy_type <- "Down"
DMR[which(DMR$value > 0),]$methy_type <- "Up"
library(GenomicRanges)
DMR_GRanges <- GRanges(seqnames = DMR$seqnames,
                       ranges = IRanges(start = DMR$start, end = DMR$end))
peak <- read.csv("D:/ccRCC_scATAC/ccRCC_motif_peak_gene.csv",
                 header=T, sep = ",")
peak <- peak[which(peak$motif %in% c("FOSL1_142",
                                     "RUNX1_733",
                                     "NFIB_741",
                                     "NFIC_740",
                                     "CEBPB_140",
                                     "HNF4A_662",
                                     "ID4_75")),]
x <- as.data.frame(strsplit(peak$peakName,":"))
y <- as.data.frame(strsplit(as.character(x[2,]),"-"))
colnames(y) <- colnames(x)
z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
peak_GRanges <- GRanges(seqnames = z$seqnames,
                        ranges = IRanges(start = z$start, end = z$end))
peak_methy <- findOverlaps(peak_GRanges, DMR_GRanges)
peak_methy <- as.data.frame(peak_methy)

z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "_")
peak_methy$peak <- z[peak_methy$queryHits]

z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
zz <- z
z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "-")
peak_methy$peak2 <- z[peak_methy$queryHits]



z <- paste(paste(DMR$seqnames,DMR$start,sep = ":"),DMR$end,sep = "-")
peak_methy$methy <- z[peak_methy$subjectHits]

DMR2 <- DMR[peak_methy$subjectHits,]
peak2 <- peak[which(peak$peakName %in% peak_methy$peak2),]
CNV <- read.table("./inferCNV_XY/out_result_cnv_ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
                  header = T, sep = "\t")
CNV<- CNV %>% dplyr::mutate(CNV=recode(state,
                                       "1"="Complete loss",
                                       "2"="One copy_loss",
                                       "3"="Neutral",
                                       "4"="One copy_gain",
                                       "5"="Two copies_gain",
                                       "6"=">Two copies_gain"))
CNV <- CNV[which(CNV$CNV != "Neutral"),]
CNV <- CNV[grepl("Tumor", CNV$cell_group_name),4:7]
CNV <- distinct(CNV)
CNV$name <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV_GRanges <- GRanges(seqnames = CNV$chr,
                       ranges = IRanges(start = CNV$start, end = CNV$end))
peak_CNV <- findOverlaps(peak_GRanges, CNV_GRanges)
peak_CNV <- as.data.frame(peak_CNV)

peak_CNV$peak <- paste(paste(zz$seqnames,zz$start,sep = ":"),zz$end,sep = "_")[peak_CNV$queryHits]
peak_CNV$CNV <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")[peak_CNV$subjectHits]
peak_CNV2 <- peak_CNV[,c(3,4)]
peak_CNV2 <- distinct(peak_CNV2)
peak_CNV2 <- peak_CNV2[which(peak_CNV2$peak %in% peak_methy$peak),]
CNV$CNV2 <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV2 <- CNV[CNV$CNV2 %in% peak_CNV2$CNV,]

write.table(peak_CNV2,"ccRCC_peak_cnv.txt",quote = F,row.names = F,sep = "\t")
write.table(peak2[which(peak2$peakName == "chrX:49043000-49043499"),],"ccRCC_peak_network.txt",quote = F,row.names = F,sep = "\t")
write.table(DMR2[which(DMR2$seqnames == "chrX"),],"ccRCC_DMR_cnv.txt",quote = F,row.names = F,sep = "\t")

load("D:/Urinary system tumors/work/2_Methylation_new/PCa_Methylation.RData")

pd <- fread("D:/Urinary system tumors/work/2_Methylation/PCa/TCGA-PRAD.GDC_phenotype.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(betaData),]
table(pd$gleason_score) 
pd$sample <- "Tumor"
pd$sample[grep("-11",pd$submitter_id.samples)] <- "Normal"
library(ChAMP)
myLoad=champ.filter(beta = betaData ,pd = pd, filterXY=F) 
dim(myLoad$beta)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm)
QC.GUI(beta=myNorm,arraytype="450K") 
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
hist(num.na) 
table(num.na) 
myNorm <- myNorm[,which(num.na < 250000)] 
pdf("PCa_Methy_number_NA.pdf", width = 4, height = 4)
hist(num.na)
dev.off()


library("FactoMineR")
library("factoextra") 
dat <- t(myNorm)
pd <- myLoad$pd[colnames(myNorm),] 

group_list=pd$sample
table(group_list)
dat.pca <- PCA(dat, graph = FALSE) 
pdf("PCa_Methy_PCA.pdf", width = 4, height = 4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE,
             legend.title = "Groups")
dev.off()

cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)

pheatmap(myNorm[cg,row.names(ac)[order(ac$group)]],show_colnames =F,show_rownames = F,cluster_cols = F,
         annotation_col=ac,filename = 'PCa_heatmap_top1000_sd.png')
dev.off()

pdf("PCa_Methy_Cor.pdf", width = 4, height = 4)
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F,filename = 'PCa_Methy_Cor.png')
dev.off()


group_list=pd$group
group_list=pd$sample
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
PCa_myDMR <- champ.DMR(beta = myNorm,pheno=group_list,method="Bumphunter")
DMR.GUI(DMR=PCa_myDMR,beta=myNorm,pheno=group_list)
write.table(PCa_myDMR[["BumphunterDMR"]],"PCa_DMR.txt",sep = "\t",quote = F)

DMR <- read.table("PCa_DMR.txt",header = T,sep = "\t")
DMR$methy_type <- "Down"
DMR[which(DMR$value > 0),]$methy_type <- "Up"
library(GenomicRanges)
DMR_GRanges <- GRanges(seqnames = DMR$seqnames,
                       ranges = IRanges(start = DMR$start, end = DMR$end))

peak <- read.csv("D:/PCa_ATAC/PCa_motif_peak_gene.csv",
                 header=T, sep = ",")
x <- as.data.frame(strsplit(peak$peakName,":"))
y <- as.data.frame(strsplit(as.character(x[2,]),"-"))
colnames(y) <- colnames(x)
z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
zz <- z
peak_GRanges <- GRanges(seqnames = z$seqnames,
                        ranges = IRanges(start = z$start, end = z$end))
peak_methy <- findOverlaps(peak_GRanges, DMR_GRanges)
peak_methy <- as.data.frame(peak_methy)

z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "_")
peak_methy$peak <- z[peak_methy$queryHits]
z <- paste(paste(DMR$seqnames,DMR$start,sep = ":"),DMR$end,sep = "-")
peak_methy$methy <- z[peak_methy$subjectHits]

DMR2 <- DMR[peak_methy$subjectHits,]


z <- as.data.frame(t(rbind(x[1,],y[c(1,2),])))
z$V2 <- as.numeric(z$V2)
z$V3 <- as.numeric(z$V3)
View(z)
z <- distinct(z)
colnames(z) <- c("seqnames", "start", "end")
z <- paste(paste(z$seqnames,z$start,sep = ":"),z$end,sep = "-")
peak_methy$peak2 <- z[peak_methy$queryHits]

peak2 <- peak[which(peak$peakName %in% peak_methy$peak2),]
CNV <- read.table("./inferCNV_XY/out_result_cnv_PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
                  header = T, sep = "\t")
CNV<- CNV %>% dplyr::mutate(CNV=recode(state,
                                       "1"="Complete loss",
                                       "2"="One copy_loss",
                                       "3"="Neutral",
                                       "4"="One copy_gain",
                                       "5"="Two copies_gain",
                                       "6"=">Two copies_gain"))
CNV <- CNV[which(CNV$CNV != "Neutral"),]
CNV <- CNV[grepl("Tumor", CNV$cell_group_name),4:7]
CNV <- distinct(CNV)
CNV$name <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV_GRanges <- GRanges(seqnames = CNV$chr,
                       ranges = IRanges(start = CNV$start, end = CNV$end))
peak_CNV <- findOverlaps(peak_GRanges, CNV_GRanges)
peak_CNV <- as.data.frame(peak_CNV)

peak_CNV$peak <- paste(paste(zz$seqnames,zz$start,sep = ":"),zz$end,sep = "_")[peak_CNV$queryHits]
peak_CNV$CNV <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")[peak_CNV$subjectHits]
peak_CNV2 <- peak_CNV[,c(3,4)]
peak_CNV2 <- distinct(peak_CNV2)
peak_CNV2 <- peak_CNV2[which(peak_CNV2$peak %in% peak_methy$peak),]
CNV$CNV2 <- paste(paste(CNV$chr,CNV$start,sep = ":"),CNV$end,sep = "-")
CNV2 <- CNV[CNV$CNV2 %in% peak_CNV2$CNV,]

#show of scATAC ----
##PCa
load("D:\\PCa_ATAC\\CancerTrajector.RData")
p4 <- plotTrajectoryHeatmap(trajPM, 
                            labelMarkers = peak_methy$peak,#read peak_methy.txt of PCa
                            labelTop = 0,
                            varCutOff = 0,
                            maxFeatures = 170000,
                            pal = paletteContinuous(set = "solarExtra"), 
                            colorColumns = TRUE)
mat <- trajPM@assays@data@listData[["mat"]]
plotPDF(p4, name = "peak_DMR.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)

#Select genes and regions based on Supplementary Table S4
#peak:chr21:44353500-44353999
markerGenes  <- c("KCNK15","TFAP2A")
#SNP chr21:44353851-44353851
#CNV chr21:17593653-46665124
highlight_GRanges <- GRanges(seqnames = c("chr20","chr6"),
                             ranges = IRanges(start = c(44747000,10419500), end = c(44747499,10419999) ))

color <- c("G3" = "#D8EBD3",
           "G4" = "#228A87")
p <- plotBrowserTrack(
  ArchRProj = proj_CancerTrajector_subset, 
  groupBy = "group", 
  geneSymbol = markerGenes, 
  upstream = 20000,
  downstream = 20000,
  loops = getPeak2GeneLinks(proj_CancerTrajector_subset),
  highlight = highlight_GRanges,
  pal = color
)
#grid::grid.newpage()
#grid::grid.draw(p$TRPM2)
plotPDF(plotList = p, 
        name = "gene.pdf", 
        ArchRProj = proj_CancerTrajector_subset, 
        addDOC = FALSE, width = 5, height = 5)
##ccRCC
load("./ccRCC_CancerTrajector.RData")
#Select genes and regions based on Supplementary Table S4
#peak:chr21:44353500-44353999
markerGenes  <- c("DSP")
#SNP chr21:44353851-44353851
#CNV chr21:17593653-46665124
highlight_GRanges <- GRanges(seqnames = c("chr6"),
                             ranges = IRanges(start = c(7541442,7580337,7583785,7570555), end = c(7542126,7580337,7583785,7570555) ))

color <- c("Normal" = "#1CB6C0",
           'I'      = "#F8CED7",
           "II"     = "#D66E71",
           "III"    = "#AA1F24")
p <- plotBrowserTrack(
  ArchRProj = proj_epi, 
  groupBy = "WHO", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj_epi),
  highlight = highlight_GRanges,
  pal = color
)
#grid::grid.newpage()
#grid::grid.draw(p$TRPM2)
plotPDF(plotList = p, 
        name = "gene_DSP.pdf", 
        ArchRProj = proj_epi, 
        addDOC = FALSE, width = 5, height = 5)