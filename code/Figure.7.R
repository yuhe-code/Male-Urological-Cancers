# Drug filter ----
#drug2cell
library(ggpubr)
#library(car)
library(outliers)
library(ggbreak)


# 4.2 Seurat to h5ad
#BC-----

sce <- readRDS("……/BC_tumor.rds")

DimPlot(sce,reduction = "umap",group.by = "cell_type_second")


reticulate::py_discover_config() 
reticulate::use_condaenv("……/miniconda3/envs/Drug2cell")
time.R2py = system.time({
  sceasy::convertFormat(sce, from="seurat", to="anndata",
                        outFile='……/sce_BC.h5ad')
})
print(time.R2py)


#ccRCC-----

sce <- readRDS("……/ccRCC_tumor.rds")

DimPlot(sce,reduction = "umap",group.by = "cell_type_second")


reticulate::py_discover_config() 
reticulate::use_condaenv("……/miniconda3/envs/Drug2cell")
time.R2py = system.time({
  sceasy::convertFormat(sce, from="seurat", to="anndata",
                        outFile='……/sce_ccRCC.h5ad')
})
print(time.R2py)


#PCa-----
sce <- readRDS("……/PCa_tumor.rds")

DimPlot(sce,reduction = "umap",group.by = "cell_type_second")


reticulate::py_discover_config() 
reticulate::use_condaenv("……/miniconda3/envs/Drug2cell")
time.R2py = system.time({
  sceasy::convertFormat(sce, from="seurat", to="anndata",
                        outFile='……/sce_PCa.h5ad')
})
print(time.R2py)


#python
import scanpy as sc
import pandas as pd
import drug2cell as d2c
#BC----
adata_BC = sc.read_h5ad("……/sce_BC.h5ad")
sc.settings.set_figure_params(dpi=150)
sc.pl.umap(adata_BC, color="cell_type_second")
adata_BC

d2c.score(adata_BC, use_raw=True)


sc.tl.rank_genes_groups(adata_BC.uns['drug2cell'], method="wilcoxon", groupby="cell_type_second")

sc.pl.rank_genes_groups_dotplot(adata_BC.uns['drug2cell'], swap_axes=True, dendrogram=False, n_genes=20)


adata_BC.uns['drug2cell'].var['genes'].to_csv('(……/BC_drug_targets.csv')


celltype_BC=adata_BC.uns['drug2cell'].obs['cell_type_second'].unique().tolist()
degs_BC = sc.get.rank_genes_groups_df(adata_BC.uns['drug2cell'],group= celltype_BC)
degs_BC.to_csv('(……/BC_drug_score.csv')

sc.pl.rank_genes_groups_dotplot(adata_BC.uns['drug2cell'],
                                swap_axes=True,
                                dendrogram=False, 
                                n_genes=5,
                                values_to_plot="scores")

sc.write("……/sce_BC.h5ad", adata_BC)


#ccRCC----
adata_ccRCC = sc.read_h5ad("……/sce_ccRCC.h5ad")
sc.pl.umap(adata_ccRCC, color="cell_type_second")
adata_ccRCC


d2c.score(adata_ccRCC, use_raw=True)


sc.tl.rank_genes_groups(adata_ccRCC.uns['drug2cell'], method="wilcoxon", groupby="cell_type_second")


sc.pl.rank_genes_groups_dotplot(adata_ccRCC.uns['drug2cell'], swap_axes=True, dendrogram=False, n_genes=20)


adata_ccRCC.uns['drug2cell'].var['genes'].to_csv("……/ccRCC_drug_targets.csv")


celltype_ccRCC=adata_ccRCC.uns['drug2cell'].obs['cell_type_second'].unique().tolist()
degs_ccRCC = sc.get.rank_genes_groups_df(adata_ccRCC.uns['drug2cell'],group= celltype_ccRCC)
degs_ccRCC.to_csv('……/ccRCC_drug_score.csv')


sc.pl.rank_genes_groups_dotplot(adata_ccRCC.uns['drug2cell'],
                                swap_axes=True,
                                dendrogram=False, 
                                n_genes=5,
                                values_to_plot="scores")

sc.write("……/sce_ccRCC.h5ad", adata_ccRCC)

#PCa----

adata_PCa = sc.read_h5ad("……/sce_PCa.h5ad")

sc.pl.umap(adata_PCa, color="cell_type_second")

adata_PCa

d2c.score(adata_PCa, use_raw=True)


sc.tl.rank_genes_groups(adata_PCa.uns['drug2cell'], method="wilcoxon", groupby="cell_type_second")


sc.pl.rank_genes_groups_dotplot(adata_PCa.uns['drug2cell'], swap_axes=True, dendrogram=False, n_genes=20)



adata_PCa.uns['drug2cell'].var['genes'].to_csv('……/PCa_drug_targets.csv')

celltype_PCa=adata_PCa.uns['drug2cell'].obs['cell_type_second'].unique().tolist()
degs_PCa = sc.get.rank_genes_groups_df(adata_PCa.uns['drug2cell'],group= celltype_PCa)
degs_PCa.to_csv('……/PCa_drug_score.csv')

sc.pl.rank_genes_groups_dotplot(adata_PCa.uns['drug2cell'],
                                swap_axes=True,
                                dendrogram=False, 
                                n_genes=5,
                                values_to_plot="scores")

sc.write("……/sce_PCa.h5ad", adata_PCa)
#Drug filter
setwd("D:/Urinary system tumors/work/3_drug2cell")
PCa_drug_score <- read.csv("PCa_drug_score.csv",row.names = 1)
PCa_drug_targets <- read.csv("PCa_drug_targets.csv")
PCa_drug_score <- PCa_drug_score[which(PCa_drug_score$group == "Tumor cells"),]
PCa_drug_score <- PCa_drug_score[which(PCa_drug_score$pvals_adj < 0.05),]

PCa_drug_targets <- PCa_drug_targets[which(PCa_drug_targets$X %in% PCa_drug_score$names),]
###
x1 <- c()
for(i in 1:length(PCa_drug_targets$X)){
  temp <- cbind(PCa_drug_targets$X[i],
                strsplit(PCa_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}
PCa_network <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/PCa/TF_gene.txt",sep = "\t",header = T)

PCa_network$motif[which(PCa_network$motif == "ELF3_342")] <- "ELF3"
PCa_network$motif[which(PCa_network$motif == "TEAD1_796")] <- "TEAD1"
PCa_network$motif[which(PCa_network$motif == "KLF5_175")] <- "KLF5"
PCa_network$motif[which(PCa_network$motif == "FOS_137")] <- "FOS"
PCa_network$motif[which(PCa_network$motif == "FOSB_121")] <- "FOSB"
PCa_network$motif[which(PCa_network$motif == "CEBPB_140")] <- "CEBPB"
PCa_network$motif[which(PCa_network$motif == "ETS2_340")] <- "ETS2"
PCa_network$motif[which(PCa_network$motif == "KLF10_826")] <- "KLF10"

gene <- c(unique(PCa_network$geneName),unique(PCa_network$motif))
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- PCa_drug_score[PCa_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("PCa_drug.pdf",height = 5,width = 10)
write.table(x2, "./PCa_filter/PCa_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "./PCa_filter/PCa_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)

setwd("D:/Urinary system tumors/work/3_drug2cell")
ccRCC_drug_score <- read.csv("ccRCC_drug_score.csv",row.names = 1)
ccRCC_drug_targets <- read.csv("ccRCC_drug_targets.csv")
ccRCC_drug_score <- ccRCC_drug_score[which(ccRCC_drug_score$group == "Tumor cells"),]
ccRCC_drug_score <- ccRCC_drug_score[which(ccRCC_drug_score$pvals_adj < 0.05),]

ccRCC_drug_targets <- ccRCC_drug_targets[which(ccRCC_drug_targets$X %in% ccRCC_drug_score$names),]
###
x1 <- c()
for(i in 1:length(ccRCC_drug_targets$X)){
  temp <- cbind(ccRCC_drug_targets$X[i],
                strsplit(ccRCC_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}
ccRCC_network <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/ccRCC/TF_gene.txt",sep = "\t",header = T)

ccRCC_network$motif[which(ccRCC_network$motif == "FOSL1_142")] <- "FOSL1"
ccRCC_network$motif[which(ccRCC_network$motif == "RUNX1_733")] <- "RUNX1"
ccRCC_network$motif[which(ccRCC_network$motif == "NFIB_741")] <- "NFIB"
ccRCC_network$motif[which(ccRCC_network$motif == "NFIC_740")] <- "NFIC"
ccRCC_network$motif[which(ccRCC_network$motif == "CEBPB_140")] <- "CEBPB"
ccRCC_network$motif[which(ccRCC_network$motif == "HNF4A_662")] <- "HNF4A"
ccRCC_network$motif[which(ccRCC_network$motif == "ID4_75")] <- "ID4"

gene <- c(unique(ccRCC_network$geneName),unique(ccRCC_network$motif))
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- ccRCC_drug_score[ccRCC_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("ccRCC_drug.pdf",height = 5,width = 10)
write.table(x2, "./ccRCC_filter/ccRCC_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "./ccRCC_filter/ccRCC_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)

CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/filter/BC_cnv_gene.txt", header = T,sep = "\t")
CNV_gene <- na.omit(CNV_gene)
Time_gene <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/BC/2_BC_tumor_modulated_genes.txt", header = T,sep = "\t")

Time_CNV_gene <- Time_gene[unique(CNV_gene$gene),]
Time_CNV_gene <- Time_CNV_gene[order(Time_CNV_gene$morans_I,decreasing = T),]

setwd("D:/Urinary system tumors/work/3_drug2cell")
BC_drug_score <- read.csv("BC_drug_score.csv",row.names = 1)
BC_drug_targets <- read.csv("BC_drug_targets.csv")
BC_drug_score <- BC_drug_score[which(BC_drug_score$group == "Tumor cells"),]
BC_drug_score <- BC_drug_score[which(BC_drug_score$pvals_adj < 0.05),]

BC_drug_targets <- BC_drug_targets[which(BC_drug_targets$X %in% BC_drug_score$names),]
###
x1 <- c()
for(i in 1:length(BC_drug_targets$X)){
  temp <- cbind(BC_drug_targets$X[i],
                strsplit(BC_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}

gene <- Time_CNV_gene$gene_short_name
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- BC_drug_score[BC_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("BC_drug.pdf",height = 5,width = 10)
write.table(x2, "./BC_filter/BC_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "./BC_filter/BC_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)
#chrX  and chr Y
CNV_gene <- read.delim("clipboard",header = T)##table s2 BC_CNV_DEG_chrX_Y
CNV_gene <- na.omit(CNV_gene)
Time_gene <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/BC/2_BC_tumor_modulated_genes.txt", header = T,sep = "\t")

Time_CNV_gene <- Time_gene[unique(CNV_gene$gene),]
Time_CNV_gene <- Time_CNV_gene[order(Time_CNV_gene$morans_I,decreasing = T),]

BC_drug_score <- read.csv("D:/Urinary system tumors/work/3_drug2cell/BC_drug_score.csv",row.names = 1)
BC_drug_targets <- read.csv("D:/Urinary system tumors/work/3_drug2cell/BC_drug_targets.csv")
BC_drug_score <- BC_drug_score[which(BC_drug_score$group == "Tumor cells"),]
BC_drug_score <- BC_drug_score[which(BC_drug_score$pvals_adj < 0.05),]

BC_drug_targets <- BC_drug_targets[which(BC_drug_targets$X %in% BC_drug_score$names),]

x1 <- c()
for(i in 1:length(BC_drug_targets$X)){
  temp <- cbind(BC_drug_targets$X[i],
                strsplit(BC_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}

gene <- Time_CNV_gene$gene_short_name
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- BC_drug_score[BC_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("BC_drug.pdf",height = 5,width = 10)
write.table(x2, "BC_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "BC_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)
ccRCC_drug_score <- read.csv("D:/Urinary system tumors/work/3_drug2cell/ccRCC_drug_score.csv",row.names = 1)
ccRCC_drug_targets <- read.csv("D:/Urinary system tumors/work/3_drug2cell/ccRCC_drug_targets.csv")
ccRCC_drug_score <- ccRCC_drug_score[which(ccRCC_drug_score$group == "Tumor cells"),]
ccRCC_drug_score <- ccRCC_drug_score[which(ccRCC_drug_score$pvals_adj < 0.05),]

ccRCC_drug_targets <- ccRCC_drug_targets[which(ccRCC_drug_targets$X %in% ccRCC_drug_score$names),]
x1 <- c()
for(i in 1:length(ccRCC_drug_targets$X)){
  temp <- cbind(ccRCC_drug_targets$X[i],
                strsplit(ccRCC_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}
ccRCC_network <- read.table("./chrXY/ccRCC_table.txt",sep = "\t",header = T)


gene <- ccRCC_network$V1
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- ccRCC_drug_score[ccRCC_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("ccRCC_drug.pdf",height = 5,width = 10)
write.table(x2, "ccRCC_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "ccRCC_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)

PCa_drug_score <- read.csv("D:/Urinary system tumors/work/3_drug2cell/PCa_drug_score.csv",row.names = 1)
PCa_drug_targets <- read.csv("D:/Urinary system tumors/work/3_drug2cell/PCa_drug_targets.csv")
PCa_drug_score <- PCa_drug_score[which(PCa_drug_score$group == "Tumor cells"),]
PCa_drug_score <- PCa_drug_score[which(PCa_drug_score$pvals_adj < 0.05),]

PCa_drug_targets <- PCa_drug_targets[which(PCa_drug_targets$X %in% PCa_drug_score$names),]
###
x1 <- c()
for(i in 1:length(PCa_drug_targets$X)){
  temp <- cbind(PCa_drug_targets$X[i],
                strsplit(PCa_drug_targets$genes[i],",")[[1]])
  x1 <- rbind(x1,temp)
}
PCa_network <- read.table("./chrXY/PCa_table.txt",sep = "\t",header = T)

gene <- PCa_network$V1
x1 <- as.data.frame(x1)
x1 <- x1[which(x1$V2 %in% gene),]

x2 <- PCa_drug_score[PCa_drug_score$names %in% x1$V1,]
x3 <- x2[1:10,c("names","scores")]

x3$names <- factor(x3$names,levels = rev(x3$names) )
ggplot(x3,aes(names,scores))+
  geom_col(aes(fill=scores))+
  coord_flip()

ggsave("PCa_drug.pdf",height = 5,width = 10)
write.table(x2, "PCa_filter_drug_scores.txt",sep = "\t",quote = F,row.names = F)
write.table(x1, "PCa_filter_drug_targets.txt",sep = "\t",quote = F,row.names = F)
#Asgard ----
##Integrate normal samples and ecotype ----
load("D:/Urinary system tumors/work/5_Ecotype/ccRCC/ccRCC_Ecotypes.RData")
ccRCC_normal <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/1_ccRCC/2_ccRCC_normal_sce/4_ccRCC_normal_sce_after_cell_type.rds")

ccRCC_normal@meta.data$group <- "Normal"
ccRCC_E1@meta.data$group <- "E1"
ccRCC_E2@meta.data$group <- "E2"
ccRCC_E3@meta.data$group <- "E3"
ccRCC_E4@meta.data$group <- "E4"
ccRCC_E5@meta.data$group <- "E5"
ccRCC_E6@meta.data$group <- "E6"


####E1
SC.list <- list(
  ccRCC_normal = ccRCC_normal,
  ccRCC_E6 = ccRCC_E6
)
CellCycle = TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
  SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
  SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                                       nfeatures = anchor.features, verbose = FALSE)
}
SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
DefaultAssay(SC.integrated) <- "integrated"
if (CellCycle) {
  ##Cell Cycle Regression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
}
#else {
#  ##Run the standard workflow for visualization and clustering
#  SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
#  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
#}
##t-SNE and Clustering
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)


#Visualize alignment result
DimPlot(SC.integrated, reduction = "umap",group.by = "group")
#ccRCC_E1 <- SC.integrated
#ccRCC_E2 <- SC.integrated
#ccRCC_E3 <- SC.integrated
#ccRCC_E4 <- SC.integrated
#ccRCC_E5 <- SC.integrated
#ccRCC_E6 <- SC.integrated
rm(list=setdiff(ls(), c("ccRCC_E1","ccRCC_E2","ccRCC_E3",
                        "ccRCC_E4","ccRCC_E5","ccRCC_E6") ))

sce <- ccRCC_E1
DefaultAssay(sce) <- "RNA"
Gene.list <- list()
C_names <- NULL
for(i in unique(sce@meta.data$cell_type)){
  # i = unique(sce@meta.data$cell_type)[1]
  Idents(sce) <- "cell_type"
  c_cells <- subset(sce, cell_type == i)
  Idents(c_cells) <- "group"
  C_data <- FindMarkers(c_cells, ident.1 = Case, ident.2 = Control)
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
lapply(Gene.list, head)
save.image("D:/Urinary system tumors/work/6_drug_asgard/ccRCC.RData")

load("D:/Urinary system tumors/work/5_Ecotype/PCa/PCa_Ecotypes.RData")
PCa_normal <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/3_PCa/three_datasets/2_PCa_normal_sce/4_PCa_normal_sce_after_cell_type.rds")

PCa_normal@meta.data$group <- "Normal"

PCa_E1@meta.data$group <- "E1"
PCa_E2@meta.data$group <- "E2"
PCa_E3@meta.data$group <- "E3"
###E1
SC.list <- list(
  PCa_normal = PCa_normal,
  PCa_E1 = PCa_E1
)
CellCycle = TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
  SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
  SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                                       nfeatures = anchor.features, verbose = FALSE)
}
SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
DefaultAssay(SC.integrated) <- "integrated"
if (CellCycle) {
  ##Cell Cycle Regression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
}
#else {
  ##Run the standard workflow for visualization and clustering
#  SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
#  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
#}
##t-SNE and Clustering
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)


#Visualize alignment result
DimPlot(SC.integrated, reduction = "umap",group.by = "group")
#PCa_E1 <- SC.integrated
PCa_E1 <- SC.integrated


SC.list <- list(
  PCa_normal = PCa_normal,
  PCa_E3 = PCa_E3
)
CellCycle = TRUE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
  SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
  SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                                       nfeatures = anchor.features, verbose = FALSE)
}
SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
DefaultAssay(SC.integrated) <- "integrated"
if (CellCycle) {
  ##Cell Cycle Regression
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
  SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
}

##t-SNE and Clustering
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)


#Visualize alignment result
DimPlot(SC.integrated, reduction = "umap",group.by = "group")
#PCa_E3 <- SC.integrated
PCa_E3 <- SC.integrated
rm(list=setdiff(ls(), c("PCa_E3","PCa_E2","PCa_E1") ))
save.image("D:/Urinary system tumors/work/6_drug_asgard/PCa.RData")

##Run Asgard -----
library('Asgard')
PCa_Asgard <- list()
Genelist <- list()
#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="./DrugRef/prostate_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="./DrugRef/prostate_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = './DrugRef/prostate_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)


load("./Urinary system tumors/6_drug/PCa/PCa.RData")
sce <- PCa_E1
DefaultAssay(sce) <- "RNA"

E1_Asgard <- list()
C_names <- NULL
Gene.list <- list()
xxx <- unique(sce@meta.data[which(sce@meta.data$group == "E1"),]$cell_type)

for(i in xxx){
  # i = unique(sce@meta.data$cell_type)[1]
  Idents(sce) <- "cell_type"
  c_cells <- subset(sce, cell_type == i)
  Idents(c_cells) <- "group"
  C_data <- FindMarkers(c_cells, ident.1 = "E1", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
lapply(Gene.list, head)
#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")
########drug score

# Change the following two lines with the paths on your computer
gse92742_gctx_path <- "./reference/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "./reference/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"

cell_metadata <- sce@meta.data
sce@meta.data$cluster <- sce@meta.data$cell_type
sce@meta.data$celltype <- sce@meta.data$cell_type
sce@meta.data$sample<- sce@meta.data$group
Idents(sce)
Drug.score <- DrugScore(SC.integrated = sce,  
                        Gene.data = Gene.list,
                        Cell.type = names(Drug.ident.res),
                        Drug.data = Drug.ident.res,
                        FDA.drug.only = TRUE,
                        GSE92742.gctx = gse92742_gctx_path,
                        GSE70138.gctx = gse70138_gctx_path,
                        Case = "E1",
                        Tissue = "kidney"
)

#Select drug using drug socre
library(Hmisc)
Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)

E1_Asgard[["Gene.list"]] <- Gene.list
E1_Asgard[["Drug.ident.res"]] <- Drug.ident.res
E1_Asgard[["Drug.score"]] <- Drug.score
E1_Asgard[["Final.drugs"]] <- Final.drugs


sce <- PCa_E2
DefaultAssay(sce) <- "RNA"

E2_Asgard <- list()
C_names <- NULL
Gene.list <- list()
xxx <- unique(sce@meta.data[which(sce@meta.data$group == "E2"),]$cell_type)

for(i in xxx){
  # i = unique(sce@meta.data$cell_type)[1]
  Idents(sce) <- "cell_type"
  c_cells <- subset(sce, cell_type == i)
  Idents(c_cells) <- "group"
  C_data <- FindMarkers(c_cells, ident.1 = "E2", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
lapply(Gene.list, head)
#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")
########Drug score

# Change the following two lines with the paths on your computer
gse92742_gctx_path <- "./reference/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "./reference/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"

cell_metadata <- sce@meta.data
sce@meta.data$cluster <- sce@meta.data$cell_type
sce@meta.data$celltype <- sce@meta.data$cell_type
sce@meta.data$sample<- sce@meta.data$group
Idents(sce)
Drug.score <- DrugScore(SC.integrated = sce,  
                        Gene.data = Gene.list,
                        Cell.type = names(Drug.ident.res),
                        Drug.data = Drug.ident.res,
                        FDA.drug.only = TRUE,
                        GSE92742.gctx = gse92742_gctx_path,
                        GSE70138.gctx = gse70138_gctx_path,
                        Case = "E2",
                        Tissue = "kidney"
)

#Select drug using drug socre
library(Hmisc)
Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)

E2_Asgard[["Gene.list"]] <- Gene.list
E2_Asgard[["Drug.ident.res"]] <- Drug.ident.res
E2_Asgard[["Drug.score"]] <- Drug.score
E2_Asgard[["Final.drugs"]] <- Final.drugs



sce <- PCa_E3
DefaultAssay(sce) <- "RNA"

E3_Asgard <- list()
C_names <- NULL
Gene.list <- list()
xxx <- unique(sce@meta.data[which(sce@meta.data$group == "E3"),]$cell_type)

for(i in xxx){
  # i = unique(sce@meta.data$cell_type)[1]
  Idents(sce) <- "cell_type"
  c_cells <- subset(sce, cell_type == i)
  Idents(c_cells) <- "group"
  C_data <- FindMarkers(c_cells, ident.1 = "E3", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names
lapply(Gene.list, head)
#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                         drug.ref.profiles = drug.ref.profiles, 
                         repurposing.unit = "drug", 
                         connectivity = "negative", 
                         drug.type = "FDA")
########Drug score

# Change the following two lines with the paths on your computer
gse92742_gctx_path <- "./reference/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "./reference/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"

cell_metadata <- sce@meta.data
sce@meta.data$cluster <- sce@meta.data$cell_type
sce@meta.data$celltype <- sce@meta.data$cell_type
sce@meta.data$sample<- sce@meta.data$group
Idents(sce)
Drug.score <- DrugScore(SC.integrated = sce,  
                        Gene.data = Gene.list,
                        Cell.type = names(Drug.ident.res),
                        Drug.data = Drug.ident.res,
                        FDA.drug.only = TRUE,
                        GSE92742.gctx = gse92742_gctx_path,
                        GSE70138.gctx = gse70138_gctx_path,
                        Case = "E3",
                        Tissue = "kidney"
)

#Select drug using drug socre
library(Hmisc)
Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)

E3_Asgard[["Gene.list"]] <- Gene.list
E3_Asgard[["Drug.ident.res"]] <- Drug.ident.res
E3_Asgard[["Drug.score"]] <- Drug.score
E3_Asgard[["Final.drugs"]] <- Final.drugs


PCa_Asgard[["E1_Asgard"]] <- E1_Asgard
PCa_Asgard[["E2_Asgard"]] <- E2_Asgard
PCa_Asgard[["E3_Asgard"]] <- E3_Asgard
rm(list=setdiff(ls(), c("PCa_Asgard") ))

