#Epithelium monocle3 of scRNA------
library(Seurat)
#ccRCC
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")

ccRCC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor.rds")
ccRCC <-  subset(ccRCC, cell_type %in% c("Epithelium"))
ccRCC@meta.data$cell_type_second <- "Normal epithelium"
ccRCC@meta.data[which(ccRCC@meta.data$cell_type == "Epithelium"&ccRCC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"

sce_ccRCC_tumor@meta.data$cell_type_second <- NA
sce_ccRCC_tumor@meta.data[row.names(ccRCC@meta.data),]$cell_type_second <- as.character(ccRCC@meta.data$cell_type_second)

ccRCC <-  subset(sce_ccRCC_tumor, cell_type %in% "Epithelium")
rm("sce_ccRCC_tumor")

ccRCC <- NormalizeData(ccRCC) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
#ScaleData()
ccRCC <- RunPCA(ccRCC, features = VariableFeatures(object = ccRCC), verbose = F)

DimPlot(ccRCC, reduction = "pca")
DimHeatmap(ccRCC, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(ccRCC, ndims = 50)

pct <- ccRCC [["pca"]]@stdev / sum( ccRCC [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #43
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #13
pcs <- min(co1, co2)
pcs  #13
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))


library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()

pc_num= 1:13
ccRCC <- FindNeighbors(ccRCC, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 
ccRCC <- ccRCC %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)
DimPlot(ccRCC, label = T, reduction = "umap",group.by = "who")

ccRCC@meta.data$group <- ccRCC@meta.data$who
ccRCC@meta.data[which(ccRCC@meta.data$cell_type_second == "Normal epithelium"),]$group <- "Normal epithelium"

##CYTOTRACE
library(CytoTRACE2)
cytotrace2_result_sce <- cytotrace2(ccRCC, 
                                    is_seurat = T, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)
ggplot(cytotrace2_result_sce@meta.data,aes(x=group,y=CytoTRACE2_Relative))+geom_boxplot()

cytotrace2_result_sce@meta.data$CytoTRACE2_Score
cytotrace2_result_sce@meta.data$CytoTRACE2_Potency
library(monocle3)
data <- GetAssayData(cytotrace2_result_sce, assay = 'RNA', slot = 'counts')
cell_metadata <- cytotrace2_result_sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 15)

plot_pc_variance_explained(cds)


cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method="Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",resolution = 1e-06)
plot_cells(cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Potency",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")

#cds <- choose_cells(cds)
big_partition <- c(rep(1,length(cds@colData@rownames))) 
names(big_partition) <- cds@colData@rownames 
big_partition <- as.factor(big_partition)
cds@clusters$UMAP$partitions <- big_partition

cds <- learn_graph(cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)
#cds <- learn_graph(cds,use_partition = TRUE)

cds <- order_cells(cds,reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

ccRCC_epi <- cytotrace2_result_sce
ccRCC_epi_cds <- cds
rm(list=setdiff(ls(), c("ccRCC_epi","ccRCC_epi_cds") ))

shanji <- as.data.frame(pseudotime(ccRCC_epi_cds)[colnames(ccRCC_epi_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_epi_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=group,fill=group))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("#FFD2DB",#CD8Teff_CCL4
                             "#E17174",#CD8Teff_GZMK
                             "#B31620",#CD8Teff_IFIT2
                             "#01C1CC",#CD8Teff_MT1E
                             "#FFC981"#CD8Trm_ZNF683
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\2_ccRCC_pseutime_shanji.pdf"
       ,width = 5, height = 3)
dev.off()


modulated_genes <- graph_test(ccRCC_epi_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\2_ccRCC_tumor_modulated_genes.txt",quote = F,sep = "\t")
modulated_genes <- read.table("D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\2_ccRCC_tumor_modulated_genes.txt",header = T,sep = "\t")
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(ccRCC_epi_cds)[match(genes,rownames(rowData(ccRCC_epi_cds))),order(pseudotime(ccRCC_epi_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 4,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

rcl.list <- row_order(htkm) 

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(pt.matrix[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
})
reordergene <- c(clu_df[[1]]$GeneID,clu_df[[2]]$GeneID,clu_df[[3]]$GeneID,clu_df[[4]]$GeneID)

kmean_gene <- as.data.frame(do.call(rbind, clu_df))
kmean_gene_order <- kmean_gene[order(factor(kmean_gene$Cluster,levels = c("cluster1","cluster2","cluster3","cluster4"))),]
kmean_gene_order$Cluster <- factor(kmean_gene_order$Cluster, levels = c("cluster1","cluster2","cluster3","cluster4"))
pt.matrix_order <- pt.matrix[order(factor(rownames(pt.matrix),levels =(kmean_gene_order$GeneID))),]
htkm_order <- Heatmap(
  pt.matrix_order,
  name = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names  = FALSE,
  show_column_names = FALSE,
  row_names_gp  = gpar(fontsize = 6),
  #split = split,
  #km = 1,
  #row_title_rot = 0,
  cluster_rows  = T,
  row_order = kmean_gene_order$GeneID,
  row_split = kmean_gene_order$Cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices  = FALSE,
  cluster_columns = FALSE)
print(htkm_order)#2_ccRCC_Pseutime_gene 4x6


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
##cluster1

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster1"),1]

gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                toType = "ENTREZID" , 
                OrgDb = org.Hs.eg.db)
head(gene.df)

enrich.go <- enrichGO(gene = gene.df$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',  
                      ont = 'ALL',  
                      pAdjustMethod = 'BH',  
                      pvalueCutoff = 0.05,  
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

ego <- as.data.frame(enrich.go)
ego$cluster <- "cluster1"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\cluster1.txt",sep = "\t",quote = F,row.names = F)
library(dplyr)
ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)

for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:30),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#6x8

saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\clusterdataframe.rds")


plot_cells(ccRCC_epi_cds, color_cells_by = "group",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+
  scale_color_manual(values=c("#FFD2DB",
                              "#E17174",
                              "#B31620",
                              "#01C1CC",
                              "#FFC981"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\2_ccRCC_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_epi_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\2_ccRCC_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

ccRCC_epi@meta.data$group
ccRCC_epi@meta.data$CytoTRACE2_Relative
ggplot(ccRCC_epi@meta.data, aes(x=group,fill=group,y=CytoTRACE2_Relative)) + geom_boxplot()

ccRCC_epi@meta.data$X <- factor(ccRCC_epi@meta.data$group,levels = c("Normal epithelium","I","II","III"))
p1 <- ggplot(ccRCC_epi@meta.data, aes(x=X,fill=X,y=CytoTRACE2_Relative))+ 
  geom_boxplot(width=0.3,position=position_dodge(width=0.1))+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  theme(axis.text.x=element_blank())+ 
  xlab("")+ 
  ylab("")

#p1=p1+ geom_point(data = tmp, size = 0.6, shape = 21, position = position_jitterdodge()) 
p1=p1+theme(strip.text = element_text(colour = 'black'))+ theme(legend.position="none")
p1=p1+stat_summary(fun.y=median,linetype="dotted", geom="smooth",size=0.8, aes(group=0),color='gray30')+
  scale_fill_manual(values=c("#01C1CC","#FFD2DB","#E17174","#B31620"))
p1#ccRCC_cytotrace2 4x3


library(Seurat)

#BC
BC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/2_BC/2_BC_after_ident_tumor.rds")
BC <-  subset(BC, cell_type %in% c("Epithelium"))
BC@meta.data$cell_type_second <- "Normal epithelium"
BC@meta.data[which(BC@meta.data$cell_type == "Epithelium"&BC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
BC <-  subset(BC, cell_type %in% "Epithelium")
BC <- NormalizeData(BC) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
BC <- RunPCA(BC, features = VariableFeatures(object = BC), verbose = F)
DimPlot(BC, reduction = "pca")
DimHeatmap(BC, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(BC, ndims = 50)


pct <- BC [["pca"]]@stdev / sum( BC [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #43
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #13
pcs <- min(co1, co2)
pcs  #13
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))


library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()

pc_num= 1:12
BC <- FindNeighbors(BC, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 
BC <- BC %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)
DimPlot(BC, label = T, reduction = "umap",group.by = "cell_type_second")

##CYTOTRACE
library(CytoTRACE2)
cytotrace2_result_sce <- cytotrace2(BC, 
                                    is_seurat = T, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)
cytotrace2_result_sce@meta.data$CytoTRACE2_Score
cytotrace2_result_sce@meta.data$CytoTRACE2_Potency
BC_epi <- cytotrace2_result_sce
ggplot(BC_epi@meta.data,aes(x=cell_type_second,y=CytoTRACE2_Relative))+geom_boxplot()

####
library(monocle3)
data <- GetAssayData(BC_epi, assay = 'RNA', slot = 'counts')
cell_metadata <- cytotrace2_result_sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 12)

plot_pc_variance_explained(cds)


cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method="Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",resolution = 1e-06)

plot_cells(cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "cell_type_second",cell_size = 0.5)+theme(legend.position = "right")

#cds <- choose_cells(cds)
big_partition <- c(rep(1,length(cds@colData@rownames))) 
names(big_partition) <- cds@colData@rownames 
big_partition <- as.factor(big_partition)
cds@clusters$UMAP$partitions <- big_partition

cds <- learn_graph(cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)
cds <- learn_graph(cds,use_partition = TRUE)

cds <- order_cells(cds,reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
plot_cells(cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "RNA_snn_res.0.6",cell_size = 0.5)+theme(legend.position = "right")
ggplot(cytotrace2_result_sce@meta.data,aes(x=RNA_snn_res.0.6,y=CytoTRACE2_Relative))+geom_boxplot()
DimPlot(cytotrace2_result_sce,group.by = "RNA_snn_res.0.6",label = T,label.size = 5)
BC_epi_cds <- cds
rm(list=setdiff(ls(), c("BC_epi","BC_epi_cds") ))

shanji <- as.data.frame(pseudotime(BC_epi_cds)[colnames(BC_epi_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_epi_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=cell_type_second,fill=cell_type_second))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("#01C1CC",
                             "#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\2_BC_pseutime_shanji.pdf"
       ,width = 5, height = 3)
dev.off()


modulated_genes <- graph_test(BC_epi_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\2_BC_tumor_modulated_genes.txt",quote = F,sep = "\t")
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(BC_epi_cds)[match(genes,rownames(rowData(BC_epi_cds))),order(pseudotime(BC_epi_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 4,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

rcl.list <- row_order(htkm) 

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(pt.matrix[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
})
reordergene <- c(clu_df[[1]]$GeneID,clu_df[[2]]$GeneID,clu_df[[3]]$GeneID,clu_df[[4]]$GeneID)

kmean_gene <- as.data.frame(do.call(rbind, clu_df))
kmean_gene_order <- kmean_gene[order(factor(kmean_gene$Cluster,levels = c("cluster1","cluster2","cluster3","cluster4"))),]
kmean_gene_order$Cluster <- factor(kmean_gene_order$Cluster, levels = c("cluster1","cluster2","cluster3","cluster4"))
pt.matrix_order <- pt.matrix[order(factor(rownames(pt.matrix),levels =(kmean_gene_order$GeneID))),]
htkm_order <- Heatmap(
  pt.matrix_order,
  name = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names  = FALSE,
  show_column_names = FALSE,
  row_names_gp  = gpar(fontsize = 6),
  #split = split,
  #km = 1,
  #row_title_rot = 0,
  cluster_rows  = T,
  row_order = kmean_gene_order$GeneID,
  row_split = kmean_gene_order$Cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices  = FALSE,
  cluster_columns = FALSE)
print(htkm_order)#2_BC_Pseutime_gene 4x6

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
##cluster1

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster1"),1]

gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                toType = "ENTREZID" , 
                OrgDb = org.Hs.eg.db)
head(gene.df)

enrich.go <- enrichGO(gene = gene.df$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',  
                      ont = 'ALL',  
                      pAdjustMethod = 'BH',  
                      pvalueCutoff = 0.05,  
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

ego <- as.data.frame(enrich.go)
ego$cluster <- "cluster1"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\cluster1.txt",sep = "\t",quote = F,row.names = F)

library(dplyr)


ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:30),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#6x8

saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\clusterdataframe.rds")

plot_cells(BC_epi_cds, color_cells_by = "cell_type_second",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("#01C1CC",
                                                              "#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\2_BC_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_epi_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\BC\\2_BC_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()


BC_epi@meta.data$X <- factor(BC_epi@meta.data$group,levels = c("Normal epithelium","I","II","III"))
p1 <- ggplot(BC_epi@meta.data, aes(x=cell_type_second,fill=cell_type_second,y=CytoTRACE2_Relative))+ 
  geom_boxplot(width=0.3,position=position_dodge(width=0.1))+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  theme(axis.text.x=element_blank())+ 
  xlab("")+ 
  ylab("")

#p1=p1+ geom_point(data = tmp, size = 0.6, shape = 21, position = position_jitterdodge()) 
p1=p1+theme(strip.text = element_text(colour = 'black'))+ theme(legend.position="none")
p1=p1+stat_summary(fun.y=median,linetype="dotted", geom="smooth",size=0.8, aes(group=0),color='gray30')+
  scale_fill_manual(values=c("#01C1CC",
                             "#377DB8"))
p1#BC_cytotrace2 4x3

sce_PCa_tumor <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/3_PCa/three_datasets/1_PCa_tumor_sce/4_PCa_tumor_sce_after_cell_type.rds")
PCa <-  subset(PCa, cell_type %in% c("Epithelium"))
PCa@meta.data$cell_type_second <- "Normal epithelium"
PCa@meta.data[which(PCa@meta.data$cell_type == "Epithelium"&PCa@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"

sce_PCa_tumor@meta.data$cell_type_second <- NA
sce_PCa_tumor@meta.data[row.names(PCa@meta.data),]$cell_type_second <- as.character(PCa@meta.data$cell_type_second)

PCa <-  subset(sce_PCa_tumor, cell_type %in% "Epithelium")
rm("sce_PCa_tumor")

PCa <- NormalizeData(PCa) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
#ScaleData()
PCa <- RunPCA(PCa, features = VariableFeatures(object = PCa), verbose = F)

DimPlot(PCa, reduction = "pca")
DimHeatmap(PCa, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(PCa, ndims = 50)


pct <- PCa [["pca"]]@stdev / sum( PCa [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
；
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #43

co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #13
pcs <- min(co1, co2)
pcs  #13
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))


library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()

pc_num= 1:13
PCa <- FindNeighbors(PCa, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 
PCa <- PCa %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)
DimPlot(PCa, label = T, reduction = "umap",group.by = "gleason")

##CYTOTRACE
library(CytoTRACE2)
PCa <- cytotrace2(PCa,
                  is_seurat = T,
                  slot_type = "counts", 
                  species = 'human',
                  seed = 1234)
PCa@meta.data$group <- PCa@meta.data$gleason
PCa@meta.data[which(PCa@meta.data$gleason %in% c("3+3","3+4")),]$group <- "G3"
PCa@meta.data[which(PCa@meta.data$gleason %in% c("4+3","4+5")),]$group <- "G4"
PCa@meta.data[which(PCa@meta.data$cell_type_second == "Normal epithelium"),]$group <- "Normal epithelium"
ggplot(PCa@meta.data,aes(x=group,y=CytoTRACE2_Relative))+geom_boxplot()

####
library(monocle3)
data <- GetAssayData(PCa, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 13)

plot_pc_variance_explained(cds)


cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method="Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",resolution = 1e-06)

plot_cells(cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "GSM_name",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Potency",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Score",cell_size = 0.5)+theme(legend.position = "right")

#cds <- choose_cells(cds)
big_partition <- c(rep(1,length(cds@colData@rownames))) 
names(big_partition) <- cds@colData@rownames 
big_partition <- as.factor(big_partition)
cds@clusters$UMAP$partitions <- big_partition

cds <- learn_graph(cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)
#cds <- learn_graph(cds,use_partition = TRUE)

cds <- order_cells(cds,reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
plot_cells(cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")

PCa_epi <- PCa
PCa_epi_cds <- cds
rm(list=setdiff(ls(), c("PCa_epi","PCa_epi_cds") ))

PCa <-  subset(PCa_epi, group %in% c("Normal epithelium","G3","G4")) 
PCa@meta.data$group2 <- PCa@meta.data$gleason
PCa@meta.data[which(PCa@meta.data$cell_type_second == "Normal epithelium"),]$group2 <- "Normal epithelium"
####
library(monocle3)
data <- GetAssayData(PCa, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 13)

plot_pc_variance_explained(cds)


cds <- align_cds(cds, alignment_group = "orig.ident", preprocess_method = "PCA")
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method="Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",resolution = 1e-06)

plot_cells(cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")
#plot_cells(cds, color_cells_by = "GSM_name",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Potency",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(cds, color_cells_by = "CytoTRACE2_Score",cell_size = 0.5)+theme(legend.position = "right")

#cds <- choose_cells(cds)
big_partition <- c(rep(1,length(cds@colData@rownames))) 
names(big_partition) <- cds@colData@rownames 
big_partition <- as.factor(big_partition)
cds@clusters$UMAP$partitions <- big_partition

cds <- learn_graph(cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)
#cds <- learn_graph(cds,use_partition = TRUE)

cds <- order_cells(cds,reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
plot_cells(cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")

plot_cells(PCa_epi_cds, color_cells_by = "group",cell_size = 0.5)+theme(legend.position = "right")

PCa_epi_cds <- cds
rm(list=setdiff(ls(), c("PCa_epi","PCa_epi_cds") ))

shanji <- as.data.frame(pseudotime(PCa_epi_cds)[colnames(PCa_epi_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_epi_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime


library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=group2,fill=group2))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("#DAFFDC",#
"#1E8F8C",
"#105958",#
"#073535",#
"#01C1CC"#
))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\2_PCa_pseutime_shanji.pdf"
       ,width = 5, height = 3)
dev.off()


modulated_genes <- graph_test(PCa_epi_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\2_PCa_tumor_modulated_genes.txt",quote = F,sep = "\t")
modulated_genes <- read.table("D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\2_PCa_tumor_modulated_genes.txt",header = T,sep = "\t")
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(PCa_epi_cds)[match(genes,rownames(rowData(PCa_epi_cds))),order(pseudotime(PCa_epi_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 4,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

rcl.list <- row_order(htkm) 

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(pt.matrix[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
})
reordergene <- c(clu_df[[1]]$GeneID,clu_df[[2]]$GeneID,clu_df[[3]]$GeneID,clu_df[[4]]$GeneID)

kmean_gene <- as.data.frame(do.call(rbind, clu_df))
kmean_gene_order <- kmean_gene[order(factor(kmean_gene$Cluster,levels = c("cluster1","cluster2","cluster3","cluster4"))),]
kmean_gene_order$Cluster <- factor(kmean_gene_order$Cluster, levels = c("cluster1","cluster2","cluster3","cluster4"))
pt.matrix_order <- pt.matrix[order(factor(rownames(pt.matrix),levels =(kmean_gene_order$GeneID))),]
htkm_order <- Heatmap(
  pt.matrix_order,
  name = "z-score",
  col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names  = FALSE,
  show_column_names = FALSE,
  row_names_gp  = gpar(fontsize = 6),
  #split = split,
  #km = 1,
  #row_title_rot = 0,
  cluster_rows  = T,
  row_order = kmean_gene_order$GeneID,
  row_split = kmean_gene_order$Cluster,
  cluster_row_slices = FALSE,
  cluster_column_slices  = FALSE,
  cluster_columns = FALSE)
print(htkm_order)#2_PCa_Pseutime_gene 4x6

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
##cluster1

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster1"),1]

gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                toType = "ENTREZID" , 
                OrgDb = org.Hs.eg.db)
head(gene.df)

enrich.go <- enrichGO(gene = gene.df$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',  
                      ont = 'ALL',  
                      pAdjustMethod = 'BH',  
                      pvalueCutoff = 0.05,  
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

ego <- as.data.frame(enrich.go)
ego$cluster <- "cluster1"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\cluster1.txt",sep = "\t",quote = F,row.names = F)

library(dplyr)


ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:30),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#6x8

saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\clusterdataframe.rds")
plot_cells(PCa_epi_cds, color_cells_by = "group2",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("#DAFFDC",#
                                                              "#1E8F8C",
                                                              "#105958",#
                                                              "#073535",#
                                                              "#01C1CC"#
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\2_PCa_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_epi_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\2_PCa_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

PCa <-  subset(PCa_epi, group %in% c("Normal epithelium","G3","G4"))
PCa@meta.data$group2 <- PCa@meta.data$gleason
PCa@meta.data[which(PCa@meta.data$cell_type_second == "Normal epithelium"),]$group2 <- "Normal epithelium"

PCa@meta.data$group2 <- factor(PCa@meta.data$group2,levels = c("Normal epithelium","3+3","3+4","4+3","4+5"))
p1 <- ggplot(PCa@meta.data, aes(x=group2,fill=group2,y=CytoTRACE2_Relative))+ 
  geom_boxplot(width=0.3,position=position_dodge(width=0.1))+ 
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  theme(axis.text.x=element_blank())+ 
  xlab("")+ 
  ylab("")
#p1=p1+ geom_point(data = tmp, size = 0.6, shape = 21, position = position_jitterdodge()) 
p1=p1+theme(strip.text = element_text(colour = 'black'))+ theme(legend.position="none")
p1=p1+stat_summary(fun.y=median,linetype="dotted", geom="smooth",size=0.8, aes(group=0),color='gray30')+
  scale_fill_manual(values=c("#01C1CC","#DAFFDC",#
                             "#1E8F8C",
                             "#105958",#
                             "#073535"))
p1#PCa_cytotrace2 4x3
##GO ----
#Read the seurat object of epithelial cells
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_ccRCC_tumor_epi.RData")
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_PCa_tumor_epi.RData")
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_BC_tumor_epi.RData")
Idents(ccRCC_epi) <- "cell_type_second"
ccRCC_marker <- FindAllMarkers(ccRCC_epi,only.pos = T)
ccRCC_Tumor_gene <- ccRCC_marker[which(ccRCC_marker$cluster == "Tumor cells"),"gene"]
gene_entrez <- bitr(
  ccRCC_Tumor_gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",    
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE    
)
ccRCC_Tumor_ego <- ego
ccRCC_Tumor <- ccRCC_Tumor_ego@result[which(ccRCC_Tumor_ego@result$p.adjust <= 0.05),]
write.table(ccRCC_Tumor_ego@result[which(ccRCC_Tumor_ego@result$p.adjust <= 0.05),],"ccRCC_Tumor_go.txt",quote = F,sep = "\t")


PCa_epi <- subset(x = PCa_epi, group != "none")

Idents(PCa_epi) <- "cell_type_second"
PCa_marker <- FindAllMarkers(PCa_epi,only.pos = T)

PCa_Tumor_gene <- PCa_marker[which(PCa_marker$cluster == "Tumor cells"),"gene"]
gene_entrez <- bitr(
  PCa_Tumor_gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",    
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
PCa_Tumor_ego <- ego
PCa_Tumor <- PCa_Tumor_ego@result[which(PCa_Tumor_ego@result$p.adjust <= 0.05),]
write.table(PCa_Tumor_ego@result[which(PCa_Tumor_ego@result$p.adjust <= 0.05),],"PCa_Tumor_go.txt",quote = F,sep = "\t")

Idents(BC_epi) <- "cell_type_second"
BC_marker <- FindAllMarkers(BC_epi,only.pos = T)
BC_Tumor_gene <- BC_marker[which(BC_marker$cluster == "Tumor cells"),"gene"]
gene_entrez <- bitr(
  BC_Tumor_gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",     
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
BC_Tumor_ego <- ego
BC_Tumor <- BC_Tumor_ego@result[which(BC_Tumor_ego@result$p.adjust <= 0.05),]
write.table(BC_Tumor_ego@result[which(BC_Tumor_ego@result$p.adjust <= 0.05),],"BC_Tumor_go.txt",quote = F,sep = "\t")
com_item <- intersect(intersect(BC_Tumor$Description, PCa_Tumor$Description), ccRCC_Tumor$Description)
com_item <- c("cell-cell junction organization",
              "gland development",
              "response to peptide hormone",
              "cell-cell junction assembly",
              "renal system development",
              "regulation of apoptotic signaling pathway",
              "negative regulation of apoptotic signaling pathway",
              "response to hypoxia")
dotplot(BC_Tumor_ego, showCategory =c(com_item,"epithelial cell development","cell growth"))
dotplot(PCa_Tumor_ego, showCategory =c(com_item,"epithelial cell development","cell growth"))
dotplot(ccRCC_Tumor_ego, showCategory =c(com_item,"negative regulation of growth"))
##NMF -----
pbmc_rna <- readRDS("./NMF/PCa/BC_sce_anno_second.rds")
BC_sce <- readRDS("./ALL_data/BC_sce_anno_second.rds")
pbmc_rna <-subset(BC_sce, cell_type %in% c("B_cells","Epithelium"))
ccRCC_sce <- readRDS("./ALL_data/ccRCC_sce_anno_second2.rds")
pbmc_rna <- subset(ccRCC_sce,subset = cell_type %in%  c("B_cells","Epithelium") )
PCa_sce <- readRDS("./ALL_data/PCa_sce_anno_second.rds")
pbmc_rna <- subset(PCa_sce,subset = cell_type %in%  c("B_cells","Epithelium") )

gene_names <- rownames(pbmc_rna)
mt_genes <- grep("^MT-", gene_names, value = TRUE)
ribo_genes <- grep("^RPL|^RPS", gene_names, value = TRUE)
genes_to_remove <- unique(c(mt_genes, ribo_genes))
genes_to_keep <- setdiff(gene_names, genes_to_remove)
pbmc_rna <- subset(pbmc_rna, features = genes_to_keep)
pbmc_rna@meta.data$sample_cell <- paste(pbmc_rna@meta.data$orig.ident,pbmc_rna@meta.data$cell_type,sep = "_")
table_cell <- table(pbmc_rna@meta.data$sample_cell)
table_cell_50 <- names(table_cell)[as.numeric(table_cell) > 100]
Seurat_obj <- subset(pbmc_rna,subset = sample_cell %in% table_cell_50 )
Target_seurat <- Seurat_obj
DefaultAssay(Target_seurat) <- "RNA"
table(Target_seurat@meta.data$sample_cell)

run_nmf_extract <- function(sample_id, range = 2:10, n_cores = 80, output_dir = "./NMF/PCa/nmf_rds/") {
  if (!sample_id %in% Target_seurat@meta.data$sample_cell) {
    warning(paste("Sample", sample_id, "not found, skipping..."))
    return(NULL)
  }
  
  temp_cells <- rownames(Target_seurat@meta.data)[Target_seurat@meta.data$sample_cell == sample_id]
  temp_seurat <- subset(Target_seurat, cells = temp_cells)
  
  temp_seurat <- SCTransform(temp_seurat, return.only.var.genes = FALSE, verbose = FALSE)
  
  scaled_data <- as.matrix(GetAssayData(temp_seurat, assay = 'SCT', slot = 'scale.data'))
  scaled_data <- scaled_data[VariableFeatures(temp_seurat), ]
  
  scaled_data[scaled_data < 0] <- 0
  scaled_data <- scaled_data[matrixStats::rowVars(scaled_data) > 0, ]
  if (nrow(scaled_data) < 10) {
    warning(paste("Sample", sample_id, "has too few variable genes, skipping..."))
    return(NULL)
  }
  
  set.seed(123)
  res.list <- parallel::mclapply(range, function(r) {
    NMF::nmf(scaled_data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
  }, mc.cores = n_cores)
  names(res.list) <- paste0("Rank_", range)
  
  saveRDS(res.list,file = paste0(output_dir, sample_id, ".RDS"))
}

sample_list <- unique(Target_seurat@meta.data$sample_cell)
library(pbapply)
result_list <- pblapply(sample_list, function(sample_id) {
  tryCatch({
    cat("\nProcessing:", sample_id)
    run_nmf_extract(sample_id)
  }, error = function(e) {
    message(paste("\nError in", sample_id, ":", e$message))
    return(NULL)
  })
})

{
  library(fs)
  library(stringr)
  
  input_dir <- "./NMF/PCa/nmf_rds/"
  setwd(input_dir)
  
  files <- list.files(pattern = "\\.RDS$")
  
  for (file in files) {
    celltype <- str_extract(file, "(?<=_)[^_]+(?=\\.RDS)")
    dir_create(celltype)
    file_move(file, path(celltype, file))
  }
}

{
  combine_nmf_components <- function(res_list, sample_name = "sample") {
    combined_matrix <- matrix(nrow = 3000, ncol = 0)
    col_names <- character(0)
    
    base_tag <- paste0(sample_name, "_rank4_9_nruns10.RDS")
    
    for (rank_idx in seq_along(res_list)) {
      current_rank <- 1 + rank_idx
      W_matrix <- basis(res_list[[rank_idx]])
      
      if (nrow(W_matrix) != 3000) {
        message("skip ", sample_name, " Rank ", current_rank, "：num row ", nrow(W_matrix), " != 3000。")
        next
      }
      
      new_colnames <- paste(base_tag, current_rank, seq_len(ncol(W_matrix)), sep = ".")
      combined_matrix <- cbind(combined_matrix, W_matrix)
      col_names <- c(col_names, new_colnames)
    }
    
    colnames(combined_matrix) <- col_names
    return(combined_matrix)
  }
  
  setwd("./NMF/PCa/nmf_rds/")
  dir <- dir()
  
  cell_types <- "Epithelium"
  
  dir.create("./NMF/PCa/nmf_rds/result")
  dir.create("./NMF/PCa/nmf_rds/result/temp")
  
  for (cell_type in cell_types) {
    docu_celltype <- paste0("./NMF/PCa/nmf_rds/",cell_type)
    setwd(docu_celltype)
    Genes_nmf_w_basis <- list()
    docu <- dir()
    
    for (i in seq_along(docu)) {
      print(docu[i])
      temp_rds <- readRDS(docu[i])
      final_matrix <- combine_nmf_components(temp_rds, sample_name = docu[i])
      Genes_nmf_w_basis[[i]] <- final_matrix
    }
    names(Genes_nmf_w_basis) <- docu
    
    Genes_nmf_w_basis <- Genes_nmf_w_basis[sapply(Genes_nmf_w_basis, function(mat) ncol(mat) > 0)]
    
    saveRDS(Genes_nmf_w_basis, paste0("./NMF/PCa/nmf_rds/result/temp/", cell_type,"_Genes_nmf_w_basis.rds"))
  }
  
  dir.create("./NMF/PCa/nmf_rds/result/cell_basis")
  setwd("./NMF/PCa/nmf_rds/result/temp/")
  file <- dir("./NMF/PCa/nmf_rds/result/temp/")
  cell_file <-"./NMF/PCa/nmf_rds/result/cell_basis/"
  
  file_move(file, cell_file)
  setwd("./NMF/PCa/nmf_rds/result/cell_basis/")
  dir <- dir()
  
  cell_types <- "Epithelium"
  
  dir.create("./NMF/PCa/New_resdullt")
  
  for (i in 1:length(dir)) {
    cell_type <- strsplit(dir[i], split = "_")[[1]][1]
    
    if (cell_type == "Epithelium") {
      print(paste("Epithelial:", dir[i]))
      
      Genes_nmf_w_basis <- readRDS(dir[i])
      
      intra_min_parameter <- 25
      intra_max_parameter <- 10
      inter_min_parameter <- 10
      source("./NMF/robust_nmf_programs.R")
      source("./NMF/custom_magma.R")
      
      nmf_programs <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
      nmf_programs <- lapply(nmf_programs,toupper)
      
      nmf_filter_ccle <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)
      nmf_programs <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
      nmf_programs <- do.call(cbind, nmf_programs)
      
      nmf_intersect <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y))))
      
      nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect), method="average")
      nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
      nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]
      
      Min_intersect_initial <- 10
      Min_intersect_cluster <- 10
      Min_group_size <- 3
      
      Sorted_intersection <- sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)) , decreasing = TRUE)
      
      Cluster_list <- list()
      MP_list <- list()
      k <- 1
      Curr_cluster <- c()
      
      nmf_intersect_original <- nmf_intersect
      
      while (Sorted_intersection[1]>Min_group_size) {
        Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
        
        Genes_MP <- nmf_programs[,names(Sorted_intersection[1])]
        nmf_programs <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]
        Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
        NMF_history <- Genes_MP
        
        while (Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {
          Curr_cluster <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
          
          Genes_MP_temp <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)
          Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]
          
          if (length(Genes_at_border)>1){
            Genes_curr_NMF_score <- c()
            for (j in Curr_cluster) {
              curr_study <- paste(strsplit(j , "[.]")[[1]][1],"RDS", sep  = ".")
              Q <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))] ,j]
              names(Q) <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])
              Genes_curr_NMF_score <- c(Genes_curr_NMF_score, Q)
            }
            Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
            Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]
            
            Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
          } else {
            Genes_MP_temp <- names(Genes_MP_temp)[1:50]
          }
          
          NMF_history <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])
          Genes_MP <- Genes_MP_temp[1:50]
          
          nmf_programs <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]
          
          Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
        }
        
        Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
        MP_list[[paste0("MP_",k)]] <- Genes_MP
        k <- k+1
        
        nmf_intersect <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect)) , -match(Curr_cluster,colnames(nmf_intersect))]
        
        Sorted_intersection <- sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)) , decreasing = TRUE)
        
        Curr_cluster <- c()
        print(dim(nmf_intersect)[2])
      }
      
      inds_sorted <- c()
      for (j in 1:length(Cluster_list)){
        inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
      }
      inds_new <- c(inds_sorted , which(is.na(match(1:dim(nmf_intersect_original)[2],inds_sorted))))
      data <- nmf_intersect_original[inds_new,rev(inds_new)]
      
      nmf_intersect_meltI_NEW <- reshape2::melt(data)
      
      df <- bind_rows(
        lapply(names(Cluster_list), function(cluster_name) {
          data.frame(Cluster = cluster_name, Content = Cluster_list[[cluster_name]], stringsAsFactors = FALSE)
        })
      )
      
      anno_data <- df
      nmf_intersect_meltI_NEW2 <- nmf_intersect_meltI_NEW
      nmf_intersect_meltI_NEW2$Var1_cluster <- anno_data$Cluster[match(nmf_intersect_meltI_NEW2$Var1,anno_data$Content)]
      nmf_intersect_meltI_NEW2$Var2_cluster <- anno_data$Cluster[match(nmf_intersect_meltI_NEW2$Var2,anno_data$Content)]
      nmf_intersect_meltI_NEW2 <- na.omit(nmf_intersect_meltI_NEW2)
      
      top_annotation_data <- nmf_intersect_meltI_NEW2 %>% distinct(Var1, Var1_cluster) %>% mutate(y_pos = 1)
      bot_annotation_data <- nmf_intersect_meltI_NEW2 %>% distinct(Var2, Var2_cluster) %>% mutate(y_pos = 1)
      
      cluster_levels <- unique(anno_data$Cluster)
      cluster_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(cluster_levels)), cluster_levels)
      
      var1_levels <- unique(nmf_intersect_meltI_NEW2$Var1)
      var2_levels <- unique(nmf_intersect_meltI_NEW2$Var2)
      
      top_annotation_data <- nmf_intersect_meltI_NEW2 %>% distinct(Var1, Var1_cluster) %>% mutate(Var1 = factor(Var1, levels = var1_levels)) %>% arrange(Var1) %>% mutate(y_pos = 1)
      left_annotation_data <- nmf_intersect_meltI_NEW2 %>% distinct(Var2, Var2_cluster) %>% mutate(Var2 = factor(Var2, levels = var2_levels)) %>% arrange(Var2) %>% mutate(x_pos = 1)
      
      cluster_levels <- unique(c(top_annotation_data$Var1_cluster, left_annotation_data$Var2_cluster))
      library(ggsci)
      npg_colors <- pal_npg()(length(cluster_levels))
      cluster_colors <- setNames(npg_colors, cluster_levels)
      
      top_annotation <- ggplot(top_annotation_data, aes(x = Var1, y = y_pos, fill = Var1_cluster)) + geom_tile() + scale_fill_manual(values = cluster_colors) + theme_void() + theme(legend.position = "none")
      left_annotation <- ggplot(left_annotation_data, aes(x = x_pos, y = Var2, fill = Var2_cluster)) + geom_tile() + scale_fill_manual(values = cluster_colors) + theme_void() + theme(legend.position = "none")
      
      library(scales)
      p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + geom_tile() + scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111], mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") + scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111], mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") + theme(axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
      
      p2 <- ggplot(data = nmf_intersect_meltI_NEW2, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + geom_tile() + scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111], mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") + scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111], mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") + theme(axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
      
      library(aplot)
      final_plot <- p2 %>% insert_top(top_annotation, height = 0.05) %>% insert_left(left_annotation, width = 0.05)
      
      pdf(paste0("./NMF/PCa/New_resdullt/",cell_type,"_NMF.pdf"),width = 7,height = 4)
      print(final_plot)
      dev.off()
      
      pdf(paste0("./NMF/PCa/New_resdullt/",cell_type,"all_NMF.pdf"),width = 7,height = 4)
      print(p1)
      dev.off()
      
      saveRDS(MP_list, paste0("./NMF/PCa/New_resdullt/",cell_type,"_MP_list.rds"))
    } else {
      print(paste("skip non-Epithelial:", cell_type))
    }
  }
}

#Pseutime analysis of archR -----
##PCa scATAC------
proj_epi <- subsetArchRProject(
  ArchRProj = PCa_ArchR_all_proj,
  cells = cell,
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  force = TRUE
)
dim <- 30
IterativeLSI <- "IterativeLSI30"
proj_epi <- addIterativeLSI(
  ArchRProj = proj_epi,
  useMatrix = "TileMatrix", 
  name = IterativeLSI, 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)


proj_epi <- addClusters(
  input = proj_epi,
  reducedDims = IterativeLSI,
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)
proj_epi <- addUMAP(
  ArchRProj = proj_epi, 
  reducedDims = IterativeLSI, 
  name = "UMAP30", 
  nNeighbors = 30, 
  minDist = 0.5,
  metric = "cosine",
  dimsToUse = 1:30,
  force = TRUE
)
saveRDS(proj_epi,"ArchR_epi_proj.rds")
PCa_ArchR_epi_proj <- readRDS("D:/Urinary system tumors/work/6_scATAC_PCa_cluster/PCa/PCa_ArchR_epi_proj.rds")
scRNA <- readRDS('D:/Urinary system tumors/work/4_merge/PCa_sce_anno_second.rds')
scRNA= scRNA[,scRNA@meta.data$cell_type %in% c("Epithelium")]
PCa_ArchR_epi_proj <- addGeneIntegrationMatrix(
  ArchRProj   = PCa_ArchR_epi_proj,
  useMatrix   = "GeneScoreMatrix",
  matrixName  = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA       = scRNA,
  addToArrow  = FALSE,
  groupRNA    = "cell_type_second",
  nameCell    = "predictedCell_Un",
  nameGroup   = "predictedGroup_Un",
  nameScore   = "predictedScore_Un"
)
groupList <- SimpleList(
  cluster = SimpleList(
    ATAC = PCa_ArchR_epi_proj$cellNames,
    RNA = colnames(scRNA)
  )
)

PCa_ArchR_epi_proj <- addGeneIntegrationMatrix(
  ArchRProj = PCa_ArchR_epi_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scRNA,
  groupList = groupList,
  groupRNA = "cell_type_second",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",force = TRUE
)
plotEmbedding(ArchRProj = PCa_ArchR_epi_proj, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
x <- as.data.frame(str_split_fixed(PCa_ArchR_epi_proj@cellColData@rownames,"#", n = 2))
x$gleason <- 0
x[which(x$V1 == "SRR14151148"),]$gleason <- "5+5"
x[which(x$V1 == "SRR14151150"),]$gleason <- "4+4"
x[which(x$V1 == "SRR14151151"),]$gleason <- "4+4"
x[which(x$V1 == "SRR14151152"),]$gleason <- "3+3"
x[which(x$V1 == "SRR14151155"),]$gleason <- "3+3"
x[which(x$V1 == "SRR14151156"),]$gleason <- "3+4"
x[which(x$V1 == "SRR14151157"),]$gleason <- "4+4"
x[which(x$V1 == "SRR14151158"),]$gleason <- "3+3"
x[which(x$V1 == "SRR14151159"),]$gleason <- "4+4"
x[which(x$V1 == "SRR14151160"),]$gleason <- "3+3"
x[which(x$V1 == "SRR14151161"),]$gleason <- "3+3"
x[which(x$V1 == "SRR14151162"),]$gleason <- "4+4"
x[which(x$V1 == "SRR14151163"),]$gleason <- "3+4"
x[which(x$V1 == "SRR14151164"),]$gleason <- "3+4"
x[which(x$V1 == "SRR14151165"),]$gleason <- "4+3"

PCa_ArchR_epi_proj@cellColData@listData[["gleason"]] <- x$gleason
plotEmbedding(ArchRProj = PCa_ArchR_epi_proj, colorBy = "cellColData", name = "gleason", embedding = "UMAP30")

x$group <- 0
x[which(x$gleason == "3+3"),]$group <- "G3"
x[which(x$gleason == "3+4"),]$group <- "G3"
x[which(x$gleason == "4+3"),]$group <- "G3"
x[which(x$gleason == "4+4"),]$group <- "G4"
x[which(x$gleason == "5+5"),]$group <- "G5"
PCa_ArchR_epi_proj@cellColData@listData[["group"]] <- x$group

cell <- PCa_ArchR_epi_proj@cellColData@rownames[which(PCa_ArchR_epi_proj@cellColData@listData[["group"]] %in% c(c("G3","G4")))]
proj_CancerTrajector_subset <- subsetArchRProject(
  ArchRProj = PCa_ArchR_epi_proj,
  cells = cell,
  outputDirectory = "CancerTrajector_subset",
  dropCells = TRUE,
  force = TRUE
)
plotEmbedding(ArchRProj = proj_CancerTrajector_subset, colorBy = "cellColData", name = "group", embedding = "UMAP30")

trajectory <- c("G3","G4")
proj_CancerTrajector_subset <- addTrajectory(
  ArchRProj = proj_CancerTrajector_subset,
  name = "CancerTrajector",
  groupBy = "group",
  trajectory = trajectory,
  embedding = "UMAP30",
  force = TRUE
)

###
p <- plotTrajectory(proj_CancerTrajector_subset, trajectory = "CancerTrajector", colorBy = "cellColData", name = "CancerTrajector",embedding = "UMAP30")
p[[1]]


#
proj_CancerTrajector_subset <- addGroupCoverages(ArchRProj = proj_CancerTrajector_subset ,groupBy = "group",force = T) 
proj_CancerTrajector_subset <- addReproduciblePeakSet(ArchRProj= proj_CancerTrajector_subset  ,groupBy = "group", 
                                                      peakMethod = 'Tiles', method = 'p',force = T)
proj_CancerTrajector_subset <- addPeakMatrix(proj_CancerTrajector_subset ,force = T)
proj_CancerTrajector_subset <- addBgdPeaks(proj_CancerTrajector_subset,force = TRUE)
#saveRDS(proj_CancerTrajector_subset,"proj_CancerTrajector_subset.rds")
# motif
#if("Motif" %ni% names(proj_CancerTrajector_subset@peakAnnotation)){
proj_CancerTrajector_subset <- addMotifAnnotations(ArchRProj = proj_CancerTrajector_subset, motifSet = "cisbp", name = "Motif")
#}
proj_CancerTrajector_subset <- addDeviationsMatrix(
  ArchRProj = proj_CancerTrajector_subset, 
  peakAnnotation = "Motif",
  force = TRUE
)

###
trajMM  <- getTrajectory(ArchRProj = proj_CancerTrajector_subset, name = "CancerTrajector", useMatrix = "MotifMatrix", log2Norm = FALSE, trajectoryLabel = "group")
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)))
plotPDF(p1, name = "pseutime_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)

trajGSM <- getTrajectory(ArchRProj = proj_CancerTrajector_subset, name = "CancerTrajector", useMatrix = "GeneScoreMatrix", log2Norm = TRUE, trajectoryLabel = "group")
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGSM)$label)))
plotPDF(p2, name = "pseutime_gene_score.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)


trajGIM <- getTrajectory(ArchRProj = proj_CancerTrajector_subset, name = "CancerTrajector", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE, trajectoryLabel = "group")
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGIM)$label)))
plotPDF(p3, name = "pseutime_gene_expr.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)

trajPM  <- getTrajectory(ArchRProj = proj_CancerTrajector_subset, name = "CancerTrajector", useMatrix = "PeakMatrix", log2Norm = TRUE, trajectoryLabel = "group")
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)))
plotPDF(p4, name = "pseutime_peak.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)


##GENE SCORE AND MOTIF
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "score_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)


##GENE INTER AND MOTIF
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"),force = T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "expr_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)

##Peak2Peak
proj_CancerTrajector_subset <- addCoAccessibility(
  ArchRProj = proj_CancerTrajector_subset,
  dimsToUse = 1:30,
  reducedDims = "IterativeLSI30"
)
cA <- getCoAccessibility(
  ArchRProj = proj_CancerTrajector_subset,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA_loop <- getCoAccessibility(
  ArchRProj = proj_CancerTrajector_subset,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = T
)
cA_data <- as.data.frame(cA)
cA_loop_data <- as.data.frame(cA_loop)

proj_CancerTrajector_subset <- addPeak2GeneLinks(
  ArchRProj = proj_CancerTrajector_subset,
  reducedDims = "IterativeLSI30"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = proj_CancerTrajector_subset,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2geneDF <-  p2g
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]

annot <- readRDS(metadata(p2geneDF)$seATAC)
p2geneDF$peakType <- annot@rowRanges$peakType[p2geneDF$idxATAC]
p2g.df.obs <- as.data.frame(p2geneDF)
p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]
p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$Correlation > quantile(p2g.df.obs$Correlation,0.5)),]
TF_motif_pair1 <- data.frame(TFs =trajGIM2@elementMetadata@listData[["name"]],
                             motifs = trajMM2@elementMetadata@listData[["name"]])

TF_motif_pair2 <- data.frame(TFs =trajGSM2_2@elementMetadata@listData[["name"]],
                             motifs = trajMM2_2@elementMetadata@listData[["name"]])
intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)
TF_motif_pair <- TF_motif_pair1[which(TF_motif_pair1$TFs %in% intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)),]

#peaks
Trajectory_genes <- trajGIM@elementMetadata@listData[["name"]]
#genes
Trajectory_peaks <- as.data.frame(trajPM@elementMetadata@listData)
Trajectory_peaks <- Trajectory_peaks[,-2]
peaks <- paste(paste(Trajectory_peaks$seqnames,Trajectory_peaks$start,sep = ":"),Trajectory_peaks$end,sep = "-")

p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$geneName %in% Trajectory_genes),]
p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$peakName %in% peaks ),]

##GENE SCORE AND MOTIF
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
TF_motif_pair2 <- data.frame(TFs =trajGSM2@elementMetadata@listData[["name"]],
                             motifs = trajMM2@elementMetadata@listData[["name"]])


#GENE INTER AND MOTIF
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))

TF_motif_pair1 <- data.frame(TFs =trajGIM2@elementMetadata@listData[["name"]],
                             motifs = trajMM2@elementMetadata@listData[["name"]])


intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)
TF_motif_pair <- TF_motif_pair1[which(TF_motif_pair1$TFs %in% intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)),]

#peaks
Trajectory_genes <- trajGIM@elementMetadata@listData[["name"]]
#genes
Trajectory_peaks <- as.data.frame(trajPM@elementMetadata@listData)
Trajectory_peaks <- Trajectory_peaks[,-2]
peaks <- paste(paste(Trajectory_peaks$seqnames,Trajectory_peaks$start,sep = ":"),Trajectory_peaks$end,sep = "-")

p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$geneName %in% Trajectory_genes),]
p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$peakName %in% peaks ),]


peakset=as.data.frame(getPeakSet(proj_CancerTrajector_subset),row.names=NULL)
anno_match <- getMatches(ArchRProj = proj_CancerTrajector_subset, name = "Motif")
#anno_match <- readRDS("D:/Urinary system tumors/work/6_scATAC/ArchRSubset/Annotations/Motif-Matches-In-Peaks.rds")
all_motif_match <- anno_match@assays@data@listData[["matches"]]
rownames(all_motif_match) <- paste0(peakset$seqnames,":",peakset$start,"-",peakset$end)
all_motif_match_data_frame <- reshape2::melt(as.matrix(all_motif_match))

motif_match_data_frame <- all_motif_match_data_frame[all_motif_match_data_frame$Var2 %in% TF_motif_pair$motifs,]
motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$value=="TRUE"),]

motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$Var1 %in% p2g.df.obs$peakName),]
motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$Var2 %in% TF_motif_pair$motifs),]
motif_match_data_frame <- motif_match_data_frame[,-3]
colnames(motif_match_data_frame) <- c("peakName","motif")

##motif_peak_gene
motif_peak_gene <- left_join(p2g.df.obs, motif_match_data_frame, by="peakName")
write.csv(motif_peak_gene,"PCa_motif_peak_gene.csv",row.names = F)

color <- c("G3" = "#D8EBD3",
           "G4" = "#228A87")
plotEmbedding(ArchRProj = proj_CancerTrajector_subset, colorBy = "cellColData", name = "group",pal = color ,embedding = "UMAP30")
p <- plotTrajectory(proj_CancerTrajector_subset, trajectory = "CancerTrajector", colorBy = "cellColData", name = "group",embedding = "UMAP30")

intersect(corGIM_MM[[1]]$name1,corGSM_MM[[1]]$name1)
intersect(corGIM_MM[[1]]$name2,corGSM_MM[[1]]$name2)

ht3 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "gene_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)




##GENE SCORE AND MOTIF
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "score_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 6, height = 8)


##GENE INTER AND MOTIF
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"),force = T,
                             varCutOff = 0, rowOrder = rowOrder)

intersect(corGIM_MM[[1]]$name1,corGSM_MM[[1]]$name1)
intersect(corGIM_MM[[1]]$name2,corGSM_MM[[1]]$name2)

name1 <- c("chr1:ELF3","chr11:TEAD1","chr13:KLF5","chr14:FOS",
           "chr19:FOSB","chr20:CEBPB","chr21:ETS2","chr8:KLF10")
name2 <- c("z:ELF3_342","z:TEAD1_796","z:KLF5_175","z:FOS_137","z:FOSB_121",       
           "z:CEBPB_140","z:ETS2_340","z:KLF10_826" )
trajGSM2 <- trajGSM[name1, ]
trajMM2 <- trajMM[name2, ]
trajGIM2 <- trajGIM[name1, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht3 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2+ht3, name = "expr_score_motif.pdf", ArchRProj = proj_CancerTrajector_subset, addDOC = FALSE, width = 16, height = 8)

pwm <- getPeakAnnotation(proj_CancerTrajector_subset, "Motif")$motifs[["KLF5_175"]]
pwm


PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
ppm


library(ggseqlogo)
ggseqlogo(ppm, method = "bits")

##ccRCC scATAC------
proj  <- loadArchRProject(outputDirectory)
cell <- proj@cellColData@rownames[which(grepl("Epithelium",proj@cellColData@listData[["predictedGroup_Un"]]))]

proj_epi <- subsetArchRProject(
  ArchRProj       = proj,
  cells           = cell,
  outputDirectory = outputDirectory_epi,
  dropCells       = FALSE,
  force           = TRUE)

### cell_type_second
Idents(scRNA) <- scRNA$cell_type
table(Idents(scRNA))
scRNA_epi <- subset(scRNA, idents = 'Epithelium')
Idents(scRNA_epi) <- scRNA_epi$cell_type_second
table(Idents(scRNA_epi))



proj_epi <- addGeneIntegrationMatrix(
  ArchRProj   = proj_epi, 
  useMatrix   = "GeneScoreMatrix",
  matrixName  = "GeneIntegrationMatrix_second",
  reducedDims = "Harmony",
  seRNA      = scRNA_epi,
  addToArrow = TRUE,
  groupRNA   = "cell_type_second",
  nameCell   = "predictedCell_second_Un",
  nameGroup  = "predictedGroup_second_Un",
  nameScore  = "predictedScore_second_Un",
  force = TRUE)



proj_epi@cellColData$WHO <- proj_epi@cellColData$Sample
head(proj_epi@cellColData)
proj_epi@cellColData$WHO <- gsub( 'RCC81','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC84','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC86','I' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC87','III',proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC94','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC96','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub( 'RCC99','I' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC100','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC101','I' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC103','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC104','III',proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC106','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC112','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC113','I' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC114','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC115','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC116','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC119','II' ,proj_epi@cellColData$WHO)
proj_epi@cellColData$WHO <- gsub('RCC120','II' ,proj_epi@cellColData$WHO)

proj_epi@cellColData$WHO[which(grepl("Normal epithelium",proj_epi@cellColData@listData[["predictedGroup_second_Un"]]))] <- "Normal"
table(proj_epi@cellColData$WHO)


### Peak calling
proj_epi <- addGroupCoverages(ArchRProj        = proj_epi, groupBy = "predictedGroup_second_Un")
proj_epi <- addReproduciblePeakSet(ArchRProj   = proj_epi, groupBy = "predictedGroup_second_Un", peakMethod = 'Tiles',method = 'p')
proj_epi <- addPeakMatrix(proj_epi)

markersPeaks <- getMarkerFeatures(
  ArchRProj  = proj_epi, 
  useMatrix  = "PeakMatrix", 
  groupBy    = "predictedGroup_Un",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")
saveRDS(markersPeaks,"ccRCC_markersPeaks_epi.rds")
markerList_GR <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

proj_epi <- addMotifAnnotations(ArchRProj = proj_epi, motifSet  = "cisbp", name = "Motif",force = TRUE)

### peakAnnoEnrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker       = markersPeaks,
  ArchRProj      = proj_epi,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# Peak2GeneLinkages分析
getAvailableMatrices(proj_epi)
proj_epi <- addCoAccessibility(ArchRProj = proj_epi,reducedDims = "Harmony")
proj_epi <- addPeak2GeneLinks(ArchRProj  = proj_epi,reducedDims = "Harmony")

saveArchRProject(ArchRProj = proj_epi, outputDirectory = outputDirectory_epi, load = TRUE)

proj_epi <- addIterativeLSI(ArchRProj        = proj_epi, 
                            useMatrix        = "TileMatrix", 
                            name             = "IterativeLSI", 
                            outlierQuantiles = NULL,
                            force            = TRUE)

proj_epi <- addHarmony(ArchRProj   = proj_epi,
                       reducedDims = "IterativeLSI",
                       name        = "Harmony",
                       groupBy     = "Sample",
                       force       = TRUE)

proj_epi <- addClusters(input       = proj_epi,
                        reducedDims = "Harmony",
                        method      = "Seurat",
                        name        = "Clusters",
                        resolution  = 0.6,
                        maxClusters = 12,
                        force       = TRUE)
proj_epi <- addImputeWeights(proj_epi)
proj_epi <- addUMAP(proj_epi,reducedDims = "Harmony",name = "UMAPHarmony",nNeighbors = 30,minDist = 0.5,metric = "cosine",force = TRUE)

p <- plotEmbedding(proj_epi,   colorBy = "cellColData",   name = "predictedGroup_second_Un", embedding ='UMAPHarmony')
plotPDF(p, name = "UMAP-Harmony-RNA-Integration.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(proj_epi,   colorBy = "cellColData",   name = "WHO", embedding ='UMAPHarmony')
plotPDF(p, name = "UMAP-Harmony-RNA-Integration_WHO.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)



proj_epi <- addPeakMatrix(proj_epi)

# getMarkerFeatures

table(proj_epi@cellColData$WHO)
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_epi, 
  useMatrix  = "PeakMatrix",
  groupBy    = "WHO",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  useGroups  = "Tumor cells",
  bgdGroups  = "Normal epithelium")
saveRDS(markerTest,"F:/Mingjie/Project_Urinary/01.Data/04.OutputData_Epi/scATAC_Data_markersPeaks_NormalTumor.rds")

saveArchRProject(ArchRProj = proj_epi, outputDirectory = outputDirectory_epi, load = TRUE)

proj_epi <- loadArchRProject(outputDirectory_epi)
getAvailableMatrices(proj_epi)

# trajectory <- c('Normal','I','II','III')
trajectory <- c('II','III')

proj_epi <- addTrajectory(
  ArchRProj   = proj_epi,
  name        = "CancerTrajector",
  groupBy     = "WHO",
  trajectory  = trajectory,
  reducedDims = "IterativeLSI",
  embedding   = "UMAPHarmony",
  force       = TRUE)

proj_epi@cellColData$WHO <- as.character(proj_epi@cellColData$WHO)
table(proj_epi@cellColData$WHO)

color <- c('I'   = "#F8CED7",
           "II"  = "#D66E71",
           "III" = "#AA1F24")
p1 <- plotEmbedding(ArchRProj = proj_epi, colorBy = "cellColData", name = "WHO",pal = color ,embedding = "UMAPHarmony")
p2 <- plotTrajectory(proj_epi, trajectory = "CancerTrajector", colorBy = "cellColData", name = "CancerTrajector", embedding = "UMAPHarmony")
p3 <- plotTrajectory(proj_epi, trajectory = 'CancerTrajector', colorBy = "ColData",     name = 'CancerTrajector', continuousSet = "greenBlue",embedding = "UMAPHarmony")
plotPDF(p1,p2,p3, name = "Trajectory-CancerTrajector.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 5, height = 5)

# proj_epi <- addArchRAnnotations(ArchRProj = proj_epi, db = "ArchR", collection = "ATAC" ,force = TRUE)
proj_epi <- addMotifAnnotations(ArchRProj = proj_epi, motifSet = "cisbp", name = "Motif",force = TRUE)
proj_epi <- addBgdPeaks(        ArchRProj = proj_epi, method = "ArchR", force = TRUE)

proj_epi <- addDeviationsMatrix(
  ArchRProj      = proj_epi,
  peakAnnotation = "Motif",
  bgdPeaks       = getBgdPeaks(proj_epi, method = "ArchR"),
  matrixName     = NULL,
  out            = c("z", "deviations"),
  binarize       = FALSE,
  version        = 2,
  threads        = getArchRThreads(),
  verbose        = TRUE,
  parallelParam  = NULL,
  force          = FALSE,
  logFile        = createLogFile("addDeviationsMatrix"))

proj_epi <- addPeak2GeneLinks(  ArchRProj = proj_epi, reducedDims = "IterativeLSI")
proj_epi <- addImputeWeights(   ArchRProj = proj_epi)
proj_epi <- addGroupCoverages(  ArchRProj = proj_epi, groupBy = "WHO",force = TRUE)


trajMM  <- getTrajectory(ArchRProj = proj_epi, name = "CancerTrajector", useMatrix = "MotifMatrix", log2Norm = FALSE, trajectoryLabel = "WHO")
p1      <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajMM)$label)))
plotPDF(p1, name = "pseutime_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)

trajGSM <- getTrajectory(ArchRProj = proj_epi, name = "CancerTrajector", useMatrix = "GeneScoreMatrix", log2Norm = TRUE, trajectoryLabel = "WHO")
p2 <- plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGSM)$label)))
plotPDF(p2, name = "pseutime_score.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)

trajGIM <- getTrajectory(ArchRProj = proj_epi, name = "CancerTrajector", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE, trajectoryLabel = "WHO")
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajGIM)$label)))
plotPDF(p3, name = "pseutime_expr.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)

trajPM  <- getTrajectory(ArchRProj = proj_epi, name = "CancerTrajector", useMatrix = "PeakMatrix", log2Norm = TRUE, trajectoryLabel = "WHO")
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"), colorColumns = TRUE, columnPal = paletteDiscrete(values = unique(colData(trajPM)$label)))
plotPDF(p4, name = "pseutime_peak.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)


##GEME SCORE AND MOTIF
corGSM_MM    <- correlateTrajectories(trajGSM, trajMM)
trajGSM2     <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2      <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder    <- match(rownames(combinedMat), rownames(trajGSM2))
ht1         <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2         <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "score_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)


##GEME INTER AND MOTIF
corGIM_MM   <- correlateTrajectories(trajGIM, trajMM)
trajGIM2    <- trajGIM[corGIM_MM[[1]]$name1, ]
trajGSM2_2  <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2     <- trajMM[corGIM_MM[[1]]$name2, ]
trajMM2_2   <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1      <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2      <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"),force = T, varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "expr_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)

##Peak2Peak

proj_epi <- addCoAccessibility(
  ArchRProj   = proj_epi,
  dimsToUse   = 1:30,
  reducedDims = "IterativeLSI")

cA <- getCoAccessibility(
  ArchRProj   = proj_epi,
  corCutOff   = 0.5,
  resolution  = 1,
  returnLoops = FALSE)

cA_loop <- getCoAccessibility(
  ArchRProj   = proj_epi,
  corCutOff   = 0.5,
  resolution  = 1,
  returnLoops = T)

cA_data      <- as.data.frame(cA)
cA_loop_data <- as.data.frame(cA_loop)

### Peak2Gene 
proj_epi <- addPeak2GeneLinks(
  ArchRProj = proj_epi,
  reducedDims = "IterativeLSI")

p2g <- getPeak2GeneLinks(
  ArchRProj = proj_epi,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE)

p2geneDF <-  p2g
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]

annot <- readRDS(metadata(p2geneDF)$seATAC)
p2geneDF$peakType <- annot@rowRanges$peakType[p2geneDF$idxATAC]
p2g.df.obs <- as.data.frame(p2geneDF)
p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]
p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$Correlation > quantile(p2g.df.obs$Correlation,0.5)),]

TF_motif_pair1 <- data.frame(TFs    = trajGIM2@elementMetadata@listData[["name"]],
                             motifs = trajMM2@elementMetadata@listData[["name"]])

TF_motif_pair2 <- data.frame(TFs    = trajGSM2_2@elementMetadata@listData[["name"]],
                             motifs = trajMM2_2@elementMetadata@listData[["name"]])

intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)

TF_motif_pair <- TF_motif_pair1[which(TF_motif_pair1$TFs %in% intersect(TF_motif_pair1$TFs,TF_motif_pair2$TFs)),]

#peaks
Trajectory_genes <- trajGIM@elementMetadata@listData[["name"]]
#genes
Trajectory_peaks <- as.data.frame(trajPM@elementMetadata@listData)
Trajectory_peaks <- Trajectory_peaks[,-2]
peaks <- paste(paste(Trajectory_peaks$seqnames,Trajectory_peaks$start,sep = ":"),Trajectory_peaks$end,sep = "-")

p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$geneName %in% Trajectory_genes),]
p2g.df.obs <- p2g.df.obs[which(p2g.df.obs$peakName %in% peaks ),]


#
peakset=as.data.frame(getPeakSet(proj_epi),row.names=NULL)
anno_match <- getMatches(ArchRProj = proj_epi, name = "Motif")
#anno_match <- readRDS("D:/Urinary system tumors/work/6_scATAC/ArchRSubset/Annotations/Motif-Matches-In-Peaks.rds")
all_motif_match <- anno_match@assays@data@listData[["matches"]]
rownames(all_motif_match) <- paste0(peakset$seqnames,":",peakset$start,"-",peakset$end)
all_motif_match_data_frame <- reshape2::melt(as.matrix(all_motif_match))

motif_match_data_frame <- all_motif_match_data_frame[all_motif_match_data_frame$Var2 %in% TF_motif_pair$motifs,]
motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$value=="TRUE"),]

motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$Var1 %in% p2g.df.obs$peakName),]
motif_match_data_frame <- motif_match_data_frame[which(motif_match_data_frame$Var2 %in% TF_motif_pair$motifs),]
motif_match_data_frame <- motif_match_data_frame[,-3]
colnames(motif_match_data_frame) <- c("peakName","motif")

##motif_peak_gene
library(dplyr)
motif_peak_gene <- left_join(p2g.df.obs, motif_match_data_frame, by="peakName")
write.csv(motif_peak_gene,"PCa_motif_peak_gene.csv",row.names = F)

intersect(corGIM_MM[[1]]$name1,corGSM_MM[[1]]$name1)
intersect(corGIM_MM[[1]]$name2,corGSM_MM[[1]]$name2)

ht3 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "score_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)


##GEME SCORE AND MOTIF
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
trajGSM2  <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2   <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "score_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 6, height = 8)


##GEME INTER AND MOTIF
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"),force = T,
                             varCutOff = 0, rowOrder = rowOrder)

intersect(corGIM_MM[[1]]$name1,corGSM_MM[[1]]$name1)
intersect(corGIM_MM[[1]]$name2,corGSM_MM[[1]]$name2)

name1 <- c( "chr11:FOSL1", "chr19:NFIC",  "chr20:HNF4A", "chr20:CEBPB", "chr21:RUNX1", "chr6:ID4", "chr9:NFIB")
name2 <- c( "z:FOSL1_142", "z:NFIC_740",  "z:HNF4A_662", "z:CEBPB_140", "z:RUNX1_733", "z:ID4_75", "z:NFIB_741")
trajGSM2 <- trajGSM[name1, ]
trajMM2 <- trajMM[name2, ]
trajGIM2 <- trajGIM[name1, ]
trajCombined <- trajGSM2
assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht3 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"),force=T,
                             varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2+ht3, name = "expr_score_motif.pdf", ArchRProj = proj_epi, addDOC = FALSE, width = 10, height = 8)

pwm <- getPeakAnnotation(proj_epi, "Motif")$motifs[["ID4_75"]]
#pwm


PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
ppm


library(ggseqlogo)
ggseqlogo(ppm, method = "bits")
ggsave("F:\\Mingjie\\Project_Urinary\\01.Data\\04.OutputData_Epi\\Plots\\Trajectory_gene\\ID4_75.pdf")

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
saveArchRProject(ArchRProj = proj_epi, outputDirectory = outputDirectory_epi, load = TRUE)

#TF-gene network ------
#PCa
CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/filter/PCa_cnv_gene.txt", header = T,sep = "\t")
CNV_gene <- na.omit(CNV_gene)
Time_gene <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/PCa/2_PCa_tumor_modulated_genes.txt", header = T,sep = "\t")

Time_CNV_gene <- Time_gene[unique(CNV_gene$gene),]

TF_gene <- read.csv("D:/PCa_ATAC/PCa_motif_peak_gene.csv", header = T)
#TF_gene <-distinct(TF_gene[,c(7,10)])
TF_gene <- TF_gene[which(TF_gene$geneName %in% Time_CNV_gene$gene_short_name),]
TF_gene <- na.omit(TF_gene)
write.table(unique(c("FOSB",
                     "FOS",
                     "KLF5",
                     "TEAD1",
                     "ETS2",
                     "ELF3",
                     "KLF10",
                     "CEBPB",
                     TF_gene$geneName)),"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\PPI.txt",quote = F,sep = "\t",row.names = F)
write.table(TF_gene,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\TF_gene_sup.txt",quote = F,sep = "\t",row.names = F)

TF <- c("FOSB",
        "FOS",
        "KLF5",
        "TEAD1",
        "ETS2",
        "ELF3",
        "KLF10",
        "CEBPB")
colnames(TF_gene) <- c("Target","Source")
TF_gene[which(TF_gene$Source == "FOSB_121"),]$Source <- "FOSB"
TF_gene[which(TF_gene$Source == "FOS_137"),]$Source <- "FOS"
TF_gene[which(TF_gene$Source == "KLF5_175"),]$Source <- "KLF5"
TF_gene[which(TF_gene$Source == "TEAD1_796"),]$Source <- "TEAD1"
TF_gene[which(TF_gene$Source == "ETS2_340"),]$Source <- "ETS2"
TF_gene[which(TF_gene$Source == "ELF3_342"),]$Source <- "ELF3"
TF_gene[which(TF_gene$Source == "KLF10_826"),]$Source <- "KLF10"
TF_gene[which(TF_gene$Source == "CEBPB_140"),]$Source <- "CEBPB"

TF_gene$Direction <- "Direction"
#Through string URLhttps://string-db.org/
#Enter the genes in the ppi. txt file to obtain the ppi network
PPI <- read.delim("clipboard",header = T)#
PPI <- PPI[,c(1,2)]
colnames(PPI) <- c("Target","Source")
PPI$Direction<- "none"
x <- rbind(TF_gene,PPI)
y <- rbind(cbind(TF_gene$Target,"gene"),cbind(TF_gene$Source,"TF"))

library(igraph)
c(as.character(x$Target),as.character(x$Source)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices


colnames(vertices) <- c("node", "n")
vertices <- vertices[order(vertices$n,decreasing = T),]
vertices <- vertices[1:108,]
x <- x[which(x$Target %in% vertices$node & x$Source %in% vertices$node),]
y <- as.data.frame(y[which(y[,1] %in% vertices$node),])
y <- distinct(y)

x$merge <- paste(x$Target,x$Source)
sort(table(x$merge),decreasing = T)
write.table(x,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\Network_Figre3D.txt.txt",quote = F,sep = "\t",row.names = F)
write.table(y,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\PCa\\Network_Figre3D.txt_table.txt",quote = F,sep = "\t",row.names = F)
###Building TF-gene and PPI network by Cytoscape
#ccRCC
CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/filter/ccRCC_cnv_gene.txt", header = T,sep = "\t")
CNV_gene <- na.omit(CNV_gene)
Time_gene <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/ccRCC/2_ccRCC_tumor_modulated_genes.txt", header = T,sep = "\t")

Time_CNV_gene <- Time_gene[unique(CNV_gene$gene),]
TF_gene <- read.csv("D:\\ccRCC_scATAC/ccRCC_motif_peak_gene.csv", header = T)
TF_gene <- TF_gene[which(TF_gene$geneName %in% Time_CNV_gene$gene_short_name),]
TF_gene <- na.omit(TF_gene)
TF_gene <- TF_gene[which(TF_gene$motif %in% c("FOSL1_142",
                                              "RUNX1_733",
                                              "NFIB_741",
                                              "NFIC_740",
                                              "CEBPB_140",
                                              "HNF4A_662",
                                              "ID4_75")),]

TF_gene <-distinct(TF_gene[,c(7,10)])
TF_gene <- TF_gene[which(TF_gene$geneName %in% Time_CNV_gene$gene_short_name),]
TF_gene <- na.omit(TF_gene)


colnames(TF_gene) <- c("Target","Source")
TF_gene[which(TF_gene$Source == "FOSL1_142"),]$Source <- "FOSL1"
TF_gene[which(TF_gene$Source == "RUNX1_733"),]$Source <- "RUNX1"
TF_gene[which(TF_gene$Source == "NFIB_741"),]$Source <- "NFIB"
TF_gene[which(TF_gene$Source == "NFIC_740"),]$Source <- "NFIC"
TF_gene[which(TF_gene$Source == "CEBPB_140"),]$Source <- "CEBPB"
TF_gene[which(TF_gene$Source == "HNF4A_662"),]$Source <- "HNF4A"
TF_gene[which(TF_gene$Source == "ID4_75"),]$Source <- "ID4"


write.table(unique(c(TF_gene$Target,TF_gene$Source)),"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\PPI.txt",quote = F,sep = "\t",row.names = F)
write.table(TF_gene,"D:\\Urinary system tumors\\work\\3_Pseutime\\1_tumor_epi\\ccRCC\\TF_gene.txt",quote = F,sep = "\t",row.names = F)


colnames(TF_gene) <- c("Target","Source")
TF_gene$Direction <- "Direction"
#Through string URLhttps://string-db.org/
#Enter the genes in the ppi. txt file to obtain the ppi network
PPI <- read.delim("clipboard",header = T)
PPI <- PPI[,c(1,2)]
colnames(PPI) <- c("Target","Source")
PPI$Direction<- "none"
x <- rbind(TF_gene,PPI)
y <- rbind(cbind(TF_gene$Target,"gene"),cbind(TF_gene$Source,"TF"))
library(igraph)
c(as.character(x$Target),as.character(x$Source)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")
vertices <- vertices[order(vertices$n,decreasing = T),]
vertices <- vertices[1:107,]
x <- x[which(x$Target %in% vertices$node & x$Source %in% vertices$node),]
y <- as.data.frame(y[which(y[,1] %in% vertices$node),])
y <- distinct(y)
write.table(x,"Network_Figre3D.txt",quote = F,sep = "\t",row.names = F)
write.table(y,"Network_Figre3D_table.txt",quote = F,sep = "\t",row.names = F)