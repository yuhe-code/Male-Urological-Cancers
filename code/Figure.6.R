#Run CYTOTRACE ----
library(Seurat)
T_cell <- readRDS("D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno.rds")
Myeloid <- readRDS("D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_sce_anno.rds")
Mast <- readRDS("D:\\Urinary system tumors\\work\\4_Mast_cells\\2_Mast_cells_anno\\Mast_cells_sce_anno.rds")
Fibroblast <- readRDS("D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\Fibroblast_sce_anno.rds")
Endothelium <- readRDS("D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_sce_anno.rds")
B_cells <- readRDS("D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_sce_anno.rds")
##CYTOTRACE
library(CytoTRACE2)
T_cell <- cytotrace2(T_cell, 
                    is_seurat = T, 
                    slot_type = "counts", 
                    species = 'human',
                    seed = 1234)

Myeloid <- cytotrace2(Myeloid, 
                     is_seurat = T, 
                     slot_type = "counts", 
                     species = 'human',
                     seed = 1234)

Mast <- cytotrace2(Mast, 
                     is_seurat = T, 
                     slot_type = "counts", 
                     species = 'human',
                     seed = 1234)

Fibroblast <- cytotrace2(Fibroblast, 
                     is_seurat = T, 
                     slot_type = "counts", 
                     species = 'human',
                     seed = 1234)

Endothelium <- cytotrace2(Endothelium, 
                     is_seurat = T, 
                     slot_type = "counts", 
                     species = 'human',
                     seed = 1234)

B_cells <- cytotrace2(B_cells, 
                     is_seurat = T, 
                     slot_type = "counts", 
                     species = 'human',
                     seed = 1234)
save.image("D:/Urinary system tumors/work/3_Pseutime/2_Other_cell_monocle3.RData")
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Other_cell_monocle3.RData")
#T cell monocle3 -----
#rm(list=setdiff(ls(), "T_cell"))
ggplot(T_cell@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
T_cell@meta.data$monocle3_cluster <- T_cell@meta.data$cluster
##ccRCC
ccRCC_CTL <-  subset(T_cell, tumor_type %in% c("ccRCC"))
ccRCC_CTL <-  subset(ccRCC_CTL, T_cell_type %in% c("CTL"))

ccRCC_CTL <- ccRCC_CTL %>% RunTSNE(dims = 1:11) %>% RunUMAP(dims = 1:11)

data <- GetAssayData(ccRCC_CTL, assay = 'RNA', slot = 'counts')
cell_metadata <- ccRCC_CTL@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
ccRCC_CTL_cds <- new_cell_data_set(data,
                                   cell_metadata = cell_metadata,
                                   gene_metadata = gene_annotation)
ccRCC_CTL_cds <- preprocess_cds(ccRCC_CTL_cds, num_dim = 11)#
plot_pc_variance_explained(ccRCC_CTL_cds)
ccRCC_CTL_cds <- reduce_dimension(ccRCC_CTL_cds,reduction_method = 'UMAP')
cds.embed <- ccRCC_CTL_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ccRCC_CTL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
ccRCC_CTL_cds@int_colData$reducedDims$UMAP <- int.embed
ccRCC_CTL_cds <- cluster_cells(ccRCC_CTL_cds,reduction_method = "UMAP")

plot_cells(ccRCC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")

#ccRCC_CTL_cds <- learn_graph(ccRCC_CTL_cds,learn_graph_control=list(ncenter=500),use_partition = TRUE)#Control the number of branch nodes
ccRCC_CTL_cds <- learn_graph(ccRCC_CTL_cds,use_partition = TRUE)
plot_cells(ccRCC_CTL_cds, color_cells_by = "CytoTRACE2_Relative")
ccRCC_CTL_cds <- order_cells(ccRCC_CTL_cds,reduction_method = "UMAP")
plot_cells(ccRCC_CTL_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
##
##PCa
PCa_CTL <-  subset(T_cell, tumor_type %in% c("PCa"))
PCa_CTL <-  subset(PCa_CTL, T_cell_type %in% c("CTL"))

PCa_CTL <- PCa_CTL %>% RunTSNE(dims = 1:11) %>% RunUMAP(dims = 1:11)
FeaturePlot(PCa_CTL,features = "CytoTRACE2_Relative")
DimPlot(PCa_CTL,group.by = "cluster")
data <- GetAssayData(PCa_CTL, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa_CTL@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
PCa_CTL_cds <- new_cell_data_set(data,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_annotation)
PCa_CTL_cds <- preprocess_cds(PCa_CTL_cds, num_dim = 11)##
plot_pc_variance_explained(PCa_CTL_cds)
PCa_CTL_cds <- reduce_dimension(PCa_CTL_cds,reduction_method = 'UMAP')
cds.embed <- PCa_CTL_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(PCa_CTL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
PCa_CTL_cds@int_colData$reducedDims$UMAP <- int.embed
PCa_CTL_cds <- cluster_cells(PCa_CTL_cds,reduction_method = "UMAP")

plot_cells(PCa_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(PCa_CTL_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


PCa_CTL_cds <- learn_graph(PCa_CTL_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#PCa_CTL_cds <- learn_graph(PCa_CTL_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(PCa_CTL_cds, color_cells_by = "CytoTRACE2_Relative")
PCa_CTL_cds <- order_cells(PCa_CTL_cds,reduction_method = "UMAP")
plot_cells(PCa_CTL_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

##BC
BC_CTL <-  subset(T_cell, tumor_type %in% c("BC"))
BC_CTL <-  subset(BC_CTL, T_cell_type %in% c("CTL"))

BC_CTL <- BC_CTL %>% RunTSNE(dims = 1:11) %>% RunUMAP(dims = 1:11)
FeaturePlot(BC_CTL,features = "CytoTRACE2_Relative")
DimPlot(BC_CTL,group.by = "cluster")
data <- GetAssayData(BC_CTL, assay = 'RNA', slot = 'counts')
cell_metadata <- BC_CTL@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
BC_CTL_cds <- new_cell_data_set(data,
                                cell_metadata = cell_metadata,
                                gene_metadata = gene_annotation)
BC_CTL_cds <- preprocess_cds(BC_CTL_cds, num_dim = 11)##
plot_pc_variance_explained(BC_CTL_cds)
BC_CTL_cds <- reduce_dimension(BC_CTL_cds,reduction_method = 'UMAP')
cds.embed <- BC_CTL_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(BC_CTL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
BC_CTL_cds@int_colData$reducedDims$UMAP <- int.embed
BC_CTL_cds <- cluster_cells(BC_CTL_cds,reduction_method = "UMAP")

plot_cells(BC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(BC_CTL_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


BC_CTL_cds <- learn_graph(BC_CTL_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#BC_CTL_cds <- learn_graph(BC_CTL_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(BC_CTL_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(BC_CTL_cds, color_cells_by = "monocle3_cluster")
BC_CTL_cds <- order_cells(BC_CTL_cds,reduction_method = "UMAP")
plot_cells(BC_CTL_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

rm(list=setdiff(ls(), c("ccRCC_CTL_cds","ccRCC_CTL",
                        "BC_CTL_cds","BC_CTL",
                        "PCa_CTL_cds","PCa_CTL") ))
##chrgene_GSEA----
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol"))
genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)
genesets_hallmark2 <- list()
for(i in unique(genesets_hallmark$gs_name)){
  genesets_hallmark2[[i]] <- genesets_hallmark[which(genesets_hallmark$gs_name == i),]$gene_symbol
}

##ccrcc T
ccRCC_CTL <-  AddModuleScore(ccRCC_CTL,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(ccRCC_CTL_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_CTL@meta.data[row.names(x),c(37:dim(ccRCC_CTL@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\ccRCC_T_cell_hallmarker_pseudotime\\ccRCC_T_cell",colnames(x)[i],".pdf"),width = 4,height = 3)
}

##PCa T 
PCa_CTL <-  AddModuleScore(PCa_CTL,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(PCa_CTL_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_CTL@meta.data[row.names(x),c(37:dim(PCa_CTL@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\PCa_T_cell_hallmarker_pseudotime\\PCa_T_cell",colnames(x)[i],".pdf"),width = 4,height = 3)
}

##BC T 
BC_CTL <-  AddModuleScore(BC_CTL,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(BC_CTL_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_CTL@meta.data[row.names(x),c(37:dim(BC_CTL@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\BC_T_cell_hallmarker_pseudotime\\BC_T_cell",colnames(x)[i],".pdf"),width = 4,height = 3)
}

rm(list=setdiff(ls(), c("ccRCC_CTL_cds","ccRCC_CTL",
                        "BC_CTL_cds","BC_CTL",
                        "PCa_CTL_cds","PCa_CTL") ))

####ccRCC PLOT
plot_cells(ccRCC_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#69DA39","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#69DA39","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_group_legend.pdf"
       ,width = 8, height = 8)
dev.off()
####BC PLOT
plot_cells(BC_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_group_legend.pdf"
       ,width = 8, height = 8)
dev.off()

####PCa PLOT
plot_cells(PCa_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_CTL_cds, color_cells_by = "monocle3_cluster",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("#F8E71C","#D4E423","#B0E02B","#8DDD32","#45D641","#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_group_legend.pdf"
       ,width = 8, height = 8)
dev.off()


##
ccRCC_CTL_modulated_genes <- graph_test(ccRCC_CTL_cds, neighbor_graph = "principal_graph", cores = 4)
BC_CTL_modulated_genes <- graph_test(BC_CTL_cds, neighbor_graph = "principal_graph", cores = 4)
PCa_CTL_modulated_genes <- graph_test(PCa_CTL_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(ccRCC_CTL_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_CTL_modulated_genes.txt",quote = F,sep = "\t")
write.table(BC_CTL_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_CTL_modulated_genes.txt",quote = F,sep = "\t")
write.table(PCa_CTL_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_CTL_modulated_genes.txt",quote = F,sep = "\t")
##ridge plot----
shanji <- as.data.frame(pseudotime(PCa_CTL_cds)[colnames(PCa_CTL_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_CTL_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CTL c1"="#F8E71C",
                             "CTL c2"="#D4E423",
                             "CTL c3"="#B0E02B",
                             "CTL c4"="#8DDD32",
                             "CTL c5"="#69DA39",
                             "CTL c6"="#45D641",
                             "CTL c7"="#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_CTL_shanji.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge-group
shanji <- as.data.frame(pseudotime(ccRCC_CTL_cds)[colnames(ccRCC_CTL_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_CTL_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CTL c1"="#F8E71C",
                             "CTL c2"="#D4E423",
                             "CTL c3"="#B0E02B",
                             "CTL c4"="#8DDD32",
                             "CTL c5"="#69DA39",
                             "CTL c6"="#45D641",
                             "CTL c7"="#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_CTL_shanji.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge-group
shanji <- as.data.frame(pseudotime(BC_CTL_cds)[colnames(BC_CTL_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_CTL_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CTL c1"="#F8E71C",
                             "CTL c2"="#D4E423",
                             "CTL c3"="#B0E02B",
                             "CTL c4"="#8DDD32",
                             "CTL c5"="#69DA39",
                             "CTL c6"="#45D641",
                             "CTL c7"="#21D348"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_CTL_shanji.pdf"
       ,width = 5, height = 3)
dev.off()




plot_cells(ccRCC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(PCa_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()



plot_cells(ccRCC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_CytoTRACE2_Relative_legend.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_CytoTRACE2_Relative_legend.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_CTL_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_CytoTRACE2_Relative_legend.pdf"
       ,width = 8, height = 8)
dev.off()




plot_cells(ccRCC_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_ccRCC_Pseutime_legend.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_PCa_Pseutime_legend.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_CTL_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\2_BC_Pseutime_legend.pdf"
       ,width = 8, height = 8)
dev.off()
#cyto-boxplot
ggplot(BC_CTL@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
ggplot(ccRCC_CTL@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
ggplot(PCa_CTL@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
##T cell gene cluster function ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
genes <- row.names(subset(PCa_CTL_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(PCa_CTL_cds)[match(genes,rownames(rowData(PCa_CTL_cds))),order(pseudotime(PCa_CTL_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\PCa_CTL_clu_df.rds")
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
print(htkm_order)#2_PCa_CTL_Pseutime_genes 4x6
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
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\PCa_CTL_cluster1.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)

##The name of the pathway is too long, so I chose the first five words of the pathway as its name
for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}


ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#CTL_PCa_cluster1 6x8

genes <- row.names(subset(ccRCC_CTL_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(ccRCC_CTL_cds)[match(genes,rownames(rowData(ccRCC_CTL_cds))),order(pseudotime(ccRCC_CTL_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\ccRCC_CTL_clu_df.rds")
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
print(htkm_order)#2_ccRCC_CTL_Pseutime_genes 4x6

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
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\ccRCC_CTL_cluster1.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again
ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)

for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}

ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))#
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#CTL_ccRCC_cluster1 6x8

###BC
genes <- row.names(subset(BC_CTL_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(BC_CTL_cds)[match(genes,rownames(rowData(BC_CTL_cds))),order(pseudotime(BC_CTL_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\BC_CTL_clu_df.rds")
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
print(htkm_order)#2_BC_CTL_Pseutime_genes 4x6

##cluster4
gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\T_cell\\CTL_BC_cluster4.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)

for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}

ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#CTL_BC_cluster4 6x8
save.image("D:/Urinary system tumors/work/3_Pseutime/2_T_cell_monocle3.RData")
#Macrphage monocle3 ----
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Other_cell_monocle3.RData")
#rm(list=setdiff(ls(), "Myeloid"))
##ccRCC
Macrophages <-  subset(Myeloid, Myeloid_type %in% c("Macrophages"))
ggplot(Myeloid@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()

ggplot(Myeloid@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
Myeloid@meta.data$monocle3_cluster <- Myeloid@meta.data$cluster
##ccRCC
ccRCC_Macr <-  subset(Myeloid, tumor_type %in% c("ccRCC"))
ccRCC_Macr <-  subset(ccRCC_Macr, Myeloid_type %in% c("Macrophages"))
ggplot(ccRCC_Macr@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()
ccRCC_Macr <- ccRCC_Macr %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(ccRCC_Macr,group.by = "monocle3_cluster",label = T,label.size = 5)

data <- GetAssayData(ccRCC_Macr, assay = 'RNA', slot = 'counts')
cell_metadata <- ccRCC_Macr@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
ccRCC_Macr_cds <- new_cell_data_set(data,
                                    cell_metadata = cell_metadata,
                                    gene_metadata = gene_annotation)
ccRCC_Macr_cds <- preprocess_cds(ccRCC_Macr_cds, num_dim = 17)
plot_pc_variance_explained(ccRCC_Macr_cds)
ccRCC_Macr_cds <- reduce_dimension(ccRCC_Macr_cds,reduction_method = 'UMAP')
cds.embed <- ccRCC_Macr_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ccRCC_Macr, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
ccRCC_Macr_cds@int_colData$reducedDims$UMAP <- int.embed
ccRCC_Macr_cds <- cluster_cells(ccRCC_Macr_cds,reduction_method = "UMAP")

plot_cells(ccRCC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")


#ccRCC_Macr_cds <- learn_graph(ccRCC_Macr_cds,learn_graph_control=list(ncenter=500),use_partition = TRUE)#Control the number of branch nodes
ccRCC_Macr_cds <- learn_graph(ccRCC_Macr_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(ccRCC_Macr_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(ccRCC_Macr_cds, color_cells_by = "monocle3_cluster")
ccRCC_Macr_cds <- order_cells(ccRCC_Macr_cds,reduction_method = "UMAP")
plot_cells(ccRCC_Macr_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
##
##PCa
PCa_Macr <-  subset(Myeloid, tumor_type %in% c("PCa"))
PCa_Macr <-  subset(PCa_Macr, Myeloid_type %in% c("Macrophages"))

PCa_Macr <- PCa_Macr %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(PCa_Macr,group.by = "monocle3_cluster",label = T,label.size = 5)

FeaturePlot(PCa_Macr,features = "CytoTRACE2_Relative")
DimPlot(PCa_Macr,group.by = "cluster")
data <- GetAssayData(PCa_Macr, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa_Macr@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
PCa_Macr_cds <- new_cell_data_set(data,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
PCa_Macr_cds <- preprocess_cds(PCa_Macr_cds, num_dim = 17)
plot_pc_variance_explained(PCa_Macr_cds)
PCa_Macr_cds <- reduce_dimension(PCa_Macr_cds,reduction_method = 'UMAP')
cds.embed <- PCa_Macr_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(PCa_Macr, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
PCa_Macr_cds@int_colData$reducedDims$UMAP <- int.embed
PCa_Macr_cds <- cluster_cells(PCa_Macr_cds,reduction_method = "UMAP")

plot_cells(PCa_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(PCa_Macr_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


PCa_Macr_cds <- learn_graph(PCa_Macr_cds,learn_graph_control=list(ncenter=100),use_partition = TRUE)#Control the number of branch nodes
#PCa_Macr_cds <- learn_graph(PCa_Macr_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(PCa_Macr_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(PCa_Macr_cds, color_cells_by = "monocle3_cluster")
PCa_Macr_cds <- order_cells(PCa_Macr_cds,reduction_method = "UMAP")
plot_cells(PCa_Macr_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

##BC
BC_Macr <-  subset(Myeloid, tumor_type %in% c("BC"))
BC_Macr <-  subset(BC_Macr, Myeloid_type %in% c("Macrophages"))

BC_Macr <- BC_Macr %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(BC_Macr,group.by = "monocle3_cluster",label = T,label.size = 5)

FeaturePlot(BC_Macr,features = "CytoTRACE2_Relative")
DimPlot(BC_Macr,group.by = "cluster")
data <- GetAssayData(BC_Macr, assay = 'RNA', slot = 'counts')
cell_metadata <- BC_Macr@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
BC_Macr_cds <- new_cell_data_set(data,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_annotation)
BC_Macr_cds <- preprocess_cds(BC_Macr_cds, num_dim = 17)
plot_pc_variance_explained(BC_Macr_cds)
BC_Macr_cds <- reduce_dimension(BC_Macr_cds,reduction_method = 'UMAP')
cds.embed <- BC_Macr_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(BC_Macr, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
BC_Macr_cds@int_colData$reducedDims$UMAP <- int.embed
BC_Macr_cds <- cluster_cells(BC_Macr_cds,reduction_method = "UMAP")

plot_cells(BC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(BC_Macr_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


BC_Macr_cds <- learn_graph(BC_Macr_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#BC_Macr_cds <- learn_graph(BC_Macr_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(BC_Macr_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(BC_Macr_cds, color_cells_by = "monocle3_cluster")
BC_Macr_cds <- order_cells(BC_Macr_cds,reduction_method = "UMAP")
plot_cells(BC_Macr_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

rm(list=setdiff(ls(), c("ccRCC_Macr_cds","ccRCC_Macr",
                        "BC_Macr_cds","BC_Macr",
                        "PCa_Macr_cds","PCa_Macr") ))
##
ccRCC_Macr_modulated_genes <- graph_test(ccRCC_Macr_cds, neighbor_graph = "principal_graph", cores = 4)
BC_Macr_modulated_genes <- graph_test(BC_Macr_cds, neighbor_graph = "principal_graph", cores = 4)
PCa_Macr_modulated_genes <- graph_test(PCa_Macr_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(ccRCC_Macr_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_ccRCC_Macr_modulated_genes.txt",quote = F,sep = "\t")
write.table(BC_Macr_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_BC_Macr_modulated_genes.txt",quote = F,sep = "\t")
write.table(PCa_Macr_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_PCa_Macr_modulated_genes.txt",quote = F,sep = "\t")

##ridge plot ----
shanji <- as.data.frame(pseudotime(PCa_Macr_cds)[colnames(PCa_Macr_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_Macr_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CD83 Macrophages c1" = "#7ED321",
                             "CD83 Macrophages c2" = "#BBDD1F",
                             "CD83 Macrophages c3" = "#F8E71C",
                             "cDC2" = "#50E3C2",
                             "cDC2" = "#84E6A4",
                             "Macrophages c1"="#4A90E2",
                             "Macrophages c10"="#6995F8",
                             "Macrophages c2"="#6783F9",
                             "Macrophages c3"="#6570FA",
                             "Macrophages c4"="#635DFB",
                             "Macrophages c5"="#614BFB",
                             "Macrophages c6"="#5F38FC",
                             "Macrophages c7"="#5F38FC",
                             "Macrophages c8"="#5D25FD",
                             "Macrophages c9"="#3700FF",
                             "Monocytes c1"="#F7C720",
                             "Monocytes c2"="#F5A623",
                             "pDCs"="#B8E986"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_PCa_Macr_shanji.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge-group
shanji <- as.data.frame(pseudotime(ccRCC_Macr_cds)[colnames(ccRCC_Macr_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Macr_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CD83 Macrophages c1" = "#7ED321",
                             "CD83 Macrophages c2" = "#BBDD1F",
                             "CD83 Macrophages c3" = "#F8E71C",
                             "cDC2" = "#50E3C2",
                             "cDC2" = "#84E6A4",
                             "Macrophages c1"="#4A90E2",
                             "Macrophages c10"="#6995F8",
                             "Macrophages c2"="#6783F9",
                             "Macrophages c3"="#6570FA",
                             "Macrophages c4"="#635DFB",
                             "Macrophages c5"="#614BFB",
                             "Macrophages c6"="#5F38FC",
                             "Macrophages c7"="#5F38FC",
                             "Macrophages c8"="#5D25FD",
                             "Macrophages c9"="#3700FF",
                             "Monocytes c1"="#F7C720",
                             "Monocytes c2"="#F5A623",
                             "pDCs"="#B8E986"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_ccRCC_Macr_shanji.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge-group
shanji <- as.data.frame(pseudotime(BC_Macr_cds)[colnames(BC_Macr_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Macr_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=monocle3_cluster,fill=monocle3_cluster))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("CD83 Macrophages c1" = "#7ED321",
                             "CD83 Macrophages c2" = "#BBDD1F",
                             "CD83 Macrophages c3" = "#F8E71C",
                             "cDC2" = "#50E3C2",
                             "cDC2" = "#84E6A4",
                             "Macrophages c1"="#4A90E2",
                             "Macrophages c10"="#6995F8",
                             "Macrophages c2"="#6783F9",
                             "Macrophages c3"="#6570FA",
                             "Macrophages c4"="#635DFB",
                             "Macrophages c5"="#614BFB",
                             "Macrophages c6"="#5F38FC",
                             "Macrophages c7"="#5F38FC",
                             "Macrophages c8"="#5D25FD",
                             "Macrophages c9"="#3700FF",
                             "Monocytes c1"="#F7C720",
                             "Monocytes c2"="#F5A623",
                             "pDCs"="#B8E986"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_BC_Macr_shanji.pdf"
       ,width = 5, height = 3)
dev.off()

plot_cells(ccRCC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_ccRCC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(PCa_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_PCa_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_BC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()



plot_cells(ccRCC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_ccRCC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_PCa_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Macr_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_BC_CytoTRACE2_Relative.pdf"
       ,width = 8, height = 8)
dev.off()




plot_cells(ccRCC_Macr_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_ccRCC_Pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Macr_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_PCa_Pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Macr_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_Macr_BC_Pseutime.pdf"
       ,width = 8, height = 8)
dev.off()
##chrgene_GSEA----
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol"))
genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)
genesets_hallmark2 <- list()
for(i in unique(genesets_hallmark$gs_name)){
  genesets_hallmark2[[i]] <- genesets_hallmark[which(genesets_hallmark$gs_name == i),]$gene_symbol
}
ccRCC_Macr <-  AddModuleScore(ccRCC_Macr,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(ccRCC_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_Macr@meta.data[row.names(x),c(37:dim(ccRCC_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\ccRCC_Macr_hallmarker_pseudotime\\ccRCC_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

PCa_Macr <-  AddModuleScore(PCa_Macr,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(PCa_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_Macr@meta.data[row.names(x),c(37:dim(PCa_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\PCa_Macr_hallmarker_pseudotime\\PCa_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

BC_Macr <-  AddModuleScore(BC_Macr,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(BC_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_Macr@meta.data[row.names(x),c(37:dim(BC_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\BC_Macr_hallmarker_pseudotime\\BC_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

###
geneset <- read.delim("clipboard",header = T)

genesets <- list("pro-inflammatory"=strsplit(geneset$gene[1],' ')[[1]],
                 "Immuneregulatory"=strsplit(geneset$gene[2],' ')[[1]],
                 "Proliferating"=strsplit(geneset$gene[3],' ')[[1]],
                 "Interferonresponsed"=strsplit(geneset$gene[4],' ')[[1]],
                 "Pro-angiogenic"=strsplit(geneset$gene[5],' ')[[1]],
                 "Lipidmetabolism"=strsplit(geneset$gene[6],' ')[[1]],
                 "Extracellular matrix remodelling"=strsplit(geneset$gene[7],' ')[[1]],
                 "Tissue-resident like"=strsplit(geneset$gene[8],' ')[[1]],
                 "M2"=c("MRC1","CD163"),
                 "M1"=c("IL1B","IL6","TNF"))

ccRCC_Macr <-  AddModuleScore(ccRCC_Macr,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(ccRCC_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_Macr@meta.data[row.names(x),c(37:dim(ccRCC_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 51:60){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\ccRCC_Macr_genesets_pseudotime\\ccRCC_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

PCa_Macr <-  AddModuleScore(PCa_Macr,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(PCa_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_Macr@meta.data[row.names(x),c(37:dim(PCa_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 51:60){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\PCa_Macr_genesets_pseudotime\\PCa_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

BC_Macr <-  AddModuleScore(BC_Macr,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(BC_Macr_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_Macr@meta.data[row.names(x),c(37:dim(BC_Macr@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 51:60){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\BC_Macr_genesets_pseudotime\\BC_Macr",colnames(x)[i],".pdf"),width = 4,height = 3)
}

ggplot(BC_Macr@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(ccRCC_Macr@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(PCa_Macr@meta.data,aes(x=cluster,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
##Macrophage gene cluster function ----
genes <- row.names(subset(PCa_Macr_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(PCa_Macr_cds)[match(genes,rownames(rowData(PCa_Macr_cds))),order(pseudotime(PCa_Macr_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\PCa_Macr_clu_df.rds")
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
print(htkm_order)#2_PCa_Macr_Pseutime_genes 4x6
##cluster1

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\PCa_Macr_cluster4.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again


ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Macr_PCa_cluster1 6x8




genes <- row.names(subset(ccRCC_Macr_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(ccRCC_Macr_cds)[match(genes,rownames(rowData(ccRCC_Macr_cds))),order(pseudotime(ccRCC_Macr_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\ccRCC_Macr_clu_df.rds")
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
print(htkm_order)#2_ccRCC_Macr_Pseutime_genes 4x6
##cluster1

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\ccRCC_Macr_cluster4.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Macr_ccRCC_cluster1 6x8

###BC

genes <- row.names(subset(BC_Macr_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(BC_Macr_cds)[match(genes,rownames(rowData(BC_Macr_cds))),order(pseudotime(BC_Macr_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\BC_Macr_clu_df.rds")
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
print(htkm_order)#2_BC_Macr_Pseutime_genes 4x6

##cluster4
gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster3"),1]

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
ego$cluster <- "cluster3"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\Macr_BC_cluster3.txt",sep = "\t",quote = F,row.names = F)
#Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Macr_BC_cluster4 6x8
save.image("D:/Urinary system tumors/work/3_Pseutime/2_Myeloid_monocle3.RData")
#Fibroblast monocle3 ----
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Other_cell_monocle3.RData")
#rm(list=setdiff(ls(), "Fibroblast"))
ggplot(Fibroblast@meta.data,aes(x=Fibroblast_type,y=CytoTRACE2_Relative))+geom_boxplot()
##ccRCC
ccRCC_Fibroblast <-  subset(Fibroblast, tumor_type %in% c("ccRCC"))

ccRCC_Fibroblast <- ccRCC_Fibroblast %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(ccRCC_Fibroblast,group.by = "Fibroblast_type",label = T,label.size = 5)

data <- GetAssayData(ccRCC_Fibroblast, assay = 'RNA', slot = 'counts')
cell_metadata <- ccRCC_Fibroblast@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
ccRCC_Fibroblast_cds <- new_cell_data_set(data,
                                          cell_metadata = cell_metadata,
                                          gene_metadata = gene_annotation)
ccRCC_Fibroblast_cds <- preprocess_cds(ccRCC_Fibroblast_cds, num_dim = 17)
plot_pc_variance_explained(ccRCC_Fibroblast_cds)
ccRCC_Fibroblast_cds <- reduce_dimension(ccRCC_Fibroblast_cds,reduction_method = 'UMAP')
cds.embed <- ccRCC_Fibroblast_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ccRCC_Fibroblast, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
ccRCC_Fibroblast_cds@int_colData$reducedDims$UMAP <- int.embed
ccRCC_Fibroblast_cds <- cluster_cells(ccRCC_Fibroblast_cds,reduction_method = "UMAP")

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")


ccRCC_Fibroblast_cds <- learn_graph(ccRCC_Fibroblast_cds,learn_graph_control=list(ncenter=500),use_partition = TRUE)#Control the number of branch nodes
#ccRCC_Fibroblast_cds <- learn_graph(ccRCC_Fibroblast_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "Fibroblast_type")
ccRCC_Fibroblast_cds <- order_cells(ccRCC_Fibroblast_cds,reduction_method = "UMAP")
plot_cells(ccRCC_Fibroblast_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
##
##PCa
PCa_Fibroblast <-  subset(Fibroblast, tumor_type %in% c("PCa"))

PCa_Fibroblast <- PCa_Fibroblast %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(PCa_Fibroblast,group.by = "Fibroblast_type",label = T,label.size = 5)

FeaturePlot(PCa_Fibroblast,features = "CytoTRACE2_Relative")
DimPlot(PCa_Fibroblast,group.by = "Fibroblast_type")
data <- GetAssayData(PCa_Fibroblast, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa_Fibroblast@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
PCa_Fibroblast_cds <- new_cell_data_set(data,
                                        cell_metadata = cell_metadata,
                                        gene_metadata = gene_annotation)
PCa_Fibroblast_cds <- preprocess_cds(PCa_Fibroblast_cds, num_dim = 17)
plot_pc_variance_explained(PCa_Fibroblast_cds)
PCa_Fibroblast_cds <- reduce_dimension(PCa_Fibroblast_cds,reduction_method = 'UMAP')
cds.embed <- PCa_Fibroblast_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(PCa_Fibroblast, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
PCa_Fibroblast_cds@int_colData$reducedDims$UMAP <- int.embed
PCa_Fibroblast_cds <- cluster_cells(PCa_Fibroblast_cds,reduction_method = "UMAP")

plot_cells(PCa_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(PCa_Fibroblast_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


PCa_Fibroblast_cds <- learn_graph(PCa_Fibroblast_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#PCa_Fibroblast_cds <- learn_graph(PCa_Fibroblast_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(PCa_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(PCa_Fibroblast_cds, color_cells_by = "Fibroblast_type")
PCa_Fibroblast_cds <- order_cells(PCa_Fibroblast_cds,reduction_method = "UMAP")
plot_cells(PCa_Fibroblast_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

##BC
BC_Fibroblast <-  subset(Fibroblast, tumor_type %in% c("BC"))

BC_Fibroblast <- BC_Fibroblast %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(BC_Fibroblast,group.by = "Fibroblast_type",label = T,label.size = 5)

FeaturePlot(BC_Fibroblast,features = "CytoTRACE2_Relative")

data <- GetAssayData(BC_Fibroblast, assay = 'RNA', slot = 'counts')
cell_metadata <- BC_Fibroblast@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
BC_Fibroblast_cds <- new_cell_data_set(data,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)
BC_Fibroblast_cds <- preprocess_cds(BC_Fibroblast_cds, num_dim = 17)
plot_pc_variance_explained(BC_Fibroblast_cds)
BC_Fibroblast_cds <- reduce_dimension(BC_Fibroblast_cds,reduction_method = 'UMAP')
cds.embed <- BC_Fibroblast_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(BC_Fibroblast, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
BC_Fibroblast_cds@int_colData$reducedDims$UMAP <- int.embed
BC_Fibroblast_cds <- cluster_cells(BC_Fibroblast_cds,reduction_method = "UMAP")

plot_cells(BC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(BC_Fibroblast_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


BC_Fibroblast_cds <- learn_graph(BC_Fibroblast_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#BC_Fibroblast_cds <- learn_graph(BC_Fibroblast_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(BC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(BC_Fibroblast_cds, color_cells_by = "Fibroblast_type")
BC_Fibroblast_cds <- order_cells(BC_Fibroblast_cds,reduction_method = "UMAP")
plot_cells(BC_Fibroblast_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

rm(list=setdiff(ls(), c("ccRCC_Fibroblast_cds","ccRCC_Fibroblast",
                        "BC_Fibroblast_cds","BC_Fibroblast",
                        "PCa_Fibroblast_cds","PCa_Fibroblast") ))
##
ccRCC_Fibroblast_modulated_genes <- graph_test(ccRCC_Fibroblast_cds, neighbor_graph = "principal_graph", cores = 4)
BC_Fibroblast_modulated_genes <- graph_test(BC_Fibroblast_cds, neighbor_graph = "principal_graph", cores = 4)
PCa_Fibroblast_modulated_genes <- graph_test(PCa_Fibroblast_cds, neighbor_graph = "principal_graph", cores = 4)
write.table(ccRCC_Fibroblast_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_ccRCC_Fibroblast_modulated_genes.txt",quote = F,sep = "\t")
write.table(BC_Fibroblast_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_BC_Fibroblast_modulated_genes.txt",quote = F,sep = "\t")
write.table(PCa_Fibroblast_modulated_genes,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Myeloid\\2_PCa_Fibroblast_modulated_genes.txt",quote = F,sep = "\t")

####
plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                              "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                               "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()
####
plot_cells(BC_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                              "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                               "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()

####
plot_cells(PCa_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                              "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Fibroblast_cds, color_cells_by = "Fibroblast_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("Myofibroblast" = "#F1A890FF", 
                                                               "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_group.pdf"
       ,width = 8, height = 8)
dev.off()

##ridge plot ----
shanji <- as.data.frame(pseudotime(PCa_Fibroblast_cds)[colnames(PCa_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=Fibroblast_type,fill=Fibroblast_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("Myofibroblast" = "#F1A890FF", 
                             "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast.pdf"
       ,width = 5, height = 3)
dev.off()
shanji <- as.data.frame(pseudotime(ccRCC_Fibroblast_cds)[colnames(ccRCC_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=Fibroblast_type,fill=Fibroblast_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("Myofibroblast" = "#F1A890FF", 
                             "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge-group
shanji <- as.data.frame(pseudotime(BC_Fibroblast_cds)[colnames(BC_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=Fibroblast_type,fill=Fibroblast_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("Myofibroblast" = "#F1A890FF", 
                             "Lipofibroblast"= "#9DD5DBFF"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_shanji.pdf"
       ,width = 5, height = 3)
dev.off()


plot_cells(PCa_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()



plot_cells(PCa_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Fibroblast_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_PCa_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_ccRCC_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Fibroblast_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\2_BC_Fibroblast_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

##chrgene_GSEA----
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol"))
genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)
genesets_hallmark2 <- list()
for(i in unique(genesets_hallmark$gs_name)){
  genesets_hallmark2[[i]] <- genesets_hallmark[which(genesets_hallmark$gs_name == i),]$gene_symbol
}

ccRCC_Fibroblast <-  AddModuleScore(ccRCC_Fibroblast,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(ccRCC_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_Fibroblast@meta.data[row.names(x),c(35:dim(ccRCC_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\ccRCC_Fibroblast_hallmarker_pseudotime\\ccRCC_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}
PCa_Fibroblast <-  AddModuleScore(PCa_Fibroblast,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(PCa_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_Fibroblast@meta.data[row.names(x),c(35:dim(PCa_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\PCa_Fibroblast_hallmarker_pseudotime\\PCa_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}

BC_Fibroblast <-  AddModuleScore(BC_Fibroblast,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(BC_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_Fibroblast@meta.data[row.names(x),c(35:dim(BC_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\BC_Fibroblast_hallmarker_pseudotime\\BC_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}

###
geneset <- read.delim("clipboard",header = T)

genesets <- list("Inflammatory fibroblasts"=strsplit(geneset$gene[1],' ')[[1]],
                 "Myofibroblast"=strsplit(geneset$gene[2],' ')[[1]],
                 "Extracellular matrix remodelling fibroblasts"=strsplit(geneset$gene[3],' ')[[1]],
                 "Pro-angiogensis fibroblasts"=strsplit(geneset$gene[4],' ')[[1]],
                 "Pro-epidermal growth fibroblasts"=strsplit(geneset$gene[5],' ')[[1]],
                 "Adipose-derived fibroblasts"=strsplit(geneset$gene[6],' ')[[1]],
                 "Proliferating fibroblasts"=strsplit(geneset$gene[7],' ')[[1]])


ccRCC_Fibroblast <-  AddModuleScore(ccRCC_Fibroblast,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(ccRCC_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_Fibroblast@meta.data[row.names(x),c(35:dim(ccRCC_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:7){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\ccRCC_Fibroblast_genesets_pseudotime\\ccRCC_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}

PCa_Fibroblast <-  AddModuleScore(PCa_Fibroblast,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(PCa_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_Fibroblast@meta.data[row.names(x),c(35:dim(PCa_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:7){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\PCa_Fibroblast_genesets_pseudotime\\PCa_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}

BC_Fibroblast <-  AddModuleScore(BC_Fibroblast,features = genesets,name = names(genesets))

x <- as.data.frame(pseudotime(BC_Fibroblast_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_Fibroblast@meta.data[row.names(x),c(35:dim(BC_Fibroblast@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:7){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\BC_Fibroblast_genesets_pseudotime\\BC_fib",colnames(x)[i],".pdf"),width = 4,height = 3)
}

ggplot(BC_Fibroblast@meta.data,aes(x=Fibroblast_type,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(ccRCC_Fibroblast@meta.data,aes(x=Fibroblast_type,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(PCa_Fibroblast@meta.data,aes(x=Fibroblast_type,y=CytoTRACE2_Relative))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
##Fibroblast gene cluster function ----
genes <- row.names(subset(PCa_Fibroblast_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(PCa_Fibroblast_cds)[match(genes,rownames(rowData(PCa_Fibroblast_cds))),order(pseudotime(PCa_Fibroblast_cds))]
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
print(htkm)#2_PCa_Fibroblast_Pseutime_genes 4x6
rcl.list <- row_order(htkm) 

clu_df <- lapply(names(rcl.list)[2:4], function(i){
  out <- data.frame(GeneID = rownames(pt.matrix[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
})
x <- as.data.frame(c("CCL2"))
colnames(x) <- "GeneID"
x$Cluster <- "cluster1"

clu_df[[4]] <- clu_df[[3]]
clu_df[[3]] <- clu_df[[2]]
clu_df[[2]] <- clu_df[[1]]
clu_df[[1]] <- x

saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\PCa_Fibroblast_clu_df.rds")
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
print(htkm_order)#2_PCa_Fibroblast_Pseutime_genes 4x6
##cluster4

gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\PCa_Fibroblast_cluster4.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again
library(dplyr)


ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Fibroblast_PCa_cluster4 6x8


genes <- row.names(subset(ccRCC_Fibroblast_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(ccRCC_Fibroblast_cds)[match(genes,rownames(rowData(ccRCC_Fibroblast_cds))),order(pseudotime(ccRCC_Fibroblast_cds))]
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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\ccRCC_clu_df.rds")

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
print(htkm_order)#2_ccRCC_Fibroblast_Pseutime_genes 4x6

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
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\ccRCC_cluster1.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again
ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#ccRCC_cluster1 6x8


genes <- row.names(subset(BC_Fibroblast_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(BC_Fibroblast_cds)[match(genes,rownames(rowData(BC_Fibroblast_cds))),order(pseudotime(BC_Fibroblast_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\BC_clu_df.rds")

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
print(htkm_order)#2_BC_Fibroblast_Pseutime_genes 4x6

##cluster4
gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Fibroblast\\Fibroblast_BC_cluster4.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Fibroblast_BC_cluster4 6x8
save.image("D:/Urinary system tumors/work/3_Pseutime/2_Fibroblast_monocle3.RData")

#Endothelium monocle3 ----
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Other_cell_monocle3.RData")
rm(list=setdiff(ls(), "Endothelium"))

ggplot(Endothelium@meta.data,aes(x=Endothelium_type,y=CytoTRACE2_Relative))+geom_boxplot()

ccRCC_Endothelium <-  subset(Endothelium, tumor_type %in% c("ccRCC"))

ccRCC_Endothelium <- ccRCC_Endothelium %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(ccRCC_Endothelium,group.by = "Endothelium_type",label = T,label.size = 5)

data <- GetAssayData(ccRCC_Endothelium, assay = 'RNA', slot = 'counts')
cell_metadata <- ccRCC_Endothelium@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
ccRCC_Endothelium_cds <- new_cell_data_set(data,
                                           cell_metadata = cell_metadata,
                                           gene_metadata = gene_annotation)
ccRCC_Endothelium_cds <- preprocess_cds(ccRCC_Endothelium_cds, num_dim = 17)
plot_pc_variance_explained(ccRCC_Endothelium_cds)
ccRCC_Endothelium_cds <- reduce_dimension(ccRCC_Endothelium_cds,reduction_method = 'UMAP')
cds.embed <- ccRCC_Endothelium_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ccRCC_Endothelium, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
ccRCC_Endothelium_cds@int_colData$reducedDims$UMAP <- int.embed
ccRCC_Endothelium_cds <- cluster_cells(ccRCC_Endothelium_cds,reduction_method = "UMAP")

plot_cells(ccRCC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")


#ccRCC_Endothelium_cds <- learn_graph(ccRCC_Endothelium_cds,learn_graph_control=list(ncenter=500),use_partition = TRUE)#Control the number of branch nodes
ccRCC_Endothelium_cds <- learn_graph(ccRCC_Endothelium_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(ccRCC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(ccRCC_Endothelium_cds, color_cells_by = "Endothelium_type")
ccRCC_Endothelium_cds <- order_cells(ccRCC_Endothelium_cds,reduction_method = "UMAP")
plot_cells(ccRCC_Endothelium_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)
##
##PCa
PCa_Endothelium <-  subset(Endothelium, tumor_type %in% c("PCa"))

PCa_Endothelium <- PCa_Endothelium %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(PCa_Endothelium,group.by = "Endothelium_type",label = T,label.size = 5)

FeaturePlot(PCa_Endothelium,features = "CytoTRACE2_Relative")
DimPlot(PCa_Endothelium,group.by = "Endothelium_type")
data <- GetAssayData(PCa_Endothelium, assay = 'RNA', slot = 'counts')
cell_metadata <- PCa_Endothelium@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
PCa_Endothelium_cds <- new_cell_data_set(data,
                                         cell_metadata = cell_metadata,
                                         gene_metadata = gene_annotation)
PCa_Endothelium_cds <- preprocess_cds(PCa_Endothelium_cds, num_dim = 17)
plot_pc_variance_explained(PCa_Endothelium_cds)
PCa_Endothelium_cds <- reduce_dimension(PCa_Endothelium_cds,reduction_method = 'UMAP')
cds.embed <- PCa_Endothelium_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(PCa_Endothelium, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
PCa_Endothelium_cds@int_colData$reducedDims$UMAP <- int.embed
PCa_Endothelium_cds <- cluster_cells(PCa_Endothelium_cds,reduction_method = "UMAP")

plot_cells(PCa_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(PCa_Endothelium_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


PCa_Endothelium_cds <- learn_graph(PCa_Endothelium_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#PCa_Endothelium_cds <- learn_graph(PCa_Endothelium_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(PCa_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(PCa_Endothelium_cds, color_cells_by = "Endothelium_type")
PCa_Endothelium_cds <- order_cells(PCa_Endothelium_cds,reduction_method = "UMAP")
plot_cells(PCa_Endothelium_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

##BC
BC_Endothelium <-  subset(Endothelium, tumor_type %in% c("BC"))

BC_Endothelium <- BC_Endothelium %>% RunTSNE(dims = 1:17) %>% RunUMAP(dims = 1:17)
DimPlot(BC_Endothelium,group.by = "Endothelium_type",label = T,label.size = 5)

FeaturePlot(BC_Endothelium,features = "CytoTRACE2_Relative")

data <- GetAssayData(BC_Endothelium, assay = 'RNA', slot = 'counts')
cell_metadata <- BC_Endothelium@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
BC_Endothelium_cds <- new_cell_data_set(data,
                                        cell_metadata = cell_metadata,
                                        gene_metadata = gene_annotation)
BC_Endothelium_cds <- preprocess_cds(BC_Endothelium_cds, num_dim = 17)
plot_pc_variance_explained(BC_Endothelium_cds)
BC_Endothelium_cds <- reduce_dimension(BC_Endothelium_cds,reduction_method = 'UMAP')
cds.embed <- BC_Endothelium_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(BC_Endothelium, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
BC_Endothelium_cds@int_colData$reducedDims$UMAP <- int.embed
BC_Endothelium_cds <- cluster_cells(BC_Endothelium_cds,reduction_method = "UMAP")

plot_cells(BC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 0.5)+theme(legend.position = "right")
plot_cells(BC_Endothelium_cds, color_cells_by = "cluster",cell_size = 0.5)+theme(legend.position = "right")


BC_Endothelium_cds <- learn_graph(BC_Endothelium_cds,learn_graph_control=list(ncenter=1000),use_partition = TRUE)#Control the number of branch nodes
#BC_Endothelium_cds <- learn_graph(BC_Endothelium_cds,use_partition = TRUE)#Control the number of branch nodes
plot_cells(BC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative")
plot_cells(BC_Endothelium_cds, color_cells_by = "Endothelium_type")
BC_Endothelium_cds <- order_cells(BC_Endothelium_cds,reduction_method = "UMAP")
plot_cells(BC_Endothelium_cds, reduction_method = "UMAP", color_cells_by = 'pseudotime', label_cell_groups=FALSE)

rm(list=setdiff(ls(), c("ccRCC_Endothelium_cds","ccRCC_Endothelium",
                        "BC_Endothelium_cds","BC_Endothelium",
                        "PCa_Endothelium_cds","PCa_Endothelium") ))

plot_cells(ccRCC_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("LECs" = "#E38623",
                                                              "VECs" = "#51BDB5",
                                                              "CapECs"="#EDE7BB",
                                                              "AECs"="#BC1A29",
                                                              "ProliferatingECs_c1"="#443972",
                                                              "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Endothelium_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(ccRCC_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("LECs" = "#E38623",
                                                               "VECs" = "#51BDB5",
                                                               "CapECs"="#EDE7BB",
                                                               "AECs"="#BC1A29",
                                                               "ProliferatingECs_c1"="#443972",
                                                               "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()
plot_cells(BC_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("LECs" = "#E38623",
                                                              "VECs" = "#51BDB5",
                                                              "CapECs"="#EDE7BB",
                                                              "AECs"="#BC1A29",
                                                              "ProliferatingECs_c1"="#443972",
                                                              "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Endothelium_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(BC_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("LECs" = "#E38623",
                                                               "VECs" = "#51BDB5",
                                                               "CapECs"="#EDE7BB",
                                                               "AECs"="#BC1A29",
                                                               "ProliferatingECs_c1"="#443972",
                                                               "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+scale_color_manual(values=c("LECs" = "#E38623",
                                                              "VECs" = "#51BDB5",
                                                              "CapECs"="#EDE7BB",
                                                              "AECs"="#BC1A29",
                                                              "ProliferatingECs_c1"="#443972",
                                                              "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Endothelium_cds, color_cells_by = "pseudotime",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium_pseutime.pdf"
       ,width = 8, height = 8)
dev.off()

plot_cells(PCa_Endothelium_cds, color_cells_by = "Endothelium_type",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "right")+scale_color_manual(values=c("LECs" = "#E38623",
                                                               "VECs" = "#51BDB5",
                                                               "CapECs"="#EDE7BB",
                                                               "AECs"="#BC1A29",
                                                               "ProliferatingECs_c1"="#443972",
                                                               "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium_group.pdf"
       ,width = 8, height = 8)
dev.off()

##ridge plot ----
shanji <- as.data.frame(pseudotime(PCa_Endothelium_cds)[colnames(PCa_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
ggplot(pdata, aes(x=pseudotime,y=Endothelium_type,fill=Endothelium_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("LECs" = "#E38623",
                             "VECs" = "#51BDB5",
                             "CapECs"="#EDE7BB",
                             "AECs"="#BC1A29",
                             "ProliferatingECs_c1"="#443972",
                             "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium.pdf"
       ,width = 5, height = 3)
dev.off()

shanji <- as.data.frame(pseudotime(ccRCC_Endothelium_cds)[colnames(ccRCC_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
ggplot(pdata, aes(x=pseudotime,y=Endothelium_type,fill=Endothelium_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("LECs" = "#E38623",
                             "VECs" = "#51BDB5",
                             "CapECs"="#EDE7BB",
                             "AECs"="#BC1A29",
                             "ProliferatingECs_c1"="#443972",
                             "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium.pdf"
       ,width = 5, height = 3)
dev.off()

##ridge group
shanji <- as.data.frame(pseudotime(BC_Endothelium_cds)[colnames(BC_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=Endothelium_type,fill=Endothelium_type))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("LECs" = "#E38623",
                             "VECs" = "#51BDB5",
                             "CapECs"="#EDE7BB",
                             "AECs"="#BC1A29",
                             "ProliferatingECs_c1"="#443972",
                             "ProliferatingECs_c2"="#377DB8"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_shanji.pdf"
       ,width = 5, height = 3)
dev.off()



plot_cells(PCa_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(ccRCC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(BC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  theme(legend.position = "none")+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()




plot_cells(PCa_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_PCa_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(ccRCC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_ccRCC_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()


plot_cells(BC_Endothelium_cds, color_cells_by = "CytoTRACE2_Relative",cell_size = 1.5,trajectory_graph_color = "black",
           label_cell_groups = F,label_leaves = F,label_branch_points = F,label_groups_by_cluster = F)+ 
  scale_color_gradient(low = "#000046",high = "#5AF9DB")
ggsave(filename = "D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\2_BC_Endothelium_CytoTRACE2.pdf"
       ,width = 8, height = 8)
dev.off()


##chrgene_GSEA----
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol"))
genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)
genesets_hallmark2 <- list()
for(i in unique(genesets_hallmark$gs_name)){
  genesets_hallmark2[[i]] <- genesets_hallmark[which(genesets_hallmark$gs_name == i),]$gene_symbol
}

ccRCC_Endothelium <-  AddModuleScore(ccRCC_Endothelium,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(ccRCC_Endothelium_cds))
colnames(x) <-"pseudotime"
x <- cbind(ccRCC_Endothelium@meta.data[row.names(x),c(35:dim(ccRCC_Endothelium@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\ccRCC_Endothelium_hallmarker_pseudotime\\ccRCC_Endo",colnames(x)[i],".pdf"),width = 4,height = 3)
}

PCa_Endothelium <-  AddModuleScore(PCa_Endothelium,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(PCa_Endothelium_cds))
colnames(x) <-"pseudotime"
x <- cbind(PCa_Endothelium@meta.data[row.names(x),c(35:dim(PCa_Endothelium@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\PCa_Endothelium_hallmarker_pseudotime\\PCa_Endo",colnames(x)[i],".pdf"),width = 4,height = 3)
}

BC_Endothelium <-  AddModuleScore(BC_Endothelium,features = genesets_hallmark2,name = names(genesets_hallmark2))

x <- as.data.frame(pseudotime(BC_Endothelium_cds))
colnames(x) <-"pseudotime"
x <- cbind(BC_Endothelium@meta.data[row.names(x),c(35:dim(BC_Endothelium@meta.data)[2],6)],x)

#x <- as.data.frame(cbind(x,t(TCells@assays$RNA$data[trackgene,row.names(x)])))

ggplot(data = x, aes(x = pseudotime, y = x[,2])) +
  theme_classic()+
  geom_smooth()

for(i in 1:50){
  ggplot(data = x, aes(x = pseudotime, y = x[,i])) +
    theme_classic()+
    geom_smooth()+ ggtitle(colnames(x)[i])
  ggsave(paste0("D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\BC_Endothelium_hallmarker_pseudotime\\BC_Endo",colnames(x)[i],".pdf"),width = 4,height = 3)
}

ggplot(BC_Endothelium@meta.data,aes(x=Endothelium_type,y=CytoTRACE2_Relative))+geom_boxplot()
ggplot(ccRCC_Endothelium@meta.data,aes(x=Endothelium_type,y=CytoTRACE2_Relative))+geom_boxplot()
ggplot(PCa_Endothelium@meta.data,aes(x=Endothelium_type,y=CytoTRACE2_Relative))+geom_boxplot()
##Endothelium gene cluster function ----
genes <- row.names(subset(PCa_Endothelium_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(PCa_Endothelium_cds)[match(genes,rownames(rowData(PCa_Endothelium_cds))),order(pseudotime(PCa_Endothelium_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\PCa_Endothelium_clu_df.rds")
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
print(htkm_order)#2_PCa_Endothelium_Pseutime_genes 4x6

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
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\PCa_Endothelium_cluster1.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Endothelium_PCa_cluster1 6x8





genes <- row.names(subset(ccRCC_Endothelium_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(ccRCC_Endothelium_cds)[match(genes,rownames(rowData(ccRCC_Endothelium_cds))),order(pseudotime(ccRCC_Endothelium_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\ccRCC_Endothelium_clu_df.rds")
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
print(htkm_order)#2_ccRCC_Endothelium_Pseutime_genes 4x6

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
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\ccRCC_Endothelium_cluster1.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

#
ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Endothelium_ccRCC_cluster1 6x8

###BC
genes <- row.names(subset(BC_Endothelium_modulated_genes, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(BC_Endothelium_cds)[match(genes,rownames(rowData(BC_Endothelium_cds))),order(pseudotime(BC_Endothelium_cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

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
saveRDS(clu_df,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\BC_Endothelium_clu_df.rds")
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
print(htkm_order)#2_BC_Endothelium_Pseutime_genes 4x6

##cluster4
gene_ID <- kmean_gene[which(kmean_gene$Cluster == "cluster4"),1]

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
ego$cluster <- "cluster4"
write.table(ego,"D:\\Urinary system tumors\\work\\3_Pseutime\\2_others\\Endothelium\\Endothelium_BC_cluster4.txt",sep = "\t",quote = F,row.names = F)
##Change the cluster to cluster2, cluster3, cluster4 in sequence, and then run this code again

#
ego_filter <- by(ego, ego$ONTOLOGY, function(x) arrange(x, x[,8])[1:10,])
ego <- rbind(ego_filter$BP,ego_filter$CC,ego_filter$MF)


for(i in 1:nrow(ego )){
  description_splite=strsplit(ego$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  ego$Description[i]=description_collapse
  ego$Description=gsub(pattern = "NA","",ego$Description)
}



ego$type_order=factor(rev(1:length(rev(ego$Description))),labels=rev(ego$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=ego, aes(x=type_order,y=Count, fill=ONTOLOGY)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) +
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()#Endothelium_BC_cluster4 6x8
save.image("D:/Urinary system tumors/work/3_Pseutime/2_Endothelium_monocle3.RData")


#Ecotype_function Heatmap----
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GOplot)
library(ggplot2)
library(reshape2)
library(pheatmap)
load("D:/Urinary system tumors/work/5_Ecotype/PCa_Ecotypes.RData")
load("D:/Urinary system tumors/work/5_Ecotype/ccRCC_Ecotypes.RData")
load("D:/Urinary system tumors/work/5_Ecotype/BC_Ecotypes.RData")
setwd("D:/Urinary system tumors/work/7_Ecotype_function")
for(i in names(PCa_Ecotypes)){
  for(z in names(PCa_Ecotypes[[i]])){
    gene_ID <- PCa_Ecotypes[[i]][[z]][["gene"]]
    gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                    toType = "ENTREZID" , 
                    OrgDb = org.Hs.eg.db)
    enrich.go <- enrichGO(gene = gene.df$ENTREZID,  
                          OrgDb = org.Hs.eg.db,  
                          keyType = 'ENTREZID',  
                          ont = 'ALL',  
                          pAdjustMethod = 'BH',  
                          pvalueCutoff = 0.05,  
                          qvalueCutoff = 0.05,  
                          readable = TRUE)
    ego <- as.data.frame(enrich.go)
    PCa_Ecotypes[[i]][[z]][["clusterProfilerGO"]] <- ego
  }
}

for(i in names(ccRCC_Ecotypes)){
  for(z in names(ccRCC_Ecotypes[[i]])){
    gene_ID <- ccRCC_Ecotypes[[i]][[z]][["gene"]]
    gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                    toType = "ENTREZID" , 
                    OrgDb = org.Hs.eg.db)
    enrich.go <- enrichGO(gene = gene.df$ENTREZID,  
                          OrgDb = org.Hs.eg.db,  
                          keyType = 'ENTREZID',  
                          ont = 'ALL',  
                          pAdjustMethod = 'BH',  
                          pvalueCutoff = 0.05,  
                          qvalueCutoff = 0.05,  
                          readable = TRUE)
    ego <- as.data.frame(enrich.go)
    ccRCC_Ecotypes[[i]][[z]][["clusterProfilerGO"]] <- ego
  }
}

for(i in names(BC_Ecotypes)){
  for(z in names(BC_Ecotypes[[i]])){
    gene_ID <- BC_Ecotypes[[i]][[z]][["gene"]]
    gene.df <- bitr(gene_ID, fromType = "SYMBOL", 
                    toType = "ENTREZID" , 
                    OrgDb = org.Hs.eg.db)
    enrich.go <- enrichGO(gene = gene.df$ENTREZID,  
                          OrgDb = org.Hs.eg.db,  
                          keyType = 'ENTREZID',  
                          ont = 'ALL',  
                          pAdjustMethod = 'BH',  
                          pvalueCutoff = 0.05,  
                          qvalueCutoff = 0.05,  
                          readable = TRUE)
    ego <- as.data.frame(enrich.go)
    BC_Ecotypes[[i]][[z]][["clusterProfilerGO"]] <- ego
  }
}

PCa_function <- c()
for(i in names(PCa_Ecotypes)){
  for(z in names(PCa_Ecotypes[[i]])){
    temp <- PCa_Ecotypes[[i]][[z]][["clusterProfilerGO"]]
    if(dim(temp)[1] > 0){
      temp$celltype <- z
      temp$Ecotypes <- i 
    }
    PCa_function <- rbind(PCa_function, temp)
  }
}

ccRCC_function <- c()
for(i in names(ccRCC_Ecotypes)){
  for(z in names(ccRCC_Ecotypes[[i]])){
    temp <- ccRCC_Ecotypes[[i]][[z]][["clusterProfilerGO"]]
    if(dim(temp)[1] > 0){
      temp$celltype <- z
      temp$Ecotypes <- i 
    }
    ccRCC_function <- rbind(ccRCC_function, temp)
  }
}
BC_function <- c()
for(i in names(BC_Ecotypes)){
  for(z in names(BC_Ecotypes[[i]])){
    temp <- BC_Ecotypes[[i]][[z]][["clusterProfilerGO"]]
    if(dim(temp)[1] > 0){
      temp$celltype <- z
      temp$Ecotypes <- i 
    }
    BC_function <- rbind(BC_function, temp)
  }
}

write.table(PCa_function, "PCa_Ecotypes_function.txt", quote = F, row.names = F,sep = "\t")
write.table(BC_function, "BC_Ecotypes_function.txt", quote = F, row.names = F,sep = "\t")
write.table(ccRCC_function, "ccRCC_Ecotypes_function.txt", quote = F, row.names = F,sep = "\t")


function_filter <- c("B cell activation", "protein folding", "T cell mediated immunity",
                     "response to endoplasmic reticulum stress", "regulation of cell-cell adhesion",
                     "miRNA transcription", "T cell differentiation", "regulation of T cell activation",
                     "endothelial cell differentiation", "renal system vasculature development", "T cell differentiation",
                     "myeloid cell differentiation", "T cell mediated immunity", "T cell mediated immunity",
                     "TP synthesis coupled electron transport","lysosome organization",
                     "muscle system process", "cytoplasmic translation", "immune response-activating signaling pathway",
                     "chemotaxis", "ATP synthesis coupled electron transport",
                     "MHC protein complex assembly", "regulation of T cell activation", "cell differentiation",
                     "regulation of T cell activation", "positive regulation of inflammatory response",
                     "protein folding", "actin filament bundle organization", "chemotaxis","cytoplasmic translation",
                     "regulation of T cell activation", "regulation of trans-synaptic signaling", "muscle tissue development",
                     "regulation of MAP kinase activity", "endothelial cell migration","ATP biosynthetic process",
                     "T cell receptor signaling pathway","regulation of T cell activation","T cell mediated immunity",
                     "regulation of angiogenesis", "ATP biosynthetic process", "regulation of T cell activation",
                     "regulation of T cell activation",
                     "zinc ion transmembrane transport", "positive regulation of T cell mediated immunity",
                     "antigen processing and presentation of exogenous antigen","cytoplasmic translation",
                     "ribosome assembly", "regulation of execution phase of apoptosis",
                     "cell-cell junction assembly","ATP synthesis coupled electron transport",
                     "regulation of muscle system process")


function_filter <- c("B cell activation","T cell mediated immunity",
                     "response to endoplasmic reticulum stress", "regulation of cell-cell adhesion",
                     "miRNA transcription", "T cell differentiation", "regulation of T cell activation",
                     "endothelial cell differentiation", "renal system vasculature development", "T cell differentiation",
                     "myeloid cell differentiation", "T cell mediated immunity", "T cell mediated immunity",
                     "muscle system process", "immune response-activating signaling pathway",
                     "MHC protein complex assembly", "regulation of T cell activation", "cell differentiation",
                     "regulation of T cell activation", "positive regulation of inflammatory response",
                     "protein folding", "actin filament bundle organization",
                     "regulation of T cell activation", "muscle tissue development",
                     "regulation of MAP kinase activity", "endothelial cell migration",
                     "T cell receptor signaling pathway","regulation of T cell activation","T cell mediated immunity",
                     "regulation of angiogenesis","regulation of T cell activation",
                     "regulation of T cell activation",
                     "positive regulation of T cell mediated immunity",
                     "antigen processing and presentation of exogenous antigen",
                     "ribosome assembly", "regulation of execution phase of apoptosis",
                     "cell-cell junction assembly",
                     "regulation of muscle system process")

function_filter2 <- intersect(intersect(ccRCC_function[which(ccRCC_function$ONTOLOGY == "BP"),]$Description,
                                        PCa_function[which(PCa_function$ONTOLOGY == "BP"),]$Description),
                              BC_function[which(BC_function$ONTOLOGY == "BP"),]$Description)

temp <- PCa_function[which(PCa_function$Description %in% function_filter),]

temp$variable <- paste(temp$Ecotypes, temp$celltype, sep = "_")
temp$p <- -log10(temp$p.adjust)

PCa_mat <- acast(temp, Description ~ variable, value.var = "p")
temp <- BC_function[which(BC_function$Description %in% function_filter),]
temp$variable <- paste(temp$Ecotypes, temp$celltype, sep = "_")
temp$p <- -log10(temp$p.adjust)


BC_mat <- acast(temp, Description ~ variable, value.var = "p")

temp <- ccRCC_function[which(ccRCC_function$Description %in% function_filter),]
temp$variable <- paste(temp$Ecotypes, temp$celltype, sep = "_")
temp$p <- -log10(temp$p.adjust)


ccRCC_mat <- acast(temp, Description ~ variable, value.var = "p")

colnames(PCa_mat) <- paste("PCa",colnames(PCa_mat),sep = "_")
colnames(BC_mat) <- paste("BC",colnames(BC_mat),sep = "_")
colnames(ccRCC_mat) <- paste("ccRCC",colnames(ccRCC_mat),sep = "_")
PCa_mat <- as.data.frame(PCa_mat)
BC_mat <- as.data.frame(BC_mat)
ccRCC_mat <- as.data.frame(ccRCC_mat)
function_filter <- unique(function_filter)
mat <- cbind(PCa_mat[function_filter,], BC_mat[function_filter,])
mat <- cbind(mat[function_filter,], ccRCC_mat[function_filter,])
row.names(mat) <- function_filter

mat[is.na(mat)] <- 0

#labels_row <- replicate(nrow(mat),"")
#labels_row[which(row.names(mat) %in% function_filter)]<- row.names(mat)[which(row.names(mat) %in% function_filter)]

mat <- t(mat)
pheatmap(mat,cluster_rows = T,cluster_cols = T,
         border_color = "NA",color = colorRampPalette(c("white","#F75800","#F75800","#F75800",
                                                        "#F75800","#F75800","#F75800",
                                                        "red" ,"#990033"))(100))#12x24

ggplot(temp,aes(x=Description,y=variable,fill=p))+ #热图绘制
  geom_raster()+scale_fill_gradient(low="#003366", high="#990033")

#Risk gene ----
##PCa Risk gene ----
load("D:/Urinary system tumors/work/5_Ecotype/PCa_Ecotypes.RData")
PCa_sce_anno_second <- readRDS("D:/Urinary system tumors/work/4_merge/PCa_sce_anno_second.rds")
cell <- c()
cell_type <- c("Myeloid","Epithelium","Fibroblast","T_cells","Endothelium","Mast_cells","B_cells" )
for(i in cell_type){
  temp <- read.table(paste0("D:/Urinary system tumors/work/5_Ecotype/output_PCa/",i,"/state_assignment.txt"))
  temp$State <- paste(i,temp$State,sep = "_")
  cell <- rbind(cell, temp)
}
PCa_sce_anno_second@meta.data$State <- NA
PCa_sce_anno_second@meta.data[str_replace_all(cell$ID,"\\.","-"),]$State <- cell$State
PCa_sce_anno_second@meta.data$Ecotypes <- NA
PCa_sce_anno_second@meta.data[which(PCa_sce_anno_second@meta.data$State %in% names(PCa_Ecotypes[["E1"]])),]$Ecotypes <- "E1"
PCa_sce_anno_second@meta.data[which(PCa_sce_anno_second@meta.data$State %in% names(PCa_Ecotypes[["E2"]])),]$Ecotypes <- "E2"
PCa_sce_anno_second@meta.data[which(PCa_sce_anno_second@meta.data$State %in% names(PCa_Ecotypes[["E3"]])),]$Ecotypes <- "E3"

pca <- subset(x = PCa_sce_anno_second, subset = Ecotypes %in% c("E1", "E2", "E3") )
pca <- subset(x = pca, subset = cell_type_second %in% c("Tumor cells"))
snp_gene <- c("TRIM31","FGFRL1","UVSSA","DLG5","SLC44A4","HOXA10")
VlnPlot(pca, features = snp_gene,group.by="Ecotypes")
##ccRCC Risk gene----
load("D:/Urinary system tumors/work/5_Ecotype/ccRCC_Ecotypes.RData")
library(Seurat)
ccRCC_sce_anno_second <- readRDS("D:/Urinary system tumors/work/4_merge/ccRCC_sce_anno_second.rds")
library(stringr)
cell <- c()
cell_type <- c("Myeloid","Epithelium","Fibroblast","T_cells","Endothelium","B_cells" )
for(i in cell_type){
  temp <- read.table(paste0("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/",i,"/state_assignment.txt"))
  temp$State <- paste(i,temp$State,sep = "_")
  cell <- rbind(cell, temp)
}
ccRCC_sce_anno_second@meta.data$State <- NA
ccRCC_sce_anno_second@meta.data[str_replace_all(cell$ID,"\\.","-"),]$State <- cell$State
ccRCC_sce_anno_second@meta.data$Ecotypes <- NA
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E1"]])),]$Ecotypes <- "E1"
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E2"]])),]$Ecotypes <- "E2"
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E3"]])),]$Ecotypes <- "E3"
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E4"]])),]$Ecotypes <- "E4"
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E5"]])),]$Ecotypes <- "E5"
ccRCC_sce_anno_second@meta.data[which(ccRCC_sce_anno_second@meta.data$State %in% names(ccRCC_Ecotypes[["E6"]])),]$Ecotypes <- "E6"
ccRCC_sce_anno_second <- readRDS("D:/Urinary system tumors/work/4_merge/ccRCC_sce_anno_second.rds")
ccRCC <- subset(x = ccRCC_sce_anno_second, subset = Ecotypes %in% c("E1", "E2", "E3", "E4", "E5", "E6") )
ccRCC <- subset(x = ccRCC, subset = cell_type_second %in% c("Tumor cells","Normal epithelium"))
snp_gene <- c("PRKCB","DSP","THRA","FLI1","HOXA10")
VlnPlot(ccRCC, features = snp_gene,group.by="Ecotypes")
VlnPlot(ccRCC, features = snp_gene,group.by="who")
ccRCC2 <- subset(x = ccRCC, subset = Ecotypes %in% c("E2") )
ccRCC2 <- subset(x = ccRCC2, subset = cell_type_second %in% c("Tumor cells"))
VlnPlot(ccRCC2, features = snp_gene,group.by="who")

#ECotype monocle3 ----
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
load("D:/Urinary system tumors/work/5_Ecotype/ccRCC/ccRCC_Ecotypes.RData")
load("D:/Urinary system tumors/work/5_Ecotype/BC/BC_Ecotypes.RData")
load("D:/Urinary system tumors/work/5_Ecotype/PCa/PCa_Ecotypes.RData")
#PCa
PCa_E <- merge(PCa_E1,PCa_E2)
PCa_E <- merge(PCa_E, PCa_E3)

#ccRCC
ccRCC_E <- merge(ccRCC_E1,ccRCC_E2)
ccRCC_E <- merge(ccRCC_E, ccRCC_E3)
ccRCC_E <- merge(ccRCC_E, ccRCC_E4)
ccRCC_E <- merge(ccRCC_E, ccRCC_E5)
ccRCC_E <- merge(ccRCC_E, ccRCC_E6)

#BC
BC_E <- merge(BC_E1,BC_E2)
BC_E <- merge(BC_E, BC_E3)
BC_E <- merge(BC_E, BC_E4)
save(PCa_E,ccRCC_E,BC_E,file = "Ecotypes_merge.RData")

load("D:/Urinary system tumors/work/3_Pseutime/2_Myeloid_monocle3.RData")
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Fibroblast_monocle3.RData")
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_Fibroblast_monocle3.RData")
load("D:/Urinary system tumors/work/3_Pseutime/2_others/2_T_cell_monocle3.RData")
##Macrphage Eco-monocle3 ------
#BC ridge
shanji <- as.data.frame(pseudotime(BC_Macr_cds)[colnames(BC_Macr_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Macr_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,BC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[39] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#3E308F",
                             "E4" = "#81CCD7"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\BC\\BC_Macr.pdf"
       ,width = 5, height = 3)
dev.off()

#ccRCC ridge
shanji <- as.data.frame(pseudotime(ccRCC_Macr_cds)[colnames(ccRCC_Macr_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Macr_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,ccRCC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[39] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4","E5","E6")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E2" = "#53A25D",
                             "E4" = "#ECAE9E",
                             "E5" = "#D1E1A0",
                             "E6" = "#5BB9E4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\ccRCC\\ccRCC_Macr.pdf"
       ,width = 5, height = 3)
dev.off()
##CTL Eco-monocle3 ------
#BC ridge
shanji <- as.data.frame(pseudotime(BC_CTL_cds)[colnames(BC_CTL_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_CTL_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,BC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[39] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#3E308F",
                             "E2" = "#2D7EA1",
                             "E3" = "#3CB371",
                             "E4" = "#81CCD7"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\BC\\BC_CTL.pdf"
       ,width = 5, height = 3)
dev.off()

#ccRCC ridge
shanji <- as.data.frame(pseudotime(ccRCC_CTL_cds)[colnames(ccRCC_CTL_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_CTL_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,ccRCC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[39] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4","E5","E6")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#E5D2DD",
                             "E2" = "#53A25D",
                             "E3" = "#F1BB72",
                             "E4" = "#ECAE9E",
                             "E5" = "#D1E1A0",
                             "E6" = "#5BB9E4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\ccRCC\\ccRCC_CTL.pdf"
       ,width = 5, height = 3)
dev.off()
##Fibroblast Eco-monocle3 ------
#PCa ridge
shanji <- as.data.frame(pseudotime(PCa_Fibroblast_cds)[colnames(PCa_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,PCa_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#84C5ED",
                             "E2" = "#EAAF25",
                             "E3" = "#8975B4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\PCa\\PCa_Fibroblast.pdf"
       ,width = 5, height = 3)
dev.off()

#BC ridge
shanji <- as.data.frame(pseudotime(BC_Fibroblast_cds)[colnames(BC_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,BC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#3E308F",
                             "E2" = "#2D7EA1",
                             "E3" = "#3CB371",
                             "E4" = "#81CCD7"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\BC\\BC_Fibroblast.pdf"
       ,width = 5, height = 3)
dev.off()

#ccRCC ridge
shanji <- as.data.frame(pseudotime(ccRCC_Fibroblast_cds)[colnames(ccRCC_Fibroblast_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Fibroblast_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,ccRCC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4","E5","E6")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#E5D2DD",
                             "E2" = "#53A25D",
                             "E3" = "#F1BB72",
                             "E4" = "#ECAE9E",
                             "E5" = "#D1E1A0",
                             "E6" = "#5BB9E4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\ccRCC\\ccRCC_Fibroblast.pdf"
       ,width = 5, height = 3)
dev.off()
##Endothelium Eco-monocle3 ------
#PCa ridge
shanji <- as.data.frame(pseudotime(PCa_Endothelium_cds)[colnames(PCa_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(PCa_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,PCa_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#84C5ED",
                             "E2" = "#EAAF25",
                             "E3" = "#8975B4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\PCa\\PCa_Endothelium.pdf"
       ,width = 5, height = 3)
dev.off()

#BC ridge
shanji <- as.data.frame(pseudotime(BC_Endothelium_cds)[colnames(BC_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(BC_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,BC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#3E308F",
                             "E2" = "#2D7EA1",
                             "E3" = "#3CB371",
                             "E4" = "#81CCD7"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\BC\\BC_Endothelium.pdf"
       ,width = 5, height = 3)
dev.off()

#ccRCC ridge
shanji <- as.data.frame(pseudotime(ccRCC_Endothelium_cds)[colnames(ccRCC_Endothelium_cds)])
colnames(shanji) <- "pseudotime"
pdata <- as.data.frame(ccRCC_Endothelium_cds@colData[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
pdata <- cbind(pdata,ccRCC_E@meta.data[row.names(pdata),"Ecotypes"])
colnames(pdata)[37] <- "Ecotypes"
pdata <- pdata[which(pdata$Ecotypes %in% c("E1","E2","E3","E4","E5","E6")),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("E1" = "#E5D2DD",
                             "E2" = "#53A25D",
                             "E3" = "#F1BB72",
                             "E4" = "#ECAE9E",
                             "E5" = "#D1E1A0",
                             "E6" = "#5BB9E4"
  ))
ggsave(filename = "D:\\Urinary system tumors\\work\\7_Ecotype_monocle3\\ccRCC\\ccRCC_Endothelium.pdf"
       ,width = 5, height = 3)
dev.off()

load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_PCa_tumor_epi.RData")
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_ccRCC_tumor_epi.RData")

shanji <- as.data.frame(pseudotime(PCa_epi_cds)[colnames(PCa_epi_cds)])
colnames(shanji) <- "pseudotime"
shanji <- cbind(shanji,"time")
shanji <- shanji[intersect(row.names(shanji),colnames(PCa_E) ),]
pdata <- as.data.frame(PCa_E@meta.data[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[which(pdata$cell_type_second == "Tumor cells"),]
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggplot(pdata, aes(x=pseudotime,y=gleason,fill=gleason))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
#ccRCC
shanji <- as.data.frame(pseudotime(ccRCC_epi_cds)[colnames(ccRCC_epi_cds)])
colnames(shanji) <- "pseudotime"
shanji <- cbind(shanji,"time")
x <- ccRCC_E@meta.data[which(ccRCC_E@meta.data$Ecotypes %in%  c("E1","E2","E3","E4","E5","E6") ),]
shanji <- shanji[intersect(row.names(shanji),row.names(x)),]
pdata <- as.data.frame(ccRCC_E@meta.data[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[which(pdata$cell_type_second == "Tumor cells"),]
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
library(ggridges)
library(ggplot2)
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggplot(pdata, aes(x=pseudotime,y=who,fill=who))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())

##Endothelium Eco-monocle3 ------
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_PCa_tumor_epi.RData")
load("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/1_ccRCC_tumor_epi.RData")

shanji <- as.data.frame(pseudotime(PCa_epi_cds)[colnames(PCa_epi_cds)])
colnames(shanji) <- "pseudotime"
shanji <- cbind(shanji,"time")
shanji <- shanji[intersect(row.names(shanji),colnames(PCa_E) ),]
pdata <- as.data.frame(PCa_E@meta.data[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
pdata <- pdata[which(pdata$cell_type_second == "Tumor cells"),]
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggplot(pdata, aes(x=pseudotime,y=gleason,fill=gleason))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
#ccRCC
shanji <- as.data.frame(pseudotime(ccRCC_epi_cds)[colnames(ccRCC_epi_cds)])
colnames(shanji) <- "pseudotime"
shanji <- cbind(shanji,"time")
x <- ccRCC_E@meta.data[which(ccRCC_E@meta.data$Ecotypes %in%  c("E1","E2","E3","E4","E5","E6") ),]
shanji <- shanji[intersect(row.names(shanji),row.names(x)),]
pdata <- as.data.frame(ccRCC_E@meta.data[row.names(shanji),])
pdata$pseudotime <- shanji$pseudotime
ggplot(pdata, aes(x=pseudotime,y=cell_type_second,fill=cell_type_second))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
#pdata <- pdata[-which(pdata$pseudotime == Inf),]
library(ggridges)
library(ggplot2)
pdata <- pdata[which(pdata$cell_type_second == "Tumor cells"),]
ggplot(pdata, aes(x=pseudotime,y=Ecotypes,fill=Ecotypes))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggplot(pdata, aes(x=pseudotime,y=cell_type_second,fill=cell_type_second))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())
ggplot(pdata, aes(x=pseudotime,y=who,fill=who))+
  geom_density_ridges(scale=1)+
  scale_y_discrete("")+
  scale_x_discrete("")+
  theme_minimal()+
  theme(panel.grid = element_blank())