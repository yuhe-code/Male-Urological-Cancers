#cell subset anno ----
##T cell -anno -----
library(Seurat)
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")
sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("T_cells"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"
sce_BC_T <-  subset(sce_BC_tumor, cell_type %in% c("T_cells"))
sce_BC_T@meta.data$tumor_type <- "BC"
sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("T_cells"))
sce_PCa_T@meta.data$tumor_type <- "PCa"
sce_T_all <- merge(sce_ccRCC_T,sce_BC_T)
sce_T_all <- merge(sce_T_all,sce_PCa_T)
table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_T_cells\\0_data\\sce_T_cell_all.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_T_cells\\0_data\\sce_T_cell_all.rds")
Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))

sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(sce, ndims = 50)

pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
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
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1

Idents(sce) <- "RNA_snn_res.0.5"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

DimPlot(sce, label = T, reduction = "tsne",cols = mycolors)
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_T_cells\\0_data\\sce_T_cell_fliter_afterdim.rds")
dim(sce)




sce <- readRDS("D:\\Urinary system tumors\\work\\4_T_cells\\0_data\\sce_T_cell_fliter_afterdim.rds")
T_marker_genes <- c("CCR7","EIF4A3","CCL20",
                    "IL23R","IL1R1","KIT",
                    "NKG7","GNLY","FCGR3A",
                    "PDCD1","CD200","ICOS",
                    "RORA","KLRB1","LTB",
                    "IL2RA","FOXP3","CTLA4",
                    "NKG7","GNLY","FCGR3A",
                    "GZMA","GZMB","CD8A",
                    "ZNF683","ITGA1","CD8A","CD79A","CD4"
)
T_marker_genes <- unique(T_marker_genes)
VlnPlot(sce, features = T_marker_genes,
        group.by = "RNA_snn_res.1",pt.size = 0)


exp_sce_T <- DotPlot(sce,features = T_marker_genes,group.by = "RNA_snn_res.1")$data
min(exp_sce_T$avg.exp.scaled)
max(exp_sce_T$avg.exp.scaled)
#min:-0.6108159  max: 2.041241
DotPlot(sce,features = T_marker_genes,group.by = "RNA_snn_res.0.7")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

sce@meta.data$T_cell_type <- recode(sce@meta.data$RNA_snn_res.0.7,
                                    "0"="CTL",
                                    "1"="CD4Th17",
                                    "2"="CTL",
                                    "3"="NK",
                                    "4"="CD4Th17",
                                    "5"="CTL",
                                    "6"="CD4Treg",
                                    "7"="NK",
                                    "8"="CTL",
                                    "9"="NK",
                                    "10"="CTL",
                                    "11"="CTL",
                                    "12"="CTL",
                                    "13"="CD4Th17",
                                    "14"="NK",
                                    "15"="ILC3",
                                    "16"="CD4Th17",
                                    "17"="B_cells",
                                    "18"="CD4Th17",
                                    "19"="unkwon",
                                    "20"="unkwon",
                                    "21"="unkwon")
DimPlot(sce, label = T, reduction = "umap",group.by = "T_cell_type")

DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.7")

sce_filter_B_cells <- subset(sce,T_cell_type %in% c("B_cells"))
sce_filter <- subset(sce,T_cell_type %in% c("ILC3","CTL","CD4Th17","NK","CD4Treg","unkwon"))

DimPlot(sce_filter, label = T, reduction = "umap",group.by = "T_cell_type")

sce_filter <- NormalizeData(sce_filter) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))

sce_filter <- RunPCA(sce_filter, features = VariableFeatures(object = sce_filter), verbose = F)

pct <- sce_filter [["pca"]]@stdev / sum( sce_filter [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  
pcs <- min(co1, co2)
pcs  
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

T_marker_genes <- c("CD4","IL23R","IL1R1","KIT",
                    "GNLY","FCGR3A",
                    "PDCD1","CD200","ICOS",
                    "RORA","KLRB1","LTB",
                    "IL2RA","FOXP3","CTLA4",
                    "GZMA","GZMB","CD8A"
)
T_marker_genes <- unique(T_marker_genes)
DotPlot(sce_filter,features = T_marker_genes,group.by = "T_cell_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 


sce_filter <- sce_filter %>% RunTSNE(dims = 1:11) %>% RunUMAP(dims = 1:11)

mycolors <- c("ILC3" = "#F3CE33",
              "CTL" = "#E1835E",
              "CD8Trm"="#BA2371",
              "CD4Treg"="#C4D9C5",
              "CD4Th17"= "#48BBB9",
              "NK"="#4B59A6",
              "unkown"="#171831")

DimPlot(sce_filter, label = F,reduction = "umap",group.by = "T_cell_type",repel = T,cols = mycolors)+ 
  theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce_filter,features = T_marker_genes,group.by = "T_cell_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()

saveRDS(sce_filter,"D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno.rds")
saveRDS(sce_filter_B_cells,"D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno_Bcells.rds")


library(Seurat)
library(ks)
library(dplyr)
T_cell <- readRDS("D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno.rds")
T_cell@meta.data$cluster <- NA
T_cell@meta.data$cluster <- recode(T_cell@meta.data$RNA_snn_res.0.7,
                                   "0"="CTL c1",
                                   "1"="CD4Th17 c1",
                                   "2"="CTL c2",
                                   "3"="NK c1",
                                   "4"="CD4Th17 c2",
                                   "5"="CTL c3",
                                   "6"="CD4Treg",
                                   "7"="NK c2",
                                   "8"="CTL c4",
                                   "9"="NK c3",
                                   "10"="CTL c5",
                                   "11"="CTL c6",
                                   "12"="CTL c7",
                                   "13"="CD4Th17 c3",
                                   "14"="NK c4",
                                   "15"="ILC3",
                                   "16"="CD4Th17 c4",
                                   "17"="B_cells",#Identify a small number of B cells
                                   "18"="CD4Th17 c5",
                                   "19"="unkwon",
                                   "20"="unkwon",
                                   "21"="unkwon")
#
saveRDS(T_cell,"D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno.rds")

data <- as.data.frame(cbind(T_cell@reductions[["umap"]]@cell.embeddings,as.character(T_cell@meta.data$cluster),
                            as.character(T_cell@meta.data$T_cell_type)))
colnames(data) <- c("umap_x","umap_y","celltype","group") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)
#
mycolors <- c("ILC3" = "#F3CE33",
              "CTL" = "#E1835E",
              "CD8Trm"="#BA2371",
              "CD4Treg"="#C4D9C5",
              "CD4Th17"= "#48BBB9",
              "NK"="#4B59A6",
              "unkown"="#171831")
##ILC3
ILC3 <- data[which(data$group == "ILC3"),c(1,2)]
colnames(ILC3) <- c("umap_x","umap_y")
ILC3$umap_x <- as.numeric(ILC3$umap_x)
ILC3$umap_y <- as.numeric(ILC3$umap_y)
ILC3 <- kde(x=ILC3)
ILC3 <- with(ILC3, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
##CTL
CTL <- data[which(data$group == "CTL"),c(1,2)]
colnames(CTL) <- c("umap_x","umap_y")
CTL$umap_x <- as.numeric(CTL$umap_x)
CTL$umap_y <- as.numeric(CTL$umap_y)
CTL <- kde(x=CTL)
CTL <- with(CTL, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
##CD4Treg
CD4Treg <- data[which(data$group == "CD4Treg"),c(1,2)]
colnames(CD4Treg) <- c("umap_x","umap_y")
CD4Treg$umap_x <- as.numeric(CD4Treg$umap_x)
CD4Treg$umap_y <- as.numeric(CD4Treg$umap_y)
CD4Treg <- kde(x=CD4Treg)
CD4Treg <- with(CD4Treg, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
##CD4Th17
CD4Th17 <- data[which(data$group == "CD4Th17"),c(1,2)]
colnames(CD4Th17) <- c("umap_x","umap_y")
CD4Th17$umap_x <- as.numeric(CD4Th17$umap_x)
CD4Th17$umap_y <- as.numeric(CD4Th17$umap_y)
CD4Th17 <- kde(x=CD4Th17)
CD4Th17 <- with(CD4Th17, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
##NK
NK <- data[which(data$group == "NK"),c(1,2)]
colnames(NK) <- c("umap_x","umap_y")
NK$umap_x <- as.numeric(NK$umap_x)
NK$umap_y <- as.numeric(NK$umap_y)
NK <- kde(x=NK)
NK <- with(NK, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
color <- c("#50E3C2","#43C9CD","#37AFD8","#2A94E3","#1D7AEE","skyblue",
           "#F8E71C","#D4E423","#B0E02B","#8DDD32","#69DA39","#45D641","#21D348","#F5A623",
           "#9013FE","#793DF5","#6166EB","#4A90E2","grey")

plotdata <- ggplot(data, aes(x = umap_x, 
                             y = umap_y,
                             color = celltype)) +
  geom_point(size = 0.001) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(color = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(values = color)
mycolors <- c("ILC3" = "#F3CE33",
              "CTL" = "#E1835E",
              "CD8Trm"="#BA2371",
              "CD4Treg"="#C4D9C5",
              "CD4Th17"= "#48BBB9",
              "NK"="#4B59A6",
              "unkown"="#171831")

con.size=1.5
cl='black'
plotdata+geom_path(aes(x, y), data=data.frame(ILC3),col="#B59304",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(CTL),col="#AF360C",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(CD4Treg),col="#688269",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(CD4Th17),col="#006B63",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(NK),col="#2E3863",linetype = 2,size=con.size)
ggplot(data, aes(x = umap_x, 
                 y = umap_y,
                 color = celltype)) +
  geom_point(size = 2)+
  scale_color_manual(values = color)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()
##Myeolid -anno -----
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")

sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("Myeloid"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"

sce_BC_T <-  subset(sce_BC_tumor, cell_type %in% c("Myeloid"))
sce_BC_T@meta.data$tumor_type <- "BC"

sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("Myeloid"))
sce_PCa_T@meta.data$tumor_type <- "PCa"

sce_T_all <- merge(sce_ccRCC_T,sce_BC_T)
sce_T_all <- merge(sce_T_all,sce_PCa_T)

table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_Myeloid\\0_data\\sce_Myeloid_all.rds")


sce <- readRDS("D:\\Urinary system tumors\\work\\4_Myeloid\\0_data\\sce_Myeloid_all.rds")
Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))

sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(sce, ndims = 50)

pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #43
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #15
pcs <- min(co1, co2)
pcs  #15
plot_df <-data.frame(pct = pct, cumu = cumu, raMacrophages_CD83 = 1:length(pct))

ggplot(plot_df,aes(cumu, pct, label = raMacrophages_CD83,color = raMacrophages_CD83 > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()
pc_num= 1:17
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1
Idents(sce) <- "RNA_snn_res.0.5"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

DimPlot(sce, label = T, reduction = "tsne",cols = mycolors)
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Myeloid\\0_data\\sce_Myeloid_fliter_afterdim.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_Myeloid\\0_data\\sce_Myeloid_fliter_afterdim.rds")
Myeloid_marker_genes <- c("CD68","CD163","CD14",
                          "LYZ","FCN1","FCGR3A",
                          "CD1C","CLEC10A","FCER1A",
                          "XCR1","CLEC9A",
                          "LAMP3","CCR7","CD83",
                          "IL3RA","LILRA4","GZMB")
Myeloid_marker_genes <- unique(Myeloid_marker_genes)
VlnPlot(sce, features = Myeloid_marker_genes,
        group.by = "RNA_snn_res.0.5",pt.size = 0)

exp_sce_T <- DotPlot(sce,features = Myeloid_marker_genes,group.by = "RNA_snn_res.0.5")$data
min(exp_sce_T$avg.exp.scaled)
max(exp_sce_T$avg.exp.scaled)

DotPlot(sce,features = Myeloid_marker_genes,group.by = "RNA_snn_res.0.3")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blaMacrophages_CD83(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 


sce@meta.data$Myeloid_type <- recode(sce@meta.data$RNA_snn_res.0.3,
                                     "0"="Macrophages",
                                     "1"="Macrophages",
                                     "2"="Macrophages",
                                     "3"="Macrophages",
                                     "4"="cDC2",
                                     "5"="Monocytes",
                                     "6"="Macrophages",
                                     "7"="Macrophages",
                                     "8"="Macrophages",
                                     "9"="Macrophages_CD83",
                                     "10"="Macrophages",
                                     "11"="cDC1",
                                     "12"="Macrophages",
                                     "13"="Macrophages",
                                     "14"="Macrophages_CD83",
                                     "15"="Macrophages_CD83",
                                     "16"="Monocytes",
                                     "17"="pDCs")

mycolors <- c("Monocytes" = "#921F58FF", 
              "cDC2"= "#CE599DFF", 
              "Macrophages"="#F1A890FF", 
              "cDC1"  = "#9DD5DBFF", 
              "Macrophages_CD83" = "#D5DED0FF", 
              "pDCs"="#68C2BFFF")
DimPlot(sce, label = T, reduction = "umap",group.by = "Myeloid_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.3")

DimPlot(sce, label = F,reduction = "umap",group.by = "Myeloid_type",repel = T,cols = mycolors)+ 
  theme(legend.position = "none",axis.title.y=element_blaMacrophages_CD83(), axis.text.y=element_blaMacrophages_CD83(),axis.title.x=element_blaMacrophages_CD83(), axis.text.x=element_blaMacrophages_CD83(),
        axis.line=element_blaMacrophages_CD83(),axis.ticks=element_blaMacrophages_CD83(),plot.title = element_blaMacrophages_CD83() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce,features = Myeloid_marker_genes,group.by = "Myeloid_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blaMacrophages_CD83(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_sce_anno.rds")

Myeloid <- readRDS("D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_sce_anno.rds")
Myeloid@meta.data$cluster <- NA
Myeloid@meta.data$cluster <- recode(Myeloid@meta.data$RNA_snn_res.0.3,
                                    "0"="Macrophages c1",
                                    "1"="Macrophages c2",
                                    "2"="Macrophages c3",
                                    "3"="Macrophages c4",
                                    "4"="cDC2",
                                    "5"="Monocytes c1",
                                    "6"="Macrophages c5",
                                    "7"="Macrophages c6",
                                    "8"="Macrophages c7",
                                    "9"="CD83 Macrophages c1",
                                    "10"="Macrophages c8",
                                    "11"="cDC1",
                                    "12"="Macrophages c9",
                                    "13"="Macrophages c10",
                                    "14"="CD83 Macrophages c2",
                                    "15"="CD83 Macrophages c3",
                                    "16"="Monocytes c2",
                                    "17"="pDCs")
saveRDS(Myeloid,"D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_sce_anno.rds")

data <- as.data.frame(cbind(Myeloid@reductions[["umap"]]@cell.embeddings,as.character(Myeloid@meta.data$cluster),
                            as.character(Myeloid@meta.data$Myeloid_type)))
colnames(data) <- c("umap_x","umap_y","celltype","group") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)
mycolors <- c("Monocytes" = "#921F58FF", 
              "cDC2"= "#CE599DFF", 
              "Macrophages"="#F1A890FF", 
              "cDC1"  = "#9DD5DBFF", 
              "Macrophages_CD83" = "#D5DED0FF", 
              "pDCs"="#68C2BFFF")
####Monocytes
Monocytes <- data[which(data$group == "Monocytes"),c(1,2)]
colnames(Monocytes) <- c("umap_x","umap_y")
Monocytes$umap_x <- as.numeric(Monocytes$umap_x)
Monocytes$umap_y <- as.numeric(Monocytes$umap_y)
Monocytes <- kde(x=Monocytes)
Monocytes <- with(Monocytes, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[2]])
####cDC2
cDC2 <- data[which(data$group == "cDC2"),c(1,2)]
colnames(cDC2) <- c("umap_x","umap_y")
cDC2$umap_x <- as.numeric(cDC2$umap_x)
cDC2$umap_y <- as.numeric(cDC2$umap_y)
cDC2 <- kde(x=cDC2)
cDC2 <- with(cDC2, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####Macrophages
Macrophages <- data[which(data$group == "Macrophages"),c(1,2)]
colnames(Macrophages) <- c("umap_x","umap_y")
Macrophages$umap_x <- as.numeric(Macrophages$umap_x)
Macrophages$umap_y <- as.numeric(Macrophages$umap_y)
Macrophages <- kde(x=Macrophages)
Macrophages <- with(Macrophages, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[3]])
####cDC1
cDC1 <- data[which(data$group == "cDC1"),c(1,2)]
colnames(cDC1) <- c("umap_x","umap_y")
cDC1$umap_x <- as.numeric(cDC1$umap_x)
cDC1$umap_y <- as.numeric(cDC1$umap_y)
cDC1 <- kde(x=cDC1)
cDC1 <- with(cDC1, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[2]])
####Macrophages_CD83
Macrophages_CD83 <- data[which(data$group == "Macrophages_CD83"),c(1,2)]
colnames(Macrophages_CD83) <- c("umap_x","umap_y")
Macrophages_CD83$umap_x <- as.numeric(Macrophages_CD83$umap_x)
Macrophages_CD83$umap_y <- as.numeric(Macrophages_CD83$umap_y)
Macrophages_CD83 <- kde(x=Macrophages_CD83)
Macrophages_CD83 <- with(Macrophages_CD83, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####pDCs
pDCs <- data[which(data$group == "pDCs"),c(1,2)]
colnames(pDCs) <- c("umap_x","umap_y")
pDCs$umap_x <- as.numeric(pDCs$umap_x)
pDCs$umap_y <- as.numeric(pDCs$umap_y)
pDCs <- kde(x=pDCs)
pDCs <- with(pDCs, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])

color <- c("CD83 Macrophages c1" = "#7ED321",
           "CD83 Macrophages c2" = "#BBDD1F",
           "CD83 Macrophages c3" = "#F8E71C",
           "cDC1" = "#50E3C2",
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
           "pDCs"="#B8E986")
plotdata <- ggplot(data, aes(x = umap_x, 
                             y = umap_y,
                             color = celltype)) +
  geom_point(size = 0.001) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(color = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(values = color)
con.size=1.5
cl='black'

mycolors <- c("Monocytes" = "#921F58FF", 
              "cDC2"= "#CE599DFF", 
              "Macrophages"="#F1A890FF", 
              "cDC1"  = "#9DD5DBFF", 
              "Macrophages_CD83" = "#D5DED0FF", 
              "pDCs"="#68C2BFFF")

plotdata+geom_path(aes(x, y), data=data.frame(Monocytes),col="#921F58FF",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(cDC2),col="#CE599DFF",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Macrophages),col="#F1A890FF",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(cDC1),col="red",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Macrophages_CD83),col="#D5DED0FF",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(pDCs),col="skyblue",linetype = 2,size=con.size)


ggplot(data, aes(x = umap_x, 
                 y = umap_y,
                 color = celltype)) +
  geom_point(size = 2)+
  scale_color_manual(values = color)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Myeloid\\2_Myeloid_anno\\Myeloid_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()
##Endo -anno -----
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")


sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("Endothelium"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"

sce_BC_T <-  subset(sce_BC_tumor, cell_type %in% c("Endothelium"))
sce_BC_T@meta.data$tumor_type <- "BC"

sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("Endothelium"))
sce_PCa_T@meta.data$tumor_type <- "PCa"


sce_T_all <- merge(sce_ccRCC_T,sce_BC_T)
sce_T_all <- merge(sce_T_all,sce_PCa_T)

table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_Endothelium\\0_data\\sce_Endothelium_all.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_Endothelium\\0_data\\sce_Endothelium_all.rds")
Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))

sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(sce, ndims = 50)
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #43
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  
pcs <- min(co1, co2)
pcs  
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()
pc_num= 1:17
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1

#tsne&umap（0_data，1_landscape_plot）
Idents(sce) <- "RNA_snn_res.0.3"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

#DimPlot(sce, label = T, reduction = "tsne",cols = mycolors)
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Endothelium\\0_data\\sce_Endothelium_fliter_afterdim.rds")
dim(sce)


sce <- readRDS("D:\\Urinary system tumors\\work\\4_Endothelium\\0_data\\sce_Endothelium_fliter_afterdim.rds")
Endothelium_marker_genes <- c("GJA5","SEMA3G","HEY1",
                              "KDR","PLVAP","RGCC",
                              "ACKR1","SELE","IL1R1",
                              "LYVE1","PROX1","CCL21",
                              "MKI67","TOP2A","STMN1")
Endothelium_marker_genes <- unique(Endothelium_marker_genes)
VlnPlot(sce, features = Endothelium_marker_genes,
        group.by = "RNA_snn_res.0.3",pt.size = 0)

exp_sce_T <- DotPlot(sce,features = Endothelium_marker_genes,group.by = "RNA_snn_res.0.5")$data
min(exp_sce_T$avg.exp.scaled)
max(exp_sce_T$avg.exp.scaled)
#min:-0.6108159  max: 2.041241

DotPlot(sce,features = Endothelium_marker_genes,group.by = "RNA_snn_res.0.3")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 


sce@meta.data$Endothelium_type <- recode(sce@meta.data$RNA_snn_res.0.3,
                                         "0"="CapECs",
                                         "1"="CapECs",
                                         "2"="VECs",
                                         "3"="VECs",
                                         "4"="CapECs",
                                         "5"="CapECs",
                                         "6"="AECs",
                                         "7"="CapECs",
                                         "8"="LECs",
                                         "9"="CapECs",
                                         "10"="CapECs",
                                         "11"="CapECs",
                                         "12"="ProliferatingECs_c1",
                                         "13"="CapECs",
                                         "14"="ProliferatingECs_c2",
                                         "15"="CapECs")

mycolors <- c("LECs" = "#E38623",
              "VECs" = "#51BDB5",
              "CapECs"="#EDE7BB",
              "AECs"="#BC1A29",
              "ProliferatingECs_c1"="#443972",
              "ProliferatingECs_c2"="#377DB8")
DimPlot(sce, label = T, reduction = "umap",group.by = "Endothelium_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.3")


sce_filter <- sce_filter %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)




DimPlot(sce, label = F,reduction = "umap",group.by = "Endothelium_type",repel = T,cols = mycolors)+ 
  theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce,features = Endothelium_marker_genes,group.by = "Endothelium_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_sce_anno.rds")

Endothelium <- readRDS("D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_sce_anno.rds")
data <- as.data.frame(cbind(Endothelium@reductions[["umap"]]@cell.embeddings,as.character(Endothelium@meta.data$Endothelium_type),
                            as.character(Endothelium@meta.data$Endothelium_type)))
colnames(data) <- c("umap_x","umap_y","celltype","group") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)
mycolors <- c("LECs" = "#E38623",
              "VECs" = "#51BDB5",
              "CapECs"="#EDE7BB",
              "AECs"="#BC1A29",
              "ProliferatingECs_c1"="#443972",
              "ProliferatingECs_c2"="#377DB8")
####LECs
LECs <- data[which(data$group == "LECs"),c(1,2)]
colnames(LECs) <- c("umap_x","umap_y")
LECs$umap_x <- as.numeric(LECs$umap_x)
LECs$umap_y <- as.numeric(LECs$umap_y)
LECs <- kde(x=LECs)
LECs <- with(LECs, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####VECs
VECs <- data[which(data$group == "VECs"),c(1,2)]
colnames(VECs) <- c("umap_x","umap_y")
VECs$umap_x <- as.numeric(VECs$umap_x)
VECs$umap_y <- as.numeric(VECs$umap_y)
VECs <- kde(x=VECs)
VECs <- with(VECs, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####CapECs
CapECs <- data[which(data$group == "CapECs"),c(1,2)]
colnames(CapECs) <- c("umap_x","umap_y")
CapECs$umap_x <- as.numeric(CapECs$umap_x)
CapECs$umap_y <- as.numeric(CapECs$umap_y)
CapECs <- kde(x=CapECs)
CapECs <- with(CapECs, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####AECs
AECs <- data[which(data$group == "AECs"),c(1,2)]
colnames(AECs) <- c("umap_x","umap_y")
AECs$umap_x <- as.numeric(AECs$umap_x)
AECs$umap_y <- as.numeric(AECs$umap_y)
AECs <- kde(x=AECs)
AECs <- with(AECs, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[4]])
####ProliferatingECs_c1
ProliferatingECs_c1 <- data[which(data$group == "ProliferatingECs_c1"),c(1,2)]
colnames(ProliferatingECs_c1) <- c("umap_x","umap_y")
ProliferatingECs_c1$umap_x <- as.numeric(ProliferatingECs_c1$umap_x)
ProliferatingECs_c1$umap_y <- as.numeric(ProliferatingECs_c1$umap_y)
ProliferatingECs_c1 <- kde(x=ProliferatingECs_c1)
ProliferatingECs_c1 <- with(ProliferatingECs_c1, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
####ProliferatingECs_c2
ProliferatingECs_c2 <- data[which(data$group == "ProliferatingECs_c2"),c(1,2)]
colnames(ProliferatingECs_c2) <- c("umap_x","umap_y")
ProliferatingECs_c2$umap_x <- as.numeric(ProliferatingECs_c2$umap_x)
ProliferatingECs_c2$umap_y <- as.numeric(ProliferatingECs_c2$umap_y)
ProliferatingECs_c2 <- kde(x=ProliferatingECs_c2)
ProliferatingECs_c2 <- with(ProliferatingECs_c2, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])

plotdata <- ggplot(data, aes(x = umap_x, 
                             y = umap_y,
                             color = celltype)) +
  geom_point(size = 0.001) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(color = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(values = mycolors)

con.size=1.5
cl='black'

mycolors <- c("LECs" = "#E38623",
              "VECs" = "#51BDB5",
              "CapECs"="#EDE7BB",
              "AECs"="#BC1A29",
              "ProliferatingECs_c1"="#443972",
              "ProliferatingECs_c2"="#377DB8")

plotdata+geom_path(aes(x, y), data=data.frame(LECs),col="#7A4515",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(VECs),col="#306863",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(CapECs),col="#9E9A7E",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(AECs),col="#70111D",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(ProliferatingECs_c1),col="#272244",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(ProliferatingECs_c2),col="#1C435E",linetype = 2,size=con.size)

ggplot(data, aes(x = umap_x, 
                 y = umap_y,
                 color = celltype)) +
  geom_point(size = 2)+
  scale_color_manual(values = mycolors)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Endothelium\\2_Endothelium_anno\\Endothelium_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()
##Fibroblast -anno -----
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")


sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("Fibroblast"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"

sce_BC_T <-  subset(sce_BC_tumor, cell_type %in% c("Fibroblast"))
sce_BC_T@meta.data$tumor_type <- "BC"

sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("Fibroblast"))
sce_PCa_T@meta.data$tumor_type <- "PCa"


sce_T_all <- merge(sce_ccRCC_T,sce_BC_T)
sce_T_all <- merge(sce_T_all,sce_PCa_T)

table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_Fibroblast\\0_data\\sce_Fibroblast_all.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_Fibroblast\\0_data\\sce_Fibroblast_all.rds")
Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))

sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(sce, ndims = 50)

pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #44

co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #17
pcs <- min(co1, co2)
pcs  #17
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))


library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()
pc_num= 1:17
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1

Idents(sce) <- "RNA_snn_res.0.3"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

DimPlot(sce, label = T, reduction = "tsne",cols = mycolors)
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Fibroblast\\0_data\\sce_Fibroblast_fliter_afterdim.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_Fibroblast\\0_data\\sce_Fibroblast_fliter_afterdim.rds")
Idents(sce) <- "RNA_snn_res.0.2"
sce.markers <- FindAllMarkers(object = sce,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0.6) 
paste(sce.markers[which(sce.markers$cluster == 0 & sce.markers$pct.1 > 0.6),]$gene,collapse = ",")
Fibroblast_marker_genes <- c("ACTA2","COL1A1","APOE")
DotPlot(sce,features = Fibroblast_marker_genes,group.by = "RNA_snn_res.0.2")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

sce@meta.data$Fibroblast_type <- recode(sce@meta.data$RNA_snn_res.0.2,
                                        "0"="Myofibroblast",
                                        "1"="Myofibroblast",
                                        "2"="Myofibroblast",
                                        "3"="Myofibroblast",
                                        "4"="Myofibroblast",
                                        "5"="Lipofibroblast",
                                        "6"="Myofibroblast",
                                        "7"="Myofibroblast",
                                        "8"="Myofibroblast",
                                        "9"="Lipofibroblast",
                                        "10"="Myofibroblast",
                                        "11"="Myofibroblast",
                                        "12"="Myofibroblast",
                                        "13"="Myofibroblast",
                                        "14"="Lipofibroblast",
                                        "15"="Myofibroblast")
DimPlot(sce, label = T, group.by = "RNA_snn_res.0.2")

mycolors <- c("Myofibroblast" = "#F1A890FF", 
              "Lipofibroblast"= "#9DD5DBFF")
DimPlot(sce, label = F,reduction = "umap",group.by = "Fibroblast_type",repel = T,cols = mycolors)+ 
  theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\Fibroblast_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce,features = Fibroblast_marker_genes,group.by = "Fibroblast_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\Fibroblast_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\Fibroblast_sce_anno.rds")
#####美化
Fibroblast <- readRDS("D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\Fibroblast_sce_anno.rds")

data <- as.data.frame(cbind(Fibroblast@reductions[["umap"]]@cell.embeddings,as.character(Fibroblast@meta.data$Fibroblast_type),
                            as.character(Fibroblast@meta.data$Fibroblast_type)))
colnames(data) <- c("umap_x","umap_y","celltype","group") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)
####Myofibroblast
Myofibroblast <- data[which(data$group == "Myofibroblast"),c(1,2)]
colnames(Myofibroblast) <- c("umap_x","umap_y")
Myofibroblast$umap_x <- as.numeric(Myofibroblast$umap_x)
Myofibroblast$umap_y <- as.numeric(Myofibroblast$umap_y)
Myofibroblast <- kde(x=Myofibroblast)
Myofibroblast <- with(Myofibroblast, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[4]])
####Lipofibroblast
Lipofibroblast <- data[which(data$group == "Lipofibroblast"),c(1,2)]
colnames(Lipofibroblast) <- c("umap_x","umap_y")
Lipofibroblast$umap_x <- as.numeric(Lipofibroblast$umap_x)
Lipofibroblast$umap_y <- as.numeric(Lipofibroblast$umap_y)
Lipofibroblast <- kde(x=Lipofibroblast)
Lipofibroblast <- with(Lipofibroblast, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])

color <- c("#B8E986","#F5A623")
plotdata <- ggplot(data, aes(x = umap_x, 
                             y = umap_y,
                             color = celltype)) +
  geom_point(size = 0.001) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(color = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(values = color)
con.size=1.5
cl='black'

plotdata+geom_path(aes(x, y), data=data.frame(Lipofibroblast),col="#597A34",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myofibroblast),col="#8C601F",linetype = 2,size=con.size)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\美化Fibroblast_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

ggplot(data, aes(x = umap_x, 
                 y = umap_y,
                 color = celltype)) +
  geom_point(size = 2)+
  scale_color_manual(values = color)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Fibroblast\\2_Fibroblast_anno\\美化Fibroblast_umap_cell_type图例.pdf"
       ,width = 5, height = 5)
dev.off()
##Mast cell -anno ----
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")


sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("Mast_cells"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"


sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("Mast_cells"))
sce_PCa_T@meta.data$tumor_type <- "PCa"


sce_T_all <- merge(sce_ccRCC_T,sce_PCa_T)

table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_Mast_cells\\0_data\\sce_Mast_cells_all.rds")

sce <- readRDS("D:\\Urinary system tumors\\work\\4_Mast_cells\\0_data\\sce_Mast_cells_all.rds")
Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(sce, ndims = 50)
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #44
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #10
pcs <- min(co1, co2)
pcs  #10
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

library(ggplot2)
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()
pc_num= 1:pcs
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1
Idents(sce) <- "RNA_snn_res.0.1"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.1")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Mast_cells\\0_data\\sce_Mast_cells_fliter_afterdim.rds")


sce <- readRDS("D:\\Urinary system tumors\\work\\4_Mast_cells\\0_data\\sce_Mast_cells_fliter_afterdim.rds")
Idents(sce) <- "RNA_snn_res.0.1"
sce.markers <- FindAllMarkers(object = sce,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0.6) 
paste(sce.markers[which(sce.markers$cluster == 0 & sce.markers$pct.1 > 0.6),]$gene,collapse = ",")

sce@meta.data$Mast_cells_type <- recode(sce@meta.data$RNA_snn_res.0.1,
                                        "0"="Mast c1",
                                        "1"="Mast c2")

DimPlot(sce, label = F,reduction = "umap",group.by = "RNA_snn_res.0.1",repel = T)

DimPlot(sce, label = F,reduction = "umap",group.by = "Mast_cells_type",repel = T)+ 
  theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Mast_cells\\2_Mast_cells_anno\\Mast_cells_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce,features = Mast_cells_marker_genes,group.by = "Mast_cells_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_Mast_cells\\2_Mast_cells_anno\\Mast_cells_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()
saveRDS(sce,"D:\\Urinary system tumors\\work\\4_Mast_cells\\2_Mast_cells_anno\\Mast_cells_sce_anno.rds")

##B cell -anno ----
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")

sce_ccRCC_T <-  subset(sce_ccRCC_tumor, cell_type %in% c("B_cells"))
sce_ccRCC_T@meta.data$tumor_type <- "ccRCC"

sce_BC_T <-  subset(sce_BC_tumor, cell_type %in% c("B_cells"))
sce_BC_T@meta.data$tumor_type <- "BC"

sce_PCa_T <- subset(sce_PCa_tumor, cell_type %in% c("B_cells"))
sce_PCa_T@meta.data$tumor_type <- "PCa"


sce_T_all <- merge(sce_ccRCC_T,sce_BC_T)
sce_T_all <- merge(sce_T_all,sce_PCa_T)

table(sce_T_all@meta.data$tumor_type)
saveRDS(sce_T_all,"D:\\Urinary system tumors\\work\\4_B_cells\\0_data\\sce_B_cells_all.rds")
sce <- readRDS("D:\\Urinary system tumors\\work\\4_B_cells\\0_data\\sce_B_cells_all.rds")
sce2 <- readRDS("D:\\Urinary system tumors\\work\\4_T_cells\\2_T_cell_anno\\T_cell_sce_anno_Bcells.rds")
sce <- merge(sce,sce2)

Idents(sce) <- "tumor_type"
sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
#ScaleData()
sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)

DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(sce, ndims = 50)
pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1  #44
co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2  #4
pcs <- min(co1, co2)
pcs  #4
plot_df <-data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
ggplot(plot_df,aes(cumu, pct, label = rank,color = rank > pcs)) + 
  geom_text()+
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]),color = "grey")+
  theme_bw()
pc_num= 1:4
sce <- FindNeighbors(sce, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 
colnames(sce@meta.data)
library(clustree)
p1 <- clustree(sce, prefix = 'RNA_snn_res.')+coord_flip()
p1
Idents(sce) <- "RNA_snn_res.0.3"
sce <- sce %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)

DimPlot(sce, label = T, reduction = "tsne",cols = mycolors)
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5")
DimPlot(sce, label = T, reduction = "tsne",group.by = "tumor_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "tumor_type")

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_B_cells\\0_data\\sce_B_cells_fliter_afterdim.rds")
dim(sce)
sce <- readRDS("D:\\Urinary system tumors\\work\\4_B_cells\\0_data\\sce_B_cells_fliter_afterdim.rds")
B_cells_marker_genes <- c("MS4A1","NR4A1","TAGLN2",
                          "JCHAIN","MZB1","SDC1"
)
B_cells_marker_genes <- unique(B_cells_marker_genes)
VlnPlot(sce, features = B_cells_marker_genes,
        group.by = "RNA_snn_res.0.5",pt.size = 0)

exp_sce_T <- DotPlot(sce,features = B_cells_marker_genes,group.by = "RNA_snn_res.0.5")$data
min(exp_sce_T$avg.exp.scaled)
max(exp_sce_T$avg.exp.scaled)
#min:-0.6108159  max: 2.041241

DotPlot(sce,features = B_cells_marker_genes,group.by = "RNA_snn_res.0.1")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

sce@meta.data$B_cells_type <- recode(sce@meta.data$RNA_snn_res.0.1,
                                     "0"="PlasmaCells",
                                     "1"="Bmem",
                                     "2"="PlasmaCells",
                                     "3"="Bmem")

mycolors <- c("Bmem" = "#B5D4AC", 
              "PlasmaCells"= "#EBA07EFF")
DimPlot(sce, label = T, reduction = "umap",group.by = "B_cells_type")
DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.1")

DimPlot(sce, label = F,reduction = "umap",group.by = "B_cells_type",repel = T,cols = mycolors)+ 
  theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
ggsave(filename = "D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

DotPlot(sce,features = B_cells_marker_genes,group.by = "B_cells_type")+
  coord_flip()+theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
  scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-2.5,2.5)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
ggsave(filename = "D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_umap_cell_type_marker.pdf"
       ,width = 5, height = 5)
dev.off()

saveRDS(sce,"D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_sce_anno.rds")
B_cell <- readRDS("D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_sce_anno.rds")

data <- as.data.frame(cbind(B_cell@reductions[["umap"]]@cell.embeddings,as.character(B_cell@meta.data$B_cells_type),
                            as.character(B_cell@meta.data$B_cells_type)))
colnames(data) <- c("umap_x","umap_y","celltype","group") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)
mycolors <- c("Bmem" = "#B5D4AC", 
              "PlasmaCells"= "#EBA07EFF")
####Bmem
Bmem <- data[which(data$group == "Bmem"),c(1,2)]
colnames(Bmem) <- c("umap_x","umap_y")
Bmem$umap_x <- as.numeric(Bmem$umap_x)
Bmem$umap_y <- as.numeric(Bmem$umap_y)
Bmem <- kde(x=Bmem)
Bmem1 <- with(Bmem, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
Bmem2 <- with(Bmem, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[2]])
####PlasmaCells
PlasmaCells <- data[which(data$group == "PlasmaCells"),c(1,2)]
colnames(PlasmaCells) <- c("umap_x","umap_y")
PlasmaCells$umap_x <- as.numeric(PlasmaCells$umap_x)
PlasmaCells$umap_y <- as.numeric(PlasmaCells$umap_y)
PlasmaCells <- kde(x=PlasmaCells)
PlasmaCells <- with(PlasmaCells, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])

color <- c("#E3541F","#F5A623")
plotdata <- ggplot(data, aes(x = umap_x, 
                             y = umap_y,
                             color = celltype)) +
  geom_point(size = 0.001) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  guides(color = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_color_manual(values = color)
con.size=1.5
cl='black'

plotdata+geom_path(aes(x, y), data=data.frame(Bmem1),col="#9B3718",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Bmem2),col="#9B3718",linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(PlasmaCells),col="#9E691B",linetype = 2,size=con.size)

ggplot(data, aes(x = umap_x, 
                 y = umap_y,
                 color = celltype)) +
  geom_point(size = 2)+
  scale_color_manual(values = color)
ggsave(filename = "D:\\Urinary system tumors\\work\\4_B_cells\\2_B_cells_anno\\B_cells_umap_cell_type.pdf"
       ,width = 5, height = 5)
dev.off()

##Merge -anno ----

sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")

B_cell <- readRDS("D:/Urinary system tumors/work/4_B_cells/2_B_cells_anno/B_cells_sce_anno.rds")
Endothelium <- readRDS("D:/Urinary system tumors/work/4_Endothelium/2_Endothelium_anno/Endothelium_sce_anno.rds")
Myeloid <- readRDS("D:/Urinary system tumors/work/4_Myeloid/2_Myeloid_anno/Myeloid_sce_anno.rds")
T_cells <- readRDS("D:/Urinary system tumors/work/4_T_cells/2_T_cell_anno/T_cell_sce_anno.rds")
Fibroblast <- readRDS("D:/Urinary system tumors/work/4_Fibroblast/2_Fibroblast_anno/Fibroblast_sce_anno.rds")
Mast_cells <- readRDS("D:/Urinary system tumors/work/4_Mast_cells/2_Mast_cells_anno/Mast_cells_sce_anno.rds")
#ccRCC
ccRCC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor.rds")
ccRCC <-  subset(ccRCC, cell_type %in% c("Epithelium"))
ccRCC@meta.data$cell_type_second <- "Normal epithelium"
ccRCC@meta.data[which(ccRCC@meta.data$cell_type == "Epithelium"&ccRCC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"


#BC
BC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/2_BC/2_BC_after_ident_tumor.rds")
BC <-  subset(BC, cell_type %in% c("Epithelium"))
BC@meta.data$cell_type_second <- "Normal epithelium"
BC@meta.data[which(BC@meta.data$cell_type == "Epithelium"&BC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"


#PCa
PCa <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/3_PCa/2_PCa_after_ident_tumor.rds")
PCa <-  subset(PCa, cell_type %in% c("Epithelium"))
PCa@meta.data$cell_type_second <- "Normal epithelium"
PCa@meta.data[which(PCa@meta.data$cell_type == "Epithelium"&PCa@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"

##
ccRCC_B_cell <- B_cell@meta.data[which(B_cell@meta.data$tumor_type == "ccRCC"),]
ccRCC_Endothelium <- Endothelium@meta.data[which(Endothelium@meta.data$tumor_type == "ccRCC"),]
ccRCC_Myeloid <- Myeloid@meta.data[which(Myeloid@meta.data$tumor_type == "ccRCC"),]
ccRCC_T_cells <- T_cells@meta.data[which(T_cells@meta.data$tumor_type == "ccRCC"),]
ccRCC_Fibroblast <- Fibroblast@meta.data[which(Fibroblast@meta.data$tumor_type == "ccRCC"),]
ccRCC_Mast_cells <- Mast_cells@meta.data[which(Mast_cells@meta.data$tumor_type == "ccRCC"),]

#ccRCC second cell type
sce_ccRCC_tumor@meta.data$cell_type_second <- NA
sce_ccRCC_tumor@meta.data[row.names(ccRCC_B_cell),]$cell_type_second <- as.character(ccRCC_B_cell$B_cells_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC_Endothelium),]$cell_type_second <- as.character(ccRCC_Endothelium$Endothelium_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC_Myeloid),]$cell_type_second <- as.character(ccRCC_Myeloid$Myeloid_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC_T_cells),]$cell_type_second <- as.character(ccRCC_T_cells$T_cell_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC_Fibroblast),]$cell_type_second <- as.character(ccRCC_Fibroblast$Fibroblast_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC_Mast_cells),]$cell_type_second <- as.character(ccRCC_Mast_cells$Mast_cells_type)
sce_ccRCC_tumor@meta.data[row.names(ccRCC@meta.data),]$cell_type_second <- as.character(ccRCC@meta.data$cell_type_second)

sce_ccRCC_tumor@meta.data$filter <- as.character(!is.na(sce_ccRCC_tumor@meta.data$cell_type_second))

saveRDS(sce_ccRCC_tumor,"D:\\Urinary system tumors\\work\\4_merge\\ccRCC_sce_anno_second.rds")

BC_B_cell <- B_cell@meta.data[which(B_cell@meta.data$tumor_type == "BC"),]
BC_Endothelium <- Endothelium@meta.data[which(Endothelium@meta.data$tumor_type == "BC"),]
BC_Myeloid <- Myeloid@meta.data[which(Myeloid@meta.data$tumor_type == "BC"),]
BC_T_cells <- T_cells@meta.data[which(T_cells@meta.data$tumor_type == "BC"),]
BC_Fibroblast <- Fibroblast@meta.data[which(Fibroblast@meta.data$tumor_type == "BC"),]
BC_Mast_cells <- Mast_cells@meta.data[which(Mast_cells@meta.data$tumor_type == "BC"),]

#BC second cell type
sce_BC_tumor@meta.data$cell_type_second <- NA
sce_BC_tumor@meta.data[row.names(BC_B_cell),]$cell_type_second <- as.character(BC_B_cell$B_cells_type)
sce_BC_tumor@meta.data[row.names(BC_Endothelium),]$cell_type_second <- as.character(BC_Endothelium$Endothelium_type)
sce_BC_tumor@meta.data[row.names(BC_Myeloid),]$cell_type_second <- as.character(BC_Myeloid$Myeloid_type)
sce_BC_tumor@meta.data[row.names(BC_T_cells),]$cell_type_second <- as.character(BC_T_cells$T_cell_type)
sce_BC_tumor@meta.data[row.names(BC_Fibroblast),]$cell_type_second <- as.character(BC_Fibroblast$Fibroblast_type)
sce_BC_tumor@meta.data[row.names(BC_Mast_cells),]$cell_type_second <- as.character(BC_Mast_cells$Mast_cells_type)
sce_BC_tumor@meta.data[row.names(BC@meta.data),]$cell_type_second <- as.character(BC@meta.data$cell_type_second)

sce_BC_tumor@meta.data$filter <- as.character(!is.na(sce_BC_tumor@meta.data$cell_type_second))


saveRDS(sce_BC_tumor,"D:\\Urinary system tumors\\work\\4_merge\\BC_sce_anno_second.rds")


PCa_B_cell <- B_cell@meta.data[which(B_cell@meta.data$tumor_type == "PCa"),]
PCa_Endothelium <- Endothelium@meta.data[which(Endothelium@meta.data$tumor_type == "PCa"),]
PCa_Myeloid <- Myeloid@meta.data[which(Myeloid@meta.data$tumor_type == "PCa"),]
PCa_T_cells <- T_cells@meta.data[which(T_cells@meta.data$tumor_type == "PCa"),]
PCa_Fibroblast <- Fibroblast@meta.data[which(Fibroblast@meta.data$tumor_type == "PCa"),]
PCa_Mast_cells <- Mast_cells@meta.data[which(Mast_cells@meta.data$tumor_type == "PCa"),]

#PCa second cell type
sce_PCa_tumor@meta.data$cell_type_second <- NA
sce_PCa_tumor@meta.data[row.names(PCa_B_cell),]$cell_type_second <- as.character(PCa_B_cell$B_cells_type)
sce_PCa_tumor@meta.data[row.names(PCa_Endothelium),]$cell_type_second <- as.character(PCa_Endothelium$Endothelium_type)
sce_PCa_tumor@meta.data[row.names(PCa_Myeloid),]$cell_type_second <- as.character(PCa_Myeloid$Myeloid_type)
sce_PCa_tumor@meta.data[row.names(PCa_T_cells),]$cell_type_second <- as.character(PCa_T_cells$T_cell_type)
sce_PCa_tumor@meta.data[row.names(PCa_Fibroblast),]$cell_type_second <- as.character(PCa_Fibroblast$Fibroblast_type)
sce_PCa_tumor@meta.data[row.names(PCa_Mast_cells),]$cell_type_second <- as.character(PCa_Mast_cells$Mast_cells_type)
sce_PCa_tumor@meta.data[row.names(PCa@meta.data),]$cell_type_second <- as.character(PCa@meta.data$cell_type_second)
sce_PCa_tumor@meta.data$filter <- as.character(!is.na(sce_PCa_tumor@meta.data$cell_type_second))
saveRDS(sce_PCa_tumor,"D:\\Urinary system tumors\\work\\4_merge\\PCa_sce_anno_second.rds")

#Extract ecotypes code -----
#Running on Linux system
#Rscript EcoTyper_discovery_scRNA.R -c config_discovery_scRNA.yml

Ecotypes <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Ecotypes/ecotypes.txt",header = T)
###B_CELL
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/B_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/B_cells/gene_info.txt")

B_cells_S01 <- list()
B_cells_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
B_cells_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

B_cells_S04 <- list()
B_cells_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
B_cells_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

B_cells <- list()
B_cells[["B_cells_S01"]] <- B_cells_S01
B_cells[["B_cells_S04"]] <- B_cells_S04

###Mast_cells
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Mast_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Mast_cells/gene_info.txt")
Mast_cells_S01 <- list()
Mast_cells_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Mast_cells_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Mast_cells_S03 <- list()
Mast_cells_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
Mast_cells_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

Mast_cells_S05 <- list()
Mast_cells_S05[["cell"]] <- cell[which(cell$State == "S05"),]$ID
Mast_cells_S05[["gene"]] <- geneinfo[which(geneinfo$State == "S05"),]$Gene

Mast_cells <- list()
Mast_cells[["Mast_cells_S01"]] <- Mast_cells_S01
Mast_cells[["Mast_cells_S03"]] <- Mast_cells_S03
Mast_cells[["Mast_cells_S05"]] <- Mast_cells_S05


###T_cells
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/T_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/T_cells/gene_info.txt")

T_cells_S03 <- list()
T_cells_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
T_cells_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

T_cells_S04 <- list()
T_cells_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
T_cells_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

T_cells_S07 <- list()
T_cells_S07[["cell"]] <- cell[which(cell$State == "S07"),]$ID
T_cells_S07[["gene"]] <- geneinfo[which(geneinfo$State == "S07"),]$Gene

T_cells_S08 <- list()
T_cells_S08[["cell"]] <- cell[which(cell$State == "S08"),]$ID
T_cells_S08[["gene"]] <- geneinfo[which(geneinfo$State == "S08"),]$Gene

T_cells <- list()
T_cells[["T_cells_S03"]] <- T_cells_S03
T_cells[["T_cells_S04"]] <- T_cells_S04
T_cells[["T_cells_S07"]] <- T_cells_S07
T_cells[["T_cells_S08"]] <- T_cells_S08

###Endothelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Endothelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Endothelium/gene_info.txt")

Endothelium_S07 <- list()
Endothelium_S07[["cell"]] <- cell[which(cell$State == "S07"),]$ID
Endothelium_S07[["gene"]] <- geneinfo[which(geneinfo$State == "S07"),]$Gene

Endothelium_S04 <- list()
Endothelium_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Endothelium_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Endothelium <- list()
Endothelium[["Endothelium_S04"]] <- Endothelium_S04
Endothelium[["Endothelium_S07"]] <- Endothelium_S07


#Epithelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Epithelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Epithelium/gene_info.txt")

Epithelium_S03 <- list()
Epithelium_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
Epithelium_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

Epithelium_S04 <- list()
Epithelium_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Epithelium_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Epithelium <- list()
Epithelium[["Epithelium_S03"]] <- Epithelium_S03
Epithelium[["Epithelium_S04"]] <- Epithelium_S04


###Myeloid
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Myeloid/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Myeloid/gene_info.txt")

Myeloid_S03 <- list()
Myeloid_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
Myeloid_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

Myeloid_S04 <- list()
Myeloid_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Myeloid_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Myeloid_S01 <- list()
Myeloid_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Myeloid_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Myeloid_S02 <- list()
Myeloid_S02[["cell"]] <- cell[which(cell$State == "S02"),]$ID
Myeloid_S02[["gene"]] <- geneinfo[which(geneinfo$State == "S02"),]$Gene

Myeloid <- list()
Myeloid[["Myeloid_S01"]] <- Myeloid_S01
Myeloid[["Myeloid_S02"]] <- Myeloid_S02
Myeloid[["Myeloid_S03"]] <- Myeloid_S03
Myeloid[["Myeloid_S04"]] <- Myeloid_S04


###Fibroblast
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Fibroblast/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_ccRCC/Fibroblast/gene_info.txt")

Fibroblast_S03 <- list()
Fibroblast_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
Fibroblast_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

Fibroblast_S04 <- list()
Fibroblast_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Fibroblast_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Fibroblast_S02 <- list()
Fibroblast_S02[["cell"]] <- cell[which(cell$State == "S02"),]$ID
Fibroblast_S02[["gene"]] <- geneinfo[which(geneinfo$State == "S02"),]$Gene

Fibroblast <- list()
Fibroblast[["Fibroblast_S02"]] <- Fibroblast_S02
Fibroblast[["Fibroblast_S03"]] <- Fibroblast_S03
Fibroblast[["Fibroblast_S04"]] <- Fibroblast_S04


Ecotypes[which(Ecotypes$Ecotype == "E1"),]$ID
#"B_cells_S01"    "Mast_cells_S05" "T_cells_S03"
E1 <- list()
E1[["B_cells_S01"]] <- B_cells[["B_cells_S01"]]
E1[["Mast_cells_S05"]] <- Mast_cells[["Mast_cells_S05"]]
E1[["T_cells_S03"]] <- T_cells[["T_cells_S03"]]

Ecotypes[which(Ecotypes$Ecotype == "E2"),]$ID
#"B_cells_S04"     "Endothelium_S07" "Epithelium_S03"  "Myeloid_S02" 
E2 <- list()
E2[["B_cells_S04"]] <- B_cells[["B_cells_S04"]]
E2[["Endothelium_S07"]] <- Endothelium[["Endothelium_S07"]]
E2[["Epithelium_S03"]] <- Epithelium[["Epithelium_S03"]]
E2[["Myeloid_S02"]] <- Myeloid[["Myeloid_S02"]]

Ecotypes[which(Ecotypes$Ecotype == "E3"),]$ID
#"Endothelium_S04"     "Fibroblast_S03" "Mast_cells_S03"  "T_cells_S04" 
E3 <- list()
E3[["Endothelium_S04"]] <- Endothelium[["Endothelium_S04"]]
E3[["Fibroblast_S03"]] <- Fibroblast[["Fibroblast_S03"]]
E3[["Mast_cells_S03"]] <- Mast_cells[["Mast_cells_S03"]]
E3[["T_cells_S04"]] <- T_cells[["T_cells_S04"]]

Ecotypes[which(Ecotypes$Ecotype == "E4"),]$ID
#"Epithelium_S04" "Myeloid_S01"    "T_cells_S07"
E4 <- list()
E4[["Epithelium_S04"]] <- Epithelium[["Epithelium_S04"]]
E4[["Myeloid_S01"]] <- Myeloid[["Myeloid_S01"]]
E4[["T_cells_S07"]] <- T_cells[["T_cells_S07"]]

Ecotypes[which(Ecotypes$Ecotype == "E5"),]$ID
# "Fibroblast_S02" "Mast_cells_S01" "Myeloid_S04" 
E5 <- list()
E5[["Fibroblast_S02"]] <- Fibroblast[["Fibroblast_S02"]]
E5[["Mast_cells_S01"]] <- Mast_cells[["Mast_cells_S01"]]
E5[["Myeloid_S04"]] <- Myeloid[["Myeloid_S04"]]

Ecotypes[which(Ecotypes$Ecotype == "E6"),]$ID
# "Fibroblast_S04" "Myeloid_S03"    "T_cells_S08" 
E6 <- list()
E6[["Fibroblast_S04"]] <- Fibroblast[["Fibroblast_S04"]]
E6[["Myeloid_S03"]] <- Myeloid[["Myeloid_S03"]]
E6[["T_cells_S08"]] <- T_cells[["T_cells_S08"]]

ccRCC_Ecotypes <- list()
ccRCC_Ecotypes[["E1"]] <- E1
ccRCC_Ecotypes[["E2"]] <- E2
ccRCC_Ecotypes[["E3"]] <- E3
ccRCC_Ecotypes[["E4"]] <- E4
ccRCC_Ecotypes[["E5"]] <- E5
ccRCC_Ecotypes[["E6"]] <- E6
rm(list=setdiff(ls(), "ccRCC_Ecotypes"))


Ecotypes <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Ecotypes/ecotypes.txt",header = T)

###B_CELL
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/B_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/B_cells/gene_info.txt")

B_cells_S01 <- list()
B_cells_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
B_cells_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

B_cells <- list()
B_cells[["B_cells_S01"]] <- B_cells_S01


###T_cells
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/T_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/T_cells/gene_info.txt")

T_cells_S01 <- list()
T_cells_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
T_cells_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

T_cells_S02 <- list()
T_cells_S02[["cell"]] <- cell[which(cell$State == "S02"),]$ID
T_cells_S02[["gene"]] <- geneinfo[which(geneinfo$State == "S02"),]$Gene

T_cells_S05 <- list()
T_cells_S05[["cell"]] <- cell[which(cell$State == "S05"),]$ID
T_cells_S05[["gene"]] <- geneinfo[which(geneinfo$State == "S05"),]$Gene


T_cells <- list()
T_cells[["T_cells_S01"]] <- T_cells_S01
T_cells[["T_cells_S02"]] <- T_cells_S02
T_cells[["T_cells_S05"]] <- T_cells_S05


###Endothelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Endothelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Endothelium/gene_info.txt")

Endothelium_S01 <- list()
Endothelium_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Endothelium_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Endothelium_S04 <- list()
Endothelium_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Endothelium_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Endothelium_S05 <- list()
Endothelium_S05[["cell"]] <- cell[which(cell$State == "S05"),]$ID
Endothelium_S05[["gene"]] <- geneinfo[which(geneinfo$State == "S05"),]$Gene

Endothelium_S06 <- list()
Endothelium_S06[["cell"]] <- cell[which(cell$State == "S06"),]$ID
Endothelium_S06[["gene"]] <- geneinfo[which(geneinfo$State == "S06"),]$Gene

Endothelium <- list()
Endothelium[["Endothelium_S01"]] <- Endothelium_S01
Endothelium[["Endothelium_S04"]] <- Endothelium_S04
Endothelium[["Endothelium_S05"]] <- Endothelium_S05
Endothelium[["Endothelium_S06"]] <- Endothelium_S06



###Epithelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Epithelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Epithelium/gene_info.txt")

Epithelium_S07 <- list()
Epithelium_S07[["cell"]] <- cell[which(cell$State == "S07"),]$ID
Epithelium_S07[["gene"]] <- geneinfo[which(geneinfo$State == "S07"),]$Gene


Epithelium <- list()
Epithelium[["Epithelium_S07"]] <- Epithelium_S07


###Myeloid
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Myeloid/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Myeloid/gene_info.txt")

Myeloid_S01 <- list()
Myeloid_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Myeloid_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Myeloid_S06 <- list()
Myeloid_S06[["cell"]] <- cell[which(cell$State == "S06"),]$ID
Myeloid_S06[["gene"]] <- geneinfo[which(geneinfo$State == "S06"),]$Gene

Myeloid <- list()
Myeloid[["Myeloid_S01"]] <- Myeloid_S01
Myeloid[["Myeloid_S06"]] <- Myeloid_S06


###Fibroblast
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Fibroblast/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_BC/Fibroblast/gene_info.txt")

Fibroblast_S09 <- list()
Fibroblast_S09[["cell"]] <- cell[which(cell$State == "S09"),]$ID
Fibroblast_S09[["gene"]] <- geneinfo[which(geneinfo$State == "S09"),]$Gene

Fibroblast_S05 <- list()
Fibroblast_S05[["cell"]] <- cell[which(cell$State == "S05"),]$ID
Fibroblast_S05[["gene"]] <- geneinfo[which(geneinfo$State == "S05"),]$Gene

Fibroblast_S02 <- list()
Fibroblast_S02[["cell"]] <- cell[which(cell$State == "S02"),]$ID
Fibroblast_S02[["gene"]] <- geneinfo[which(geneinfo$State == "S02"),]$Gene

Fibroblast <- list()
Fibroblast[["Fibroblast_S02"]] <- Fibroblast_S02
Fibroblast[["Fibroblast_S05"]] <- Fibroblast_S05
Fibroblast[["Fibroblast_S09"]] <- Fibroblast_S09


Ecotypes[which(Ecotypes$Ecotype == "E1"),]$ID
#"B_cells_S01"     "Endothelium_S05" "Myeloid_S06"     "T_cells_S05"
E1 <- list()
E1[["B_cells_S01"]] <- B_cells[["B_cells_S01"]]
E1[["Endothelium_S05"]] <- Endothelium[["Endothelium_S05"]]
E1[["Myeloid_S06"]] <- Myeloid[["Myeloid_S06"]]
E1[["T_cells_S05"]] <- T_cells[["T_cells_S05"]]

Ecotypes[which(Ecotypes$Ecotype == "E2"),]$ID
#"Endothelium_S01"     "Epithelium_S07" "Fibroblast_S09"
E2 <- list()
E2[["Endothelium_S01"]] <- Endothelium[["Endothelium_S01"]]
E2[["Epithelium_S07"]] <- Epithelium[["Epithelium_S07"]]
E2[["Fibroblast_S09"]] <- Fibroblast[["Fibroblast_S09"]]

Ecotypes[which(Ecotypes$Ecotype == "E3"),]$ID
#"Endothelium_S04" "Fibroblast_S05"  "T_cells_S01" 
E3 <- list()
E3[["Endothelium_S04"]] <- Endothelium[["Endothelium_S04"]]
E3[["Fibroblast_S05"]] <- Fibroblast[["Fibroblast_S05"]]
E3[["T_cells_S01"]] <- T_cells[["T_cells_S01"]]

Ecotypes[which(Ecotypes$Ecotype == "E4"),]$ID
# "Endothelium_S06" "Fibroblast_S02"  "Myeloid_S01"     "T_cells_S02"
E4 <- list()
E4[["Endothelium_S06"]] <- Endothelium[["Endothelium_S06"]]
E4[["Fibroblast_S02"]] <- Fibroblast[["Fibroblast_S02"]]
E4[["Myeloid_S01"]] <- Myeloid[["Myeloid_S01"]]
E4[["T_cells_S02"]] <- T_cells[["T_cells_S02"]]


BC_Ecotypes <- list()
BC_Ecotypes[["E1"]] <- E1
BC_Ecotypes[["E2"]] <- E2
BC_Ecotypes[["E3"]] <- E3
BC_Ecotypes[["E4"]] <- E4
rm(list=setdiff(ls(), "BC_Ecotypes"))

###
#########PCa
Ecotypes <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Ecotypes/ecotypes.txt",header = T)

###B_CELL
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/B_cells/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/B_cells/gene_info.txt")

B_cells_S01 <- list()
B_cells_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
B_cells_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene


B_cells <- list()
B_cells[["B_cells_S01"]] <- B_cells_S01

###Endothelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Endothelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Endothelium/gene_info.txt")

Endothelium_S01 <- list()
Endothelium_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Endothelium_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Endothelium_S04 <- list()
Endothelium_S04[["cell"]] <- cell[which(cell$State == "S04"),]$ID
Endothelium_S04[["gene"]] <- geneinfo[which(geneinfo$State == "S04"),]$Gene

Endothelium <- list()
Endothelium[["Endothelium_S01"]] <- Endothelium_S01
Endothelium[["Endothelium_S04"]] <- Endothelium_S04

###Fibroblast
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Fibroblast/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Fibroblast/gene_info.txt")

Fibroblast_S03 <- list()
Fibroblast_S03[["cell"]] <- cell[which(cell$State == "S03"),]$ID
Fibroblast_S03[["gene"]] <- geneinfo[which(geneinfo$State == "S03"),]$Gene

Fibroblast_S07 <- list()
Fibroblast_S07[["cell"]] <- cell[which(cell$State == "S07"),]$ID
Fibroblast_S07[["gene"]] <- geneinfo[which(geneinfo$State == "S07"),]$Gene

Fibroblast <- list()
Fibroblast[["Fibroblast_S07"]] <- Fibroblast_S07
Fibroblast[["Fibroblast_S03"]] <- Fibroblast_S03

###Myeloid
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Myeloid/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Myeloid/gene_info.txt")
Myeloid_S05 <- list()
Myeloid_S05[["cell"]] <- cell[which(cell$State == "S05"),]$ID
Myeloid_S05[["gene"]] <- geneinfo[which(geneinfo$State == "S05"),]$Gene

Myeloid <- list()
Myeloid[["Myeloid_S05"]] <- Myeloid_S05

###Epithelium
cell <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Epithelium/state_assignment.txt")
geneinfo <- read.table("D:/Urinary system tumors/work/5_Ecotype/output_PCa/Epithelium/gene_info.txt")
Epithelium_S01 <- list()
Epithelium_S01[["cell"]] <- cell[which(cell$State == "S01"),]$ID
Epithelium_S01[["gene"]] <- geneinfo[which(geneinfo$State == "S01"),]$Gene

Epithelium_S02 <- list()
Epithelium_S02[["cell"]] <- cell[which(cell$State == "S02"),]$ID
Epithelium_S02[["gene"]] <- geneinfo[which(geneinfo$State == "S02"),]$Gene

Epithelium_S07 <- list()
Epithelium_S07[["cell"]] <- cell[which(cell$State == "S07"),]$ID
Epithelium_S07[["gene"]] <- geneinfo[which(geneinfo$State == "S07"),]$Gene

Epithelium <- list()
Epithelium[["Epithelium_S01"]] <- Epithelium_S01
Epithelium[["Epithelium_S02"]] <- Epithelium_S02
Epithelium[["Epithelium_S07"]] <- Epithelium_S07

Ecotypes[which(Ecotypes$Ecotype == "E1"),]$ID
E1 <- list()
E1[["B_cells_S01"]] <- B_cells[["B_cells_S01"]]
E1[["Epithelium_S07"]] <- Epithelium[["Epithelium_S07"]]
E1[["Myeloid_S05"]] <- Myeloid[["Myeloid_S05"]]

Ecotypes[which(Ecotypes$Ecotype == "E2"),]$ID
E2 <- list()
E2[["Endothelium_S01"]] <- Endothelium[["Endothelium_S01"]]
E2[["Epithelium_S01"]] <- Epithelium[["Epithelium_S01"]]
E2[["Fibroblast_S07"]] <- Fibroblast[["Fibroblast_S07"]]

Ecotypes[which(Ecotypes$Ecotype == "E3"),]$ID
E3 <- list()
E3[["Endothelium_S04"]] <- Endothelium[["Endothelium_S04"]]
E3[["Epithelium_S02"]] <- Epithelium[["Epithelium_S02"]]
E3[["Fibroblast_S03"]] <- Fibroblast[["Fibroblast_S03"]]

PCa_Ecotypes <- list()
PCa_Ecotypes[["E1"]] <- E1
PCa_Ecotypes[["E2"]] <- E2
PCa_Ecotypes[["E3"]] <- E3
rm(list=setdiff(ls(), "PCa_Ecotypes"))

#CellphoneDB ----
#Each ecological type takes turns running
cpdb_file_path = './cellphonedb.zip'
meta_file_path = './ccRCC/ccRCC_E1/metadata.tsv' 
counts_file_path ='./ccRCC/ccRCC_E1/CellphoneDB_count.txt'
out_path = './ccRCC/ccRCC_E1/result_cpdb'
Running on Linux system
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
  meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
  counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
  #microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
  score_interactions = True,                       # optional: whether to score interactions or not.
  iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
  threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
  threads = 30,                                     # number of threads to use in the analysis.
  debug_seed = 42,                                 # debug randome seed. To disable >=0.
  result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
  pvalue = 0.05,                                   # P-value threshold to employ for significance.
  subsampling = False,                             # To enable subsampling the data (geometri sketching).
  subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
  subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
  subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
  separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
  output_path = out_path,                          # Path to save results.
  output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
)
#query_minimum_score = 50
from cellphonedb.utils import search_utils
search_results = search_utils.search_analysis_results(
  query_cell_types_1 = ["Bmem","CD4Th17","PlasmaCells","NK",
                        "CTL"],  # List of cells 1, will be paired to cells 2 (list or 'All').
  query_cell_types_2 = ["Bmem","CD4Th17","PlasmaCells","NK",
                        "CTL"],     # List of cells 2, will be paired to cells 1 (list or 'All').
  significant_means = cpdb_results['significant_means'],          # significant_means file generated by CellphoneDB.
  deconvoluted = cpdb_results['deconvoluted'],                    # devonvoluted file generated by CellphoneDB.
  interaction_scores = cpdb_results['interaction_scores'],        # interaction score generated by CellphoneDB.
  query_minimum_score = 50,                                       # minimum score that an interaction must have to be filtered.
  separator = '|',                                                # separator (default: |) employed to split cells (cellA|cellB).
  long_format = True                                             # converts the output into a wide table, removing non-significant interactions
)
search_results.head()
search_results.to_csv("./ccRCC/ccRCC_E1/result_cpdb/search_results.txt", sep="\t")

##cellphonedb network -----
setwd("D:/Urinary system tumors/work/5_Ecotype/cellphonedb")
shuchu <- list()
shuchu[["ccRCC_E1"]] <- ccRCC_E1
shuchu[["ccRCC_E2"]] <- ccRCC_E2
shuchu[["ccRCC_E3"]] <- ccRCC_E3
shuchu[["ccRCC_E4"]] <- ccRCC_E4
shuchu[["ccRCC_E5"]] <- ccRCC_E5
shuchu[["ccRCC_E6"]] <- ccRCC_E6
for(i in names(shuchu)){
  count_raw <- shuchu[[i]]@assays$RNA@counts
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm, paste0("ccRCC/",i, "/CellphoneDB_count.txt"), sep='\t', quote=F)
  metadata <- shuchu[[i]]@meta.data[,c("orig.ident","cell_type_second")]
  metadata$barcode <- row.names(metadata)
  metadata <- metadata %>% dplyr::select(3,2)
  write_tsv(metadata, paste0("ccRCC/",i, "/metadata.tsv"))
}

shuchu <- list()
shuchu[["BC_E1"]] <- BC_E1
shuchu[["BC_E2"]] <- BC_E2
shuchu[["BC_E3"]] <- BC_E3
shuchu[["BC_E4"]] <- BC_E4
for(i in names(shuchu)){
  count_raw <- shuchu[[i]]@assays$RNA@counts
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm, paste0("BC/",i, "/CellphoneDB_count.txt"), sep='\t', quote=F)
  metadata <- shuchu[[i]]@meta.data[,c("orig.ident","cell_type_second")]
  metadata$barcode <- row.names(metadata)
  metadata <- metadata %>% dplyr::select(3,2)
  write_tsv(metadata, paste0("BC/",i, "/metadata.tsv"))
}



shuchu <- list()
shuchu[["PCa_E1"]] <- PCa_E1
shuchu[["PCa_E2"]] <- PCa_E2
shuchu[["PCa_E3"]] <- PCa_E3
for(i in names(shuchu)){
  count_raw <- shuchu[[i]]@assays$RNA@counts
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm, paste0("PCa/",i, "/CellphoneDB_count.txt"), sep='\t', quote=F)
  metadata <- shuchu[[i]]@meta.data[,c("orig.ident","cell_type_second")]
  metadata$barcode <- row.names(metadata)
  metadata <- metadata %>% dplyr::select(3,2)
  write_tsv(metadata, paste0("PCa/",i, "/metadata.tsv"))
}

library(igraph)
library(RColorBrewer)
rm(list = ls())
setwd("D:/Urinary system tumors/work/5_Ecotype/cellphonedb/PCa")
cancer_name <- "PCa_E3"
# Read the file of CellphoneDB
tmp <- list.files(paste0(cancer_name,"/result_cpdb/"))
pvals <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "pvalues",tmp)]),header = T,check.names = F)
colnames(pvals) <- gsub("nan","None",x=colnames(pvals))
means <- read.delim(paste0(cancer_name,"/result_cpdb/",tmp[grep(pattern = "analysis_means",tmp)]),header = T,check.names = F)
colnames(means) <- gsub("nan","None",x=colnames(means))
#interaction_scores <- read.delim(paste0(cancer_name,"/",tmp[grep(pattern = "interaction_scores",tmp)]),header = T,check.names = F)
library(dplyr)
#Ligand-Receptor
pvals <- pvals %>% dplyr::filter(directionality=="Ligand-Receptor")
means <- means %>% dplyr::filter(directionality=="Ligand-Receptor")
library(ktplots)
heatmap <- plot_cpdb_heatmap(pvals = pvals, cellheight = 10, cellwidth = 10,main = "Cell-Cell Interactions between Death States")
pdf(paste0(cancer_name,"/","heatmap.pdf"),width = 50,height = 50)
print(heatmap)
dev.off()

pvals2 <- pvals[pvals[,c(14:dim(pvals)[2])] < 0.05]

inter_cell <- c()
for(i in 14:dim(pvals)[2]){
  if(length(which(pvals[,i] < 0.05)) != 0){
    temp <- means[which(pvals[,i] < 0.05), c(5,6,i)]
    colnames(temp)[3] <- "cell_cell"
    temp$cell <- colnames(pvals)[i]
    inter_cell <- rbind(inter_cell,temp)
  }
}
inter_cell <- inter_cell[which(inter_cell$gene_a != ""),]
inter_cell <- inter_cell[which(inter_cell$gene_b != ""),]

inter_cell$indexa=str_replace(inter_cell$cell,"\\|.*$","")
inter_cell$indexb=str_replace(inter_cell$cell,"^.*\\|","")

write.table(inter_cell,paste0(cancer_name,"/","cell_cellnetwork.txt"),row.names = F,quote = F,sep = "\t")


statdf <- as.data.frame(table(inter_cell$indexa,inter_cell$indexb))
colnames(statdf) <- c("indexa","indexb","number")


rankname=sort(unique(c(inter_cell$indexa,inter_cell$indexb)) )

A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count1 <- statdf[statdf$indexa == i & statdf$indexb == j,"number"]
    count2 <- statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    count=count1+count2
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}

statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] 
write.table(statdf2,paste0(cancer_name,"/","network.txt"),row.names = F,quote = F,sep = "\t")
color1=c("ILC3" = "#F3CE33",
         "CTL" = "#E1835E",
         "CD8Trm"="#BA2371",
         "CD4Treg"="#C4D9C5",
         "CD4Th17"= "#48BBB9",
         "NK"="#4B59A6",
         "unkown"="#171831",
         "Monocytes" = "#921F58FF", 
         "cDC2"= "#CE599DFF", 
         "Macrophages"="#F1A890FF", 
         "cDC1"  = "#9DD5DBFF", 
         "Macrophages_CD83" = "#D5DED0FF", 
         "pDCs"="#68C2BFFF",
         "Myofibroblast" = "#A05B08", 
         "Lipofibroblast"= "#717C66",
         "LECs" = "#E38623",
         "VECs" = "#51BDB5",
         "CapECs"="#EDE7BB",
         "AECs"="#BC1A29",
         "ProliferatingECs_c1"="#443972",
         "ProliferatingECs_c2"="#377DB8",
         "Bmem" = "#B5D4AC", 
         "PlasmaCells"= "#EBA07EFF",
         "Mast c1"= "#F6E534",
         "Mast c2"= "#C334F7",
         "Tumor cells"= "#6D3030",
         "Normal epithelium" = "#F6CDCC",
         "unkwon"="#A3A3A3")
color2=colorRampPalette(brewer.pal(9, "Greys")[3:7])(max(statdf2$number))
names(color2)=1:max(statdf2$number) 

net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 20 
E(net)$color <- color2[as.character(ifelse(E(net)$number > 20,20,E(net)$number))]
E(net)$label = E(net)$number 
E(net)$label.color <- "black" 
V(net)$label.color <- "black" 
V(net)$color <- color1[names(V(net))] 

loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
pdf(paste0(cancer_name,"/","cell_interaction.pdf"),width = 5,height = 5)
plot(net,
     edge.arrow.size = 0, 
     edge.curved = 0, 
     vertex.frame.color = "black", 
     layout = coords,
     vertex.label.cex = 0.8, 
     vertex.size = 30) 
dev.off()

#Ecotype of spatial transcriptomics-----
##Integrate multiple samples -----
PCa_Spatial <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/PCa_Spatial.rds")
PCa <- merge(PCa_Spatial[[1]],PCa_Spatial[[2]])
PCa <- merge(PCa,PCa_Spatial[[3]])
PCa <- merge(PCa,PCa_Spatial[[4]])

DefaultAssay(PCa) <- "SCT"

DefaultAssay(PCa_Spatial[[1]]) <- "SCT"
PCa_Spatial[[1]] <- FindVariableFeatures(PCa_Spatial[[1]])
DefaultAssay(PCa_Spatial[[2]]) <- "SCT"
PCa_Spatial[[2]] <- FindVariableFeatures(PCa_Spatial[[2]])
DefaultAssay(PCa_Spatial[[3]]) <- "SCT"
PCa_Spatial[[3]] <- FindVariableFeatures(PCa_Spatial[[3]])
DefaultAssay(PCa_Spatial[[4]]) <- "SCT"
PCa_Spatial[[4]] <- FindVariableFeatures(PCa_Spatial[[4]])


VariableFeatures(PCa) <-c(VariableFeatures(PCa_Spatial[[1]]), VariableFeatures(PCa_Spatial[[2]]),
                          VariableFeatures(PCa_Spatial[[3]]), VariableFeatures(PCa_Spatial[[4]]))
ElbowPlot(PCa)
PCa <- RunPCA(PCa, verbose = FALSE)
PCa <- FindNeighbors(PCa, dims = 1:15)
PCa <- FindClusters(PCa, verbose = FALSE,resolution = 0.3)
PCa <- RunUMAP(PCa, dims = 1:15)
DimPlot(PCa, reduction = "umap", group.by = c("ident", "orig.ident"))


PCa@meta.data$cluster <- as.numeric(PCa@meta.data$seurat_clusters)
PCa@meta.data$cluster <- paste0("CC",PCa@meta.data$cluster)
PCa@meta.data$group <- "HG"
PCa@meta.data[which(PCa@meta.data$orig.ident == "Tumor01"),]$group <-"LG"
PCa@meta.data[which(PCa@meta.data$orig.ident == "Tumor02"),]$group <-"LG"

df <- as.data.frame(table(PCa@meta.data[,c("group" ,"orig.ident", "cluster")]))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

library(ggplot2)
library(ggalluvial)
ggplot(df,aes(y= Freq, axis1 = group,axis2 = cluster))+
  geom_alluvium(aes(fill = cluster))+
  geom_stratum(width=1/6,fill ="black",color ="grey")+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme

ggplot(data = df, mapping = aes(x = orig.ident, y = Freq, fill = cluster,))+ 
  geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme



Idents(PCa) <- PCa@meta.data$cluster
av <- AverageExpression(PCa , 
                        assays = "SCT")
av=av[[1]] 
cg=names(tail(sort(apply(av, 1, sd)),1000)) 
pheatmap::pheatmap(cor(av[cg,])) 
library(COSG)


input_sce = PCa
table(Idents(input_sce))
pro = 'cosg_seurat_clusters'
if(T){
  
  library(COSG)
  marker_cosg <- cosg(
    input_sce,
    groups='all',
    assay='SCT',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  
  sce.Scale <- ScaleData(input_sce ,features =  top_10  )  
  
  DoHeatmap(sce.Scale, top_10)
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
}

symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
#source('~/Code/com_go_kegg_ReactomePA_human.R')
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
com_go_kegg_ReactomePA_human <- function(symbols_list ,pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  
  gcSample = lapply(symbols_list, function(y){ 
    y=as.character(na.omit(select(org.Hs.eg.db,
                                  keys = y,
                                  columns = 'ENTREZID',
                                  keytype = 'SYMBOL')[,2])
    )
    y
  })
  gcSample
  
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_kegg.pdf'),width = 10,height = 8)
  
  
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_ReactomePA.pdf'),width = 10,height = 8)
  
  
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_BP_cluster_simplified.csv'))
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_CC_cluster_simplified.csv'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_MF_cluster_simplified.csv'))
  
}
com_go_kegg_ReactomePA_human(symbols_list, pro='PCa' )

library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)



expr <- AverageExpression(PCa, assays = "SCT", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)

library(GSVA)
gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

genesets <- msigdbr(species = "Homo sapiens", category = "C2") 

genesets <- subset(genesets, gs_subcat=="CP:KEGG",
                   select = c("gs_name", "gene_symbol")) %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "KEGG_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "KEGG_gsva_res.csv", row.names = F)


genesets <- msigdbr(species = "Homo sapiens", category = "C5")
genesets <- subset(genesets,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()


genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "GO_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "GO_gsva_res.csv", row.names = F)

select <- c()
for(i in 1:dim(gsva.df)[1]){
  select <- c(select, gsva.df$CC4[i] == max(gsva.df[i,2:10]))
}
gsva.df2 <- gsva.df[select,]
gsva.df2 <- gsva.df2[order(gsva.df2$CC4,decreasing = T),]
GO <- c("GOBP_REGULATION_OF_CARDIAC_VASCULAR_SMOOTH_MUSCLE_CELL_DIFFERENTIATION",
        "GOBP_CHOLINE_METABOLIC_PROCESS" ,
        "GOBP_POSITIVE_REGULATION_OF_HYDROGEN_PEROXIDE_METABOLIC_PROCESS",
        "GOBP_LIPOXIN_METABOLIC_PROCESS",
        "GOBP_POSITIVE_REGULATION_OF_CORTICOSTEROID_HORMONE_SECRETION")

df <- gsva.df2[GO,2:11]
library(pheatmap)
pdf("GO.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()


df <- read.delim("clipboard",header = T,row.names = 1) #read gsva result file
df <- df[order(df$CC4,decreasing = T),]
KEGG <- c("KEGG_LINOLEIC_ACID_METABOLISM",
          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")
df <- df[KEGG,]

pdf("KEGG.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()



df <- read.delim("clipboard",header = T,row.names = 1) #read gsva result file
df <- df[order(df$CC4,decreasing = T),]
hallmarker <- c("HALLMARK_ANGIOGENESIS",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_HYPOXIA",
                "HALLMARK_GLYCOLYSIS")
df <- df[hallmarker,]

pdf("hallmarker.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()
saveRDS(PCa,"PCa_spaital_merge.rds")
library(Seurat)
ccRCC_Spatial <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/ccRCC_Spatial.rds")
ccRCC <- merge(ccRCC_Spatial[[1]],ccRCC_Spatial[[2]])
ccRCC <- merge(ccRCC,ccRCC_Spatial[[3]])
ccRCC <- merge(ccRCC,ccRCC_Spatial[[4]])
ccRCC <- merge(ccRCC,ccRCC_Spatial[[5]])

DefaultAssay(ccRCC) <- "SCT"

DefaultAssay(ccRCC_Spatial[[1]]) <- "SCT"
ccRCC_Spatial[[1]] <- FindVariableFeatures(ccRCC_Spatial[[1]])
DefaultAssay(ccRCC_Spatial[[2]]) <- "SCT"
ccRCC_Spatial[[2]] <- FindVariableFeatures(ccRCC_Spatial[[2]])
DefaultAssay(ccRCC_Spatial[[3]]) <- "SCT"
ccRCC_Spatial[[3]] <- FindVariableFeatures(ccRCC_Spatial[[3]])
DefaultAssay(ccRCC_Spatial[[4]]) <- "SCT"
ccRCC_Spatial[[4]] <- FindVariableFeatures(ccRCC_Spatial[[4]])
DefaultAssay(ccRCC_Spatial[[5]]) <- "SCT"
ccRCC_Spatial[[5]] <- FindVariableFeatures(ccRCC_Spatial[[5]])


VariableFeatures(ccRCC) <-unique(c(VariableFeatures(ccRCC_Spatial[[1]]), VariableFeatures(ccRCC_Spatial[[2]]),
                                   VariableFeatures(ccRCC_Spatial[[3]]), VariableFeatures(ccRCC_Spatial[[4]]),
                                   VariableFeatures(ccRCC_Spatial[[5]])))
ElbowPlot(ccRCC)
ccRCC <- RunPCA(ccRCC, verbose = FALSE)
ccRCC <- FindNeighbors(ccRCC, dims = 1:30)
ccRCC <- FindClusters(ccRCC, verbose = FALSE,resolution = 0.5)
ccRCC <- RunUMAP(ccRCC, dims = 1:30)
DimPlot(ccRCC, reduction = "umap", group.by = c("ident", "orig.ident"))


ccRCC@meta.data$cluster <- as.numeric(ccRCC@meta.data$seurat_clusters)
ccRCC@meta.data$cluster <- paste0("CC",ccRCC@meta.data$cluster)
ccRCC@meta.data$group <- "HG"
ccRCC@meta.data[which(ccRCC@meta.data$orig.ident == "LG_3"),]$group <-"LG"

df <- as.data.frame(table(ccRCC@meta.data[,c("group" , "cluster")]))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

library(ggplot2)
library(ggalluvial)
ggplot(df,aes(y= Freq, axis1 = orig.ident,axis2 = cluster))+
  geom_alluvium(aes(fill = cluster))+
  geom_stratum(width=1/6,fill ="black",color ="grey")+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme

ggplot(data = df, mapping = aes(x = orig.ident, y = Freq, fill = cluster,))+ 
  geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme




df <- as.data.frame(table(ccRCC@meta.data[,c("group" , "cluster")]))
ggplot(df,aes(y= Freq, axis1 = group,axis2 = cluster))+
  geom_alluvium(aes(fill = cluster))+
  geom_stratum(width=1/6,fill ="black",color ="grey")+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme

E1 <- unique(c(ccRCC_Ecotypes[["E1"]][["B_cells_S01"]][["gene"]],
               ccRCC_Ecotypes[["E1"]][["Mast_cells_S05"]][["gene"]],
               ccRCC_Ecotypes[["E1"]][["T_cells_S03"]][["gene"]]))

E2 <- unique(c(ccRCC_Ecotypes[["E2"]][["B_cells_S04"]][["gene"]],
               ccRCC_Ecotypes[["E2"]][["Endothelium_S07"]][["gene"]],
               ccRCC_Ecotypes[["E2"]][["Epithelium_S03"]][["gene"]],
               ccRCC_Ecotypes[["E2"]][["Myeloid_S02"]][["gene"]]))

E3 <- unique(c(ccRCC_Ecotypes[["E3"]][["Endothelium_S04"]][["gene"]],
               ccRCC_Ecotypes[["E3"]][["Fibroblast_S03"]][["gene"]],
               ccRCC_Ecotypes[["E3"]][["Mast_cells_S03"]][["gene"]],
               ccRCC_Ecotypes[["E3"]][["T_cells_S04"]][["gene"]]))

E4 <- unique(c(ccRCC_Ecotypes[["E4"]][["Epithelium_S04"]][["gene"]],
               ccRCC_Ecotypes[["E4"]][["Myeloid_S01"]][["gene"]],
               ccRCC_Ecotypes[["E4"]][["T_cells_S07"]][["gene"]]))


E5 <- unique(c(ccRCC_Ecotypes[["E5"]][["Fibroblast_S02"]][["gene"]],
               ccRCC_Ecotypes[["E5"]][["Mast_cells_S01"]][["gene"]],
               ccRCC_Ecotypes[["E5"]][["Myeloid_S04"]][["gene"]]))

E6 <- unique(c(ccRCC_Ecotypes[["E6"]][["Fibroblast_S04"]][["gene"]],
               ccRCC_Ecotypes[["E6"]][["Myeloid_S03"]][["gene"]],
               ccRCC_Ecotypes[["E6"]][["T_cells_S08"]][["gene"]]))

genelist <- list("E1" = E1,
                 "E2" = E2,
                 "E3" = E3,
                 "E4" = E4,
                 "E5" = E5,
                 "E6" = E6)

input_sce <- ccRCC
input_sce <- AddModuleScore(input_sce,features = genelist,name = names(genelist))

df <- input_sce@meta.data[,c(9:16)]

ggplot(df, aes(x=cluster, y=E44)) +geom_boxplot()
saveRDS(ccRCC,"ccRCC_spaital_merge.rds")

##GSVA -----
Idents(PCa) <- PCa@meta.data$cluster
av <- AverageExpression(PCa , 
                        assays = "SCT")
av=av[[1]] 
cg=names(tail(sort(apply(av, 1, sd)),1000)) 
pheatmap::pheatmap(cor(av[cg,])) 
library(COSG)


input_sce = PCa
table(Idents(input_sce))
pro = 'cosg_seurat_clusters'
if(T){
  
  library(COSG)
  marker_cosg <- cosg(
    input_sce,
    groups='all',
    assay='SCT',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  
  sce.Scale <- ScaleData(input_sce ,features =  top_10  )  
  
  DoHeatmap(sce.Scale, top_10)
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
}

symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
#source('~/Code/com_go_kegg_ReactomePA_human.R')
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
com_go_kegg_ReactomePA_human <- function(symbols_list ,pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  gcSample = lapply(symbols_list, function(y){ 
    y=as.character(na.omit(select(org.Hs.eg.db,
                                  keys = y,
                                  columns = 'ENTREZID',
                                  keytype = 'SYMBOL')[,2])
    )
    y
  })
  gcSample
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_kegg.pdf'),width = 10,height = 8)
  
  
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_ReactomePA.pdf'),width = 10,height = 8)
  
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_BP_cluster_simplified.csv'))
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_CC_cluster_simplified.csv'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_MF_cluster_simplified.csv'))
  
}
com_go_kegg_ReactomePA_human(symbols_list, pro='PCa' )

library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)



expr <- AverageExpression(PCa, assays = "SCT", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)

library(GSVA)
gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

genesets <- msigdbr(species = "Homo sapiens", category = "C2") 

genesets <- subset(genesets, gs_subcat=="CP:KEGG",
                   select = c("gs_name", "gene_symbol")) %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "KEGG_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "KEGG_gsva_res.csv", row.names = F)

genesets <- msigdbr(species = "Homo sapiens", category = "C5")
genesets <- subset(genesets,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()


genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "GO_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "GO_gsva_res.csv", row.names = F)

select <- c()
for(i in 1:dim(gsva.df)[1]){
  select <- c(select, gsva.df$CC4[i] == max(gsva.df[i,2:10]))
}
gsva.df2 <- gsva.df[select,]
gsva.df2 <- gsva.df2[order(gsva.df2$CC4,decreasing = T),]
GO <- c("GOBP_REGULATION_OF_CARDIAC_VASCULAR_SMOOTH_MUSCLE_CELL_DIFFERENTIATION",
        "GOBP_CHOLINE_METABOLIC_PROCESS" ,
        "GOBP_POSITIVE_REGULATION_OF_HYDROGEN_PEROXIDE_METABOLIC_PROCESS",
        "GOBP_LIPOXIN_METABOLIC_PROCESS",
        "GOBP_POSITIVE_REGULATION_OF_CORTICOSTEROID_HORMONE_SECRETION")

df <- gsva.df2[GO,2:11]


df <- read.delim("clipboard",header = T,row.names = 1) #read gsva result file
select <- c()
for(i in 1:dim(df)[1]){
  select <- c(select, df$CC5[i] == max(df[i,2:9]))
}
df <- df[select,]
df <- df[order(df$CC5,decreasing = T),]
row.names(df)[grepl("META",row.names(df))]
GO <- c("GOBP_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
        "GOBP_INSULIN_METABOLIC_PROCESS" ,
        "GOBP_CAMP_METABOLIC_PROCESS",
        "GOBP_GLYCEROL_3_PHOSPHATE_METABOLIC_PROCESS")

df <- df[GO,]


library(pheatmap)
pdf("GO.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()


df <- read.delim("clipboard",header = T,row.names = 1)#read gsva result file
df <- df[order(df$CC5,decreasing = T),]
KEGG <- c("KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",
          "KEGG_ERBB_SIGNALING_PATHWAY")
df <- df[KEGG,]

pdf("KEGG.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()



df <- read.delim("clipboard",header = T,row.names = 1)#read gsva result file
df <- df[order(df$CC5,decreasing = T),]
hallmarker <- c("HALLMARK_ANDROGEN_RESPONSE",
                "HALLMARK_HEDGEHOG_SIGNALING",
                "HALLMARK_CHOLESTEROL_HOMEOSTASIS")
df <- df[hallmarker,]

pdf("hallmarker.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()
saveRDS(PCa,"PCa_spaital_merge.rds")
Idents(ccRCC) <- ccRCC@meta.data$cluster
av <- AverageExpression(ccRCC , 
                        assays = "SCT")
av=av[[1]] 
cg=names(tail(sort(apply(av, 1, sd)),1000)) 
pheatmap::pheatmap(cor(av[cg,])) 
library(COSG)


input_sce = ccRCC
table(Idents(input_sce))
pro = 'cosg_seurat_clusters'
if(T){
  
  library(COSG)
  marker_cosg <- cosg(
    input_sce,
    groups='all',
    assay='SCT',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  
  sce.Scale <- ScaleData(input_sce ,features =  top_10  )  
  
  DoHeatmap(sce.Scale, top_10)
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
}

symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
#source('~/Code/com_go_kegg_ReactomePA_human.R')
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
com_go_kegg_ReactomePA_human <- function(symbols_list ,pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(stringr)
  
  gcSample = lapply(symbols_list, function(y){ 
    y=as.character(na.omit(select(org.Hs.eg.db,
                                  keys = y,
                                  columns = 'ENTREZID',
                                  keytype = 'SYMBOL')[,2])
    )
    y
  })
  gcSample
  
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_kegg.pdf'),width = 10,height = 8)
  
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                       pvalueCutoff=0.05)
  dotplot(xx)  + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_ReactomePA.pdf'),width = 10,height = 8)
  
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_BP_cluster_simplified.csv'))
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_CC_cluster_simplified.csv'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont     = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_MF_cluster_simplified.csv'))
  
}
com_go_kegg_ReactomePA_human(symbols_list, pro='ccRCC' )

library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)



expr <- AverageExpression(ccRCC, assays = "SCT", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] 
expr <- as.matrix(expr)


library(GSVA)
gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

genesets <- msigdbr(species = "Homo sapiens", category = "C2") 

genesets <- subset(genesets, gs_subcat=="CP:KEGG",
                   select = c("gs_name", "gene_symbol")) %>% as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "KEGG_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "KEGG_gsva_res.csv", row.names = F)

genesets <- msigdbr(species = "Homo sapiens", category = "C5")
genesets <- subset(genesets,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()


genesets <- split(genesets$gene_symbol, genesets$gs_name)

gsvapar <- gsvaParam(expr, genesets, maxDiff=TRUE)
gsva.res <- gsva(gsvapar) 
saveRDS(gsva.res, "GO_gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "GO_gsva_res.csv", row.names = F)

select <- c()
for(i in 1:dim(gsva.df)[1]){
  select <- c(select, gsva.df$CC4[i] == max(gsva.df[i,2:10]))
}
gsva.df2 <- gsva.df[select,]
gsva.df2 <- gsva.df2[order(gsva.df2$CC4,decreasing = T),]
GO <- c("GOBP_REGULATION_OF_CARDIAC_VASCULAR_SMOOTH_MUSCLE_CELL_DIFFERENTIATION",
        "GOBP_CHOLINE_METABOLIC_PROCESS" ,
        "GOBP_POSITIVE_REGULATION_OF_HYDROGEN_PEROXIDE_METABOLIC_PROCESS",
        "GOBP_LIPOXIN_METABOLIC_PROCESS",
        "GOBP_POSITIVE_REGULATION_OF_CORTICOSTEROID_HORMONE_SECRETION")

df <- gsva.df2[GO,2:11]
library(pheatmap)
pdf("GO.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()


df <- read.delim("clipboard",header = T,row.names = 1)#read gsva result file
df <- df[order(df$CC4,decreasing = T),]
KEGG <- c("KEGG_LINOLEIC_ACID_METABOLISM",
          "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")
df <- df[KEGG,]

pdf("KEGG.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()

df <- read.delim("clipboard",header = T,row.names = 1) #read gsva result file
df <- df[order(df$CC4,decreasing = T),]
hallmarker <- c("HALLMARK_ANGIOGENESIS",
                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                "HALLMARK_HYPOXIA",
                "HALLMARK_GLYCOLYSIS")
df <- df[hallmarker,]

pdf("hallmarker.pdf",height = 3,width = 6)
pheatmap(df,scale ="row",
         border_color ="white",
         color = colorRampPalette(c('#5671B6','white','#AC2928'))(100),
         show_rownames=F,
         legend_breaks=c(-2.5,-1.25,0,1.25,2.5),
         legend_labels=c("-2.5","-1.25","0","1.25","2.5"),
         treeheight_row=40,
         treeheight_col=40,
         fontsize=6,
         fontsize_col=8,
         angle_col=90,
         annotation_legend=F,
         cluster_rows = F,
         cutree_cols = 10,
         cluster_cols = F
)
dev.off()

##Ecotype score on st data -----
PCa_spaital_merge <- readRDS("D:/Urinary system tumors/work/5_Ecotype_spaital/PCa/PCa_spaital_merge.rds")
load("D:/Urinary system tumors/work/5_Ecotype/PCa_Ecotypes.RData")
E1 <- c(PCa_Ecotypes[["E1"]][["B_cells_S01"]][["gene"]],
        PCa_Ecotypes[["E1"]][["Epithelium_S07"]][["gene"]],
        PCa_Ecotypes[["E1"]][["Myeloid_S05"]][["gene"]])

E2 <- c(PCa_Ecotypes[["E2"]][["Endothelium_S01"]][["gene"]],
        PCa_Ecotypes[["E2"]][["Epithelium_S01"]][["gene"]],
        PCa_Ecotypes[["E2"]][["Fibroblast_S07"]][["gene"]])

E3 <- c(PCa_Ecotypes[["E3"]][["Endothelium_S04"]][["gene"]],
        PCa_Ecotypes[["E3"]][["Epithelium_S02"]][["gene"]],
        PCa_Ecotypes[["E3"]][["Fibroblast_S03"]][["gene"]])

E.list <- list("E1" = E1,
               "E2" = E2,
               "E3" = E3)

PCa_spaital_merge <- AddModuleScore(PCa_spaital_merge,features = E.list,name = names(E.list))

df <- PCa_spaital_merge@meta.data[,c("cluster","E11","E22","E33")]


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


ggplot(df, aes(x=cluster, y=E33,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/PCa/E3_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E22,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/PCa/E2_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E11,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/PCa/E1_spaital.pdf",height = 3,width = 4)


ccRCC_spaital_merge <- readRDS("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/ccRCC_spaital_merge.rds")
load("D:/Urinary system tumors/work/5_Ecotype/ccRCC_Ecotypes.RData")
E1 <- c(ccRCC_Ecotypes[["E1"]][["B_cells_S01"]][["gene"]],
        ccRCC_Ecotypes[["E1"]][["Mast_cells_S05"]][["gene"]],
        ccRCC_Ecotypes[["E1"]][["T_cells_S03"]][["gene"]])

E2 <- c(ccRCC_Ecotypes[["E2"]][["B_cells_S04"]][["gene"]],
        ccRCC_Ecotypes[["E2"]][["Endothelium_S07"]][["gene"]],
        ccRCC_Ecotypes[["E2"]][["Epithelium_S03"]][["gene"]],
        ccRCC_Ecotypes[["E2"]][["Myeloid_S02"]][["gene"]])

E3 <- c(ccRCC_Ecotypes[["E3"]][["Endothelium_S04"]][["gene"]],
        ccRCC_Ecotypes[["E3"]][["Fibroblast_S03"]][["gene"]],
        ccRCC_Ecotypes[["E3"]][["Mast_cells_S03"]][["gene"]],
        ccRCC_Ecotypes[["E3"]][["T_cells_S04"]][["gene"]])


E4 <- c(ccRCC_Ecotypes[["E4"]][["Epithelium_S04"]][["gene"]],
        ccRCC_Ecotypes[["E4"]][["Myeloid_S01"]][["gene"]],
        ccRCC_Ecotypes[["E4"]][["T_cells_S07"]][["gene"]])

E5 <- c(ccRCC_Ecotypes[["E5"]][["Fibroblast_S02"]][["gene"]],
        ccRCC_Ecotypes[["E5"]][["Mast_cells_S01"]][["gene"]],
        ccRCC_Ecotypes[["E5"]][["Myeloid_S04"]][["gene"]])

E6 <- c(ccRCC_Ecotypes[["E6"]][["Fibroblast_S04"]][["gene"]],
        ccRCC_Ecotypes[["E6"]][["Myeloid_S03"]][["gene"]],
        ccRCC_Ecotypes[["E6"]][["T_cells_S08"]][["gene"]])

E.list <- list("E1" = E1,
               "E2" = E2,
               "E3" = E3,
               "E4" = E4,
               "E5" = E5,
               "E6" = E6)

ccRCC_spaital_merge <- AddModuleScore(ccRCC_spaital_merge,features = E.list,name = names(E.list))

df <- ccRCC_spaital_merge@meta.data[,c("cluster","E11","E22","E33","E44","E55","E66")]


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )




ggplot(df, aes(x=cluster, y=E11,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E1_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E22,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E2_spaital.pdf",height = 3,width = 4)


ggplot(df, aes(x=cluster, y=E33,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E3_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E44,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E4_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E55,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E5_spaital.pdf",height = 3,width = 4)

ggplot(df, aes(x=cluster, y=E66,fill = cluster)) + 
  geom_boxplot()+
  scale_fill_brewer(type ="qual",palette = "Set3")+blank_theme
ggsave("D:/Urinary system tumors/work/5_Ecotype_spaital/ccRCC/E6_spaital.pdf",height = 3,width = 4)

