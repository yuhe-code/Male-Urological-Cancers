library(DoubletFinder)
library(dplyr)
library(Seurat)
library(harmony)
library(devtools)
library(SingleR)
#Seurat ------
{
  setwd(".\\1_ccRCC_tumor")
  
  assays <- dir(".\\scRNA-seq")
  
  samples_name = c('ccRCC_tumour1','ccRCC_tumour2','ccRCC_tumour3',"ccRCC_tumour4","ccRCC_tumour5",
                   "ccRCC_tumour6",'ccRCC_tumour7','ccRCC_tumour8','ccRCC_tumour9',"ccRCC_tumour10",
                   "ccRCC_tumour11",'ccRCC_tumour12','ccRCC_tumour13','ccRCC_tumour14',"ccRCC_tumour15",
                   "ccRCC_tumour16",'ccRCC_tumour17')
  GSM_name = c("GSM6290520","GSM6290522","GSM6290524","GSM6290526","GSM6290528",
               "GSM6290530","GSM6290532","GSM6290534","GSM6290536","GSM6290538",
               "GSM6290540","GSM6290542","GSM6290544","GSM6290548","GSM6290550",
               "GSM6290552","GSM6290554")
  
  dir <- paste0(".\\scRNA-seq\\",assays)
  
  scelist <- list()
  
  
  for (i in 1:length(assays)) {
    counts <- Read10X(data.dir = dir[i])
    scelist[[i]] <-CreateSeuratObject(counts, project=samples_name[i],min.cells=3, min.features = 200)
    
    scelist[[i]] <- RenameCells(scelist[[i]], add.cell.id = samples_name[i])
    scelist[[i]][["orig.ident"]] <- samples_name[i]
    scelist[[i]][["GSM_name"]] <- GSM_name[i]
    
    if(T){
      scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^MT-")
    }

    if(T){
      scelist[[i]][["percent.rb"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^RP[SL]")
    }
    if(T){
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
      HB.genes <- CaseMatch(HB.genes, rownames(scelist[[i]]))
      scelist[[i]][["percent.HB"]] <- PercentageFeatureSet(scelist[[i]], features = HB.genes)
    }
  }
  
  names(scelist) <- samples_name
  
  saveRDS(scelist,".\\ccRCC_Tumour_scelist_ori.rds") 
  
  
}

{
  
  scelist <- readRDS(".\\ccRCC_Tumour_scelist_ori.rds")
  
  vln_plot_list <- list()
  for (i in names(scelist)) {
    
    vln_plot_list[[i]]<- VlnPlot(scelist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                                 ncol = 3)
  }
  
  
  scelist_rm_double <- list()
  
  for (i in 1:length(scelist)) {
    
    sce <- subset(scelist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 300&percent.HB<5)
    sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData( )
    sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)
    

    ElbowPlot(sce, ndims = 50)
    pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1  
    co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2  
    pc_num <- min(co1, co2)
    sce <- sce %>% RunTSNE(dims = 1:pc_num) %>% RunUMAP(dims = 1:pc_num)
    sce <- FindNeighbors(sce, dims = 1:pc_num) %>% FindClusters(resolution = c(1:10/10)) 
    
    
    sweep.res.list <- paramSweep_v3(sce, PCs = 1:pc_num, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)] %>% as.character() %>% as.numeric()
    
    DoubletRate= ncol(sce)*8*1e-6
    
    homotypic.prop <- modelHomotypic(sce$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(sce)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
    
    
    sce <- doubletFinder_v3(sce, PCs= 1:pc_num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(sce@meta.data)[ncol(sce@meta.data)] <- "DF.classification"
    colnames(sce@meta.data)[ncol(sce@meta.data) - 1] <- "pANN"
    
    scelist_rm_double[[i]] <- sce
    
    print(paste0("______!!!!sample: ",i,"!!!!!______"))
    
    
  }
  
  
  sce <- merge(scelist_rm_double[[1]], scelist_rm_double[2:length(scelist_rm_double)])
  
  Idents(sce)="DF.classification"
  
  sce <- subset(sce,idents = "Singlet")
  
  saveRDS(sce, ".\\1_ccRCC_tumor_sce_afterQC.rds")
  
  
  sce_harmony <- sce
  
  sce_harmony <- NormalizeData(sce_harmony) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData( )
  sce_harmony <- RunPCA(sce_harmony, features = VariableFeatures(object = sce_harmony), verbose = F)
  
  sce_harmony <- RunHarmony(sce_harmony, group.by.vars = "orig.ident",
                            assay.use = "RNA", 
                            max.iter.harmony = 20, 
                            lambda = 1,
                            plot_convergence = TRUE)
  
  
  harmony_embeddings  <-  Embeddings(sce_harmony, 'harmony')
  DimPlot(object = sce_harmony, reduction = "harmony", group.by = "orig.ident",raster=FALSE)
  VlnPlot(object = sce_harmony, features = "harmony_1", group.by = "orig.ident",raster=FALSE)
  
  ElbowPlot(sce_harmony, ndims = 50,reduction = "harmony")

  pct <- sce_harmony [["pca"]]@stdev / sum( sce_harmony [["pca"]]@stdev)* 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1  
  co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2  
  pca_num <- min(co1,co2)
  
  sce_harmony <- FindNeighbors(sce_harmony,reduction = "harmony") %>% FindClusters(resolution = c(1:10/10))
  sce_harmony <- sce_harmony %>% RunTSNE(reduction = "harmony", dims = 1:pca_num) %>% RunUMAP(reduction = "harmony", dims = 1:pca_num) 
  
  colnames(sce_harmony@meta.data)
  library(clustree)
  p1 <- clustree(sce_harmony, prefix = 'RNA_snn_res.')+coord_flip()
  p1
  
  Idents(sce_harmony) <- "RNA_snn_res.0.6"
  DimPlot(sce_harmony, label = T, reduction = "tsne",raster=FALSE)
  DimPlot(sce_harmony, label = T, reduction = "umap",raster=FALSE)
  
  
  
  saveRDS(sce_harmony, ".\\2_ccRCC_tumor_sce_aftercluster.rds")

  Idents(sce_harmony) <- "RNA_snn_res.0.6"
  testdata <- GetAssayData(sce_harmony, slot="data")
  clusters <- sce_harmony@meta.data$RNA_snn_res.0.6
  
  ref <- readRDS(".\\singleR_human_ref.rds")
  
  cellpred_main <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                           clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  plotScoreHeatmap(cellpred_main,
                   max.labels = 8,
                   clusters = sce_harmony$RNA_snn_res.0.6,
                   order.by = c("labels"),
                   show.labels = F,
                   cluster_cols = T,
                   cluster_rows = T,
                   show_colnames = T)
  clusterAnn_main=as.data.frame(cellpred_main)
  celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype_main=cellpred_main$labels, stringsAsFactors = FALSE)
  
  sce_harmony@meta.data$singleR_ann ="NA"
  
  for(i in 1:nrow(celltype_main)){
    sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.6 == celltype_main$ClusterID[i]),'singleR_ann'] <- celltype_main$celltype[i]
  }
  
  
  unique(sce_harmony@meta.data[,c("RNA_snn_res.0.6","singleR_ann")])
  DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_ann",raster=FALSE)
  DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_ann",raster=FALSE)

  
  cellpred_fine <- SingleR(test = testdata, ref = ref, labels = ref$label.fine,
                           clusters = clusters,
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  plotScoreHeatmap(cellpred_fine,
                   max.labels = 8,
                   order.by = c("labels"),
                   show.labels = F,
                   cluster_cols = T,
                   cluster_rows = T,
                   show_colnames = T)
  clusterAnn_fine=as.data.frame(cellpred_fine)
  celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype_fine=cellpred_fine$labels, stringsAsFactors = FALSE)
  sce_harmony@meta.data$singleR_fine_ann ="NA"
  for(i in 1:nrow(celltype_fine)){
    sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.6 == celltype_fine$ClusterID[i]),'singleR_fine_ann'] <- celltype_fine$celltype[i]
  }
  unique(sce_harmony@meta.data[,c("RNA_snn_res.0.6","singleR_fine_ann")])
  
  DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_fine_ann",repel = T,raster=FALSE)
  DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_fine_ann",repel = T,raster=FALSE)
  
  
  saveRDS(sce_harmony,".\\3_ccRCC_tumor_sce_after_singleR_harmony.rds")
 
  Idents(sce_harmony) <- "RNA_snn_res.0.6"
  sce.markers <- FindAllMarkers(object = sce_harmony,
                                only.pos = T, 
                                min.pct = 0.25,
                                logfc.threshold = 0.5)
  sce.markers$pct=sce.markers$pct.1/sce.markers$pct.2
  
  sce.markers_used <- sce.markers %>%
    dplyr::filter(p_val_adj<0.05&pct.1>pct.2&pct>1.5&avg_log2FC>1) %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  sce <- readRDS("3_ccRCC_tumor_sce_after_singleR_harmony.rds")
  sce@meta.data$seurat_clusters <- sce@meta.data$RNA_snn_res.0.6
  sce@meta.data$seurat_clusters <- factor(sce@meta.data$seurat_clusters,levels=c(0:(length(unique(sce@meta.data$RNA_snn_res.0.6))-1)))
  library(lattice)
  library(ggplot2)
  feature_gene<-list(T_cells=c('CD3E','CD2','CD3D'),
                     B_cells=c('MS4A1','CD79A','JCHAIN'),
                     Plasma=c('IGHA1','IGHG1'),
                     Myeloid=c('CD68','C1QA','LYZ','CD14'),
                     Mast_cells=c('FCER1A', 'SLC18A2', 'TPSB2'),
                     Cycling_cells=c('MKI67','TOP2A'),
                     Endothelium=c('VWF','RAMP2','PECAM1','CLDN5'),
                     Epithelium=c('EPCAM','KRT18','CD24','KRT19'),
                     Mesenchymal_cells=c('ACTA2','PDGFRB'),
                     Fibroblast=c('COL1A1','TAGLN','DCN')
  )
  feature_gene2<-as.vector(unlist(feature_gene))
  p <- DotPlot(object = sce,group.by = 'seurat_clusters',
               features =feature_gene2,cols = c("#bababa", "#ca0020"))+
    scale_y_discrete(position = "left") +
    scale_fill_manual( breaks = rev(names(colors)),values = colors)+
    coord_flip()+
    theme_bw()
  p
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(19,25)]<-"B_cells"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(0,1,12)]<-"T_cells"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(5,15)]<-"Fibroblast"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(4,7,9,10,14,21,24)]<-"Myeloid"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(2,3,8,11,17,20,26)]<-"Epithelium"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(6,13,18,22,23)]<-"Endothelium"
  sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(16)]<-"Mast_cells"
  saveRDS(object = sce,file = "./manual/4_ccRCC_tumor_sce_after_cell_type.rds")
  p1 <- DimPlot(sce,group.by = "cell_type",label = T, raster = F, reduction = "umap")
  p2 <- DimPlot(sce,group.by = "cell_type",label = T, raster = F, reduction = "tsne")
  ggsave(filename = "./manual/tumor_manual_cell_type_umap.pdf",p1,width = 8,height = 6)
  ggsave(filename = "./manual/tumor_manual_cell_type_tsne.pdf",p2,width = 8,height = 6)
  ggsave(filename = "./manual/tumor_manual_cell_type_anno.pdf",p,width = 8,height = 8)
  DimPlot(sce,group.by = "seurat_clusters",label = T, raster = F, reduction = "umap")
}

##BC
{
  
  assays <- dir(".\\BC_data")
  
  for (i in assays) {
    sample_count <- read.table(paste0(".\\BC_data\\",i),
                               header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    index <- order(rowMeans(sample_count[,-1]),decreasing = T)
    exprSet_ordered <- sample_count[index,]
    keep <- !duplicated(exprSet_ordered$Name)
    exprSet_max <- exprSet_ordered[keep,]
    rownames(exprSet_max) <- exprSet_max[,1]
    exprSet_max <- exprSet_max[,-1]
    sce <- CreateSeuratObject(counts = exprSet_max)
    saveRDS(sce,paste0(".\\BC_seurat\\",i,".rds"))
    print(i)
    rm(sample_count,exprSet_max,sce)
    gc()
    
  }
  
  
}


{
  
  assays <- dir(".\\BC_seurat\\")
  
  dir <- paste0(".\\BC_seurat\\", assays)
  scelist <- list()
  for (i in 1:length(dir)) {
    scelist[[i]] <- readRDS(dir[i])
  }
  
  samples_name = c('BC_tumour1','BC_tumour2','BC_tumour3',"BC_tumour4","BC_tumour5","BC_tumour6")
  GSM_name = c("GSM6919778","GSM6919782","GSM6919784","GSM6919786","GSM7898117","GSM7898118")
  names(scelist) <- samples_name
  
  for(i in 1:length(scelist)){
    scelist[[i]] <- RenameCells(scelist[[i]], add.cell.id = samples_name[i])
    scelist[[i]][["orig.ident"]] <- samples_name[i]
    scelist[[i]][["GSM_name"]] <- GSM_name[i]
    scelist[[i]][["sample"]] <- GSM_name[i]
    
    if(T){
      scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^MT-")
    }
    if(T){
      scelist[[i]][["percent.rb"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^RP[SL]")
    }
    if(T){
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
      HB.genes <- CaseMatch(HB.genes, rownames(scelist[[i]]))
      scelist[[i]][["percent.HB"]] <- PercentageFeatureSet(scelist[[i]], features = HB.genes)
    }
  }
  
  vln_plot_list <- list()
  for (i in names(scelist)) {
    
    vln_plot_list[[i]]<- VlnPlot(scelist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                                 ncol = 3)
    
  }
  
  
  saveRDS(scelist,".\\2_BC_seuratlist\\BC_Tumour_scelist_ori.rds")
  
  
}

{
  scelist <- readRDS(".\\2_BC_seuratlist\\BC_Tumour_scelist_ori.rds")
  
  
  
  scelist_rm_double <- list()
  
  for (i in 1:length(scelist)) {
    
    sce <- subset(scelist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 300&percent.HB<5)
    sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData( )
    sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)
    
    ElbowPlot(sce, ndims = 50)
    pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1  
    co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2  
    pc_num <- min(co1, co2)
    sce <- sce %>% RunTSNE(dims = 1:pc_num) %>% RunUMAP(dims = 1:pc_num)
    sce <- FindNeighbors(sce, dims = 1:pc_num) %>% FindClusters(resolution = c(1:10/10)) 
    
    sweep.res.list <- paramSweep_v3(sce, PCs = 1:pc_num, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn_pbmc <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)] %>% as.character() %>% as.numeric()
    
    DoubletRate= ncol(sce)*8*1e-6
    
    homotypic.prop <- modelHomotypic(sce$seurat_clusters)
    nExp_poi <- round(DoubletRate*ncol(sce)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
    
    
    sce <- doubletFinder_v3(sce, PCs= 1:pc_num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
    colnames(sce@meta.data)[ncol(sce@meta.data)] <- "DF.classification"
    colnames(sce@meta.data)[ncol(sce@meta.data) - 1] <- "pANN"
    
    scelist_rm_double[[i]] <- sce
    
    print(paste0("______!!!!sample: ",i,"!!!!!______"))
    
    
  }
  
  
  sce <- merge(scelist_rm_double[[1]], scelist_rm_double[2:length(scelist_rm_double)])
  
  Idents(sce)="DF.classification"
  
  sce <- subset(sce,idents = "Singlet")
  
  saveRDS(sce, ".\\1_BC_Tumor_sce_afterQC.rds")
  
  
  sce_harmony <- sce
  
  sce_harmony <- NormalizeData(sce_harmony) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
  sce_harmony <- RunPCA(sce_harmony, features = VariableFeatures(object = sce_harmony), verbose = F)     
  
  sce_harmony <- RunHarmony(sce_harmony, group.by.vars = "orig.ident",
                            assay.use = "RNA", 
                            max.iter.harmony = 20, 
                            lambda = 1,
                            plot_convergence = TRUE)
  
  
  harmony_embeddings  <-  Embeddings(sce_harmony, 'harmony')
  DimPlot(object = sce_harmony, reduction = "harmony", group.by = "orig.ident")
  VlnPlot(object = sce_harmony, features = "harmony_1", group.by = "orig.ident")
  
  ElbowPlot(sce_harmony, ndims = 50,reduction = "pca")
  
  pct <- sce_harmony [["pca"]]@stdev / sum( sce_harmony [["pca"]]@stdev)* 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1  
  co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2 
  pc_num <- min(co1, co2)
  
  sce_harmony <- FindNeighbors(sce_harmony,reduction = "harmony") %>% FindClusters(resolution = c(1:10/10))
  sce_harmony <- sce_harmony %>% RunTSNE(reduction = "harmony", dims = 1:pc_num) %>% RunUMAP(reduction = "harmony", dims = 1:pc_num) 
  
  colnames(sce_harmony@meta.data)
  library(clustree)
  p1 <- clustree(sce_harmony, prefix = 'RNA_snn_res.')+coord_flip()
  p1
  
  Idents(sce_harmony) <- "RNA_snn_res.0.5"
  DimPlot(sce_harmony, label = T, reduction = "tsne")
  DimPlot(sce_harmony, label = T, reduction = "umap")
  
  
  
  saveRDS(sce_harmony, "2_BC_Tumor_sce_aftercluster.rds")
  
  Idents(sce_harmony) <- "RNA_snn_res.0.5"
  testdata <- GetAssayData(sce_harmony, slot="data")
  clusters <- sce_harmony@meta.data$RNA_snn_res.0.5
  
  ref <- readRDS("singleR_human_ref.rds")
  
  cellpred_main <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                           clusters = clusters, 
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  plotScoreHeatmap(cellpred_main,
                   max.labels = 8,
                   clusters = sce_harmony$RNA_snn_res.0.5,
                   order.by = c("labels"),
                   show.labels = F,
                   cluster_cols = T,
                   cluster_rows = T,
                   show_colnames = T)
  clusterAnn_main=as.data.frame(cellpred_main)
  celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype_main=cellpred_main$labels, stringsAsFactors = FALSE)
  
  sce_harmony@meta.data$singleR_ann ="NA"
  
  for(i in 1:nrow(celltype_main)){
    sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.5 == celltype_main$ClusterID[i]),'singleR_ann'] <- celltype_main$celltype[i]
  }
  
  
  unique(sce_harmony@meta.data[,c("RNA_snn_res.0.5","singleR_ann")])
  DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_ann")
  DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_ann")

  
  cellpred_fine <- SingleR(test = testdata, ref = ref, labels = ref$label.fine,
                           clusters = clusters,
                           assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  plotScoreHeatmap(cellpred_fine,
                   max.labels = 8,
                   #clusters = sce_harmony$RNA_snn_res.0.2,
                   order.by = c("labels"),
                   show.labels = F,
                   cluster_cols = T,
                   cluster_rows = T,
                   show_colnames = T)
  clusterAnn_fine=as.data.frame(cellpred_fine)
  celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype_fine=cellpred_fine$labels, stringsAsFactors = FALSE)
  sce_harmony@meta.data$singleR_fine_ann ="NA"
  for(i in 1:nrow(celltype_fine)){
    sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.5 == celltype_fine$ClusterID[i]),'singleR_fine_ann'] <- celltype_fine$celltype[i]
  }
  unique(sce_harmony@meta.data[,c("RNA_snn_res.0.5","singleR_fine_ann")])
  
  DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_fine_ann",repel = T)
  DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_fine_ann",repel = T)

  
  saveRDS(sce_harmony,"3_BC_Tumor_sce_after_singleR.rds")
  
  Idents(sce_harmony) <- "RNA_snn_res.0.5"
  sce.markers <- FindAllMarkers(object = sce_harmony,
                                only.pos = T, 
                                min.pct = 0.25,
                                logfc.threshold = 0.5)
  sce.markers$pct=sce.markers$pct.1/sce.markers$pct.2
  
  saveRDS(sce.markers,"./sce_BC_T_markers_large_up.rds")
  write.table(sce.markers,"./sce_BC_T_markers_large_up.xls",sep = "\t",quote = F)
  write.table(sce.markers,"./sce_BC_T_markers_large_up.txt",sep = "\t",quote = F)
  
  sce.markers_used <- sce.markers %>%
    dplyr::filter(p_val_adj<0.05&pct.1>pct.2&pct>1.5&avg_log2FC>1) %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  saveRDS(sce.markers,"./sce_BC_T_markers_fliter_up.rds")
  write.table(sce.markers_used,"./sce_BC_T_markers_fliter_up.xls",sep = "\t",quote = F)
  write.table(sce.markers_used,"./sce_BC_T_markers_fliter_up.txt",sep = "\t",quote = F)
  
  
  
  sce <- readRDS(".\\3_BC_Tumor_sce_after_singleR.rds")
  DimPlot(sce, label = T, reduction = "umap",group.by = "RNA_snn_res.0.5",repel = T)
  DimPlot(sce, label = T, reduction = "umap",group.by = "singleR_ann",repel = T)
  
  sce_markers_used <- read.table(".\\sce_BC_T_markers_fliter_up.txt")
  
  
  
  feature_gene<-list(T_cells=c('CD3E','CD2','CD3D'),
                     B_cells=c('MS4A1','CD79A','JCHAIN'),
                     Plasma=c('IGHA1','IGHG1'),
                     Myeloid=c('CD68','C1QA','LYZ','CD14'),
                     Mast_cells=c('FCER1A', 'SLC18A2', 'TPSB2'),
                     Cycling_cells=c('MKI67','TOP2A'),
                     Endothelium=c('VWF','RAMP2','PECAM1','CLDN5'),
                     Epithelium=c('EPCAM','KRT18','CD24','KRT19'),
                     Mesenchymal_cells=c('ACTA2','PDGFRB'),
                     Fibroblast=c('COL1A1','TAGLN','DCN')
  )
  feature_gene2<-as.vector(unlist(feature_gene))
  DotPlot(object = sce,group.by = 'RNA_snn_res.0.5',
          features =feature_gene2,cols = c("#bababa", "#ca0020"))+
    scale_y_discrete(position = "left") +
    scale_fill_manual( breaks = rev(names(colors)),values = colors)+
    coord_flip()+
    theme_bw()
  
  
  
}
Idents(sce) <- "RNA_snn_res.0.5"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(0,2,15)]<-"T_cells"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(1,9,17,21)]<-"Fibroblast"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(3,5,6,10,14,20)]<-"Epithelium"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(4,8,16,18,22)]<-"Endothelium"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(7,12,19)]<-"Myeloid"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(11,13)]<-"B_cells"
sce@meta.data$cell_type[sce@meta.data$RNA_snn_res.0.5%in%c(23)]<-"Unknow"
DimPlot(sce, label = T, reduction = "umap",group.by = "cell_type",repel = T)
DimPlot(sce, label = T, reduction = "tsne",group.by = "cell_type",repel = T)
saveRDS(object = sce,file = ".\\4_BC_tumor_sce_after_cell_type.rds")
p1 <- DimPlot(sce,group.by = "cell_type",label = T, raster = F, reduction = "umap")
p2 <- DimPlot(sce,group.by = "cell_type",label = T, raster = F, reduction = "tsne")
setwd(".\\1_BC_tumor_sce_ann\\manual")
ggsave(filename = "tumor_manual_cell_type_umap.pdf",p1,width = 8,height = 6)
ggsave(filename = "tumor_manual_cell_type_tsne.pdf",p2,width = 8,height = 6)

#PCa
{
  
  #1_GSE176031_17_8
  setwd(".\\3_PCa")
  
  assays <- dir(".\\1_GSE176031_17_8\\1_PCa_tumor_GSE176031")
  
  
  scelist <- list()
  
  for (i in assays) {
    sample_count <- read.table(paste0(".\\3_PCa\\1_GSE176031_17_8\\1_PCa_tumor_GSE176031\\",i),
                               row.names=1,head=T )
    
    sce <- CreateSeuratObject(counts = sample_count)
    scelist[[i]] <- sce
    print(i)
    rm(sample_count,sce)
    gc()
    
  }
  
  samples_name = c('PCa_tumour1','PCa_tumour2','PCa_tumour3',"PCa_tumour4","PCa_tumour5",
                   "PCa_tumour6",'PCa_tumour7','PCa_tumour8','PCa_tumour9',"PCa_tumour10",
                   "PCa_tumour11","PCa_tumour12","PCa_tumour13","PCa_tumour14","PCa_tumour15",
                   "PCa_tumour16","PCa_tumour17")
  GSM_name = c('GSM5353224','GSM5353225','GSM5353226',"GSM5353227","GSM5353228",
               "GSM5353229",'GSM5353232','GSM5353233','GSM5353236',"GSM5353237",
               "GSM5353240","GSM5353243","GSM5353244","GSM5353245","GSM5353246",
               "GSM5353247","GSM5353248")
  names(scelist) <- samples_name
  saveRDS(scelist,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist\\1_PCa_tumor_GSE176031_scelist_ori.rds")
  
  
  
  #2_GSE193337_4_4
  
  setwd(".\\3_PCa")
  
  assays <- dir(".\\2_GSE193337_4_4\\1_PCa_tumor_GSE193337")
  scelist <- list()
  
  for (i in assays) {
    sample_count <- Read10X(paste0(".\\3_PCa\\2_GSE193337_4_4\\1_PCa_tumor_GSE193337\\",i))
    sce <- CreateSeuratObject(counts = sample_count)
    scelist[[i]] <- sce
    print(i)
    rm(sample_count,sce)
    gc()
    
  }
  
  samples_name = c('PCa_tumour18','PCa_tumour19','PCa_tumour20',"PCa_tumour21")
  GSM_name = c('GSM5793828','GSM5793829','GSM5793831',"GSM5793832")
  names(scelist) <- samples_name
  saveRDS(scelist,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist\\2_PCa_tumor_GSE193337_scelist_ori.rds")
  
  
  #3_GSE157703_2
  setwd(".\\3_PCa")
  
  assays <- dir(".\\3_GSE157703_2")
  
  
  scelist <- list()
  
  for (i in assays) {
    sample_count <- read.table(paste0(".\\3_PCa\\3_GSE157703_2\\",i),
                               row.names=1,head=T )
    
    sce <- CreateSeuratObject(counts = sample_count)
    scelist[[i]] <- sce
    print(i)
    rm(sample_count,sce)
    gc()
    
  }
  
  samples_name = c('PCa_tumour22','PCa_tumour23')
  GSM_name = c('GSM4773521','GSM4773522')
  names(scelist) <- samples_name
  saveRDS(scelist,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist\\3_PCa_tumor_GSE157703_scelist_ori.rds")
  
  
}


{
  
  setwd(".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist")
  a_PCa_tumor_GSE176031_scelist <- readRDS(".\\1_PCa_tumor_GSE176031_scelist_ori.rds")
  b_PCa_tumor_GSE193337_scelist <- readRDS(".\\2_PCa_tumor_GSE193337_scelist_ori.rds")
  c_PCa_tumor_GSE157703_scelist <- readRDS(".\\3_PCa_tumor_GSE157703_scelist_ori.rds")
  
  for (i in 1:length(a_PCa_tumor_GSE176031_scelist)) {
    a_PCa_tumor_GSE176031_scelist[[i]][["GSE"]] <- "GSE176031"
    
  }
  
  for (i in 1:length(b_PCa_tumor_GSE193337_scelist)) {
    b_PCa_tumor_GSE193337_scelist[[i]][["GSE"]] <- "GSE193337"
    
  }
  
  for (i in 1:length(c_PCa_tumor_GSE157703_scelist)) {
    c_PCa_tumor_GSE157703_scelist[[i]][["GSE"]] <- "GSE157703"
    
  }
  
  PCa_tumor_list <- c(a_PCa_tumor_GSE176031_scelist,b_PCa_tumor_GSE193337_scelist,c_PCa_tumor_GSE157703_scelist)
  
  saveRDS(PCa_tumor_list,".\\PCa_tumor_all_scelist_ori.rds")
  
  
  
  
  scelist <- readRDS(".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist\\PCa_tumor_all_scelist_ori.rds")
  
  
  samples_name = c('PCa_tumour1','PCa_tumour2','PCa_tumour3',"PCa_tumour4","PCa_tumour5",
                   "PCa_tumour6",'PCa_tumour7','PCa_tumour8','PCa_tumour9',"PCa_tumour10",
                   "PCa_tumour11","PCa_tumour12","PCa_tumour13","PCa_tumour14","PCa_tumour15",
                   "PCa_tumour16","PCa_tumour17",'PCa_tumour18','PCa_tumour19','PCa_tumour20',
                   "PCa_tumour21",'PCa_tumour22','PCa_tumour23')
  
  GSM_name = c('GSM5353224','GSM5353225','GSM5353226',"GSM5353227","GSM5353228",
               "GSM5353229",'GSM5353232','GSM5353233','GSM5353236',"GSM5353237",
               "GSM5353240","GSM5353243","GSM5353244","GSM5353245","GSM5353246",
               "GSM5353247","GSM5353248",'GSM5793828','GSM5793829','GSM5793831',
               "GSM5793832",'GSM4773521','GSM4773522')
  
  
  
  
  for(i in 1:length(scelist)){
    scelist[[i]] <- RenameCells(scelist[[i]], add.cell.id = samples_name[i])
    scelist[[i]][["orig.ident"]] <- samples_name[i]
    scelist[[i]][["GSM_name"]] <- GSM_name[i]
    
    if(T){
      scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^MT-")
    }
    if(T){
      scelist[[i]][["percent.rb"]] <- PercentageFeatureSet(scelist[[i]],pattern = "^RP[SL]")
    }
    if(T){
      HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
      HB.genes <- CaseMatch(HB.genes, rownames(scelist[[i]]))
      scelist[[i]][["percent.HB"]] <- PercentageFeatureSet(scelist[[i]], features = HB.genes)
    }
  }
  
  vln_plot_list <- list()
  for (i in names(scelist)) {
    
    vln_plot_list[[i]]<- VlnPlot(scelist[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                                 ncol = 3)
    
  }
  
  
  saveRDS(scelist,".\\PCa_tumor_all_scelist_ori.rds")

  
  {
    scelist <- readRDS(".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_seuratlist\\PCa_tumor_all_scelist_ori.rds")
    
  
    
    
    scelist_rm_double <- list()
    
    for (i in 1:length(scelist)) {
      
      sce <- subset(scelist[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10 & nCount_RNA > 300&percent.HB<5)
      sce <- NormalizeData(sce) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData( )
      sce <- RunPCA(sce, features = VariableFeatures(object = sce), verbose = F)
      
      ElbowPlot(sce, ndims = 50)
      pct <- sce [["pca"]]@stdev / sum( sce [["pca"]]@stdev)* 100
      cumu <- cumsum(pct)
      co1 <- which(cumu > 90 & pct < 5)[1]
      co1  
      co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
      co2  
      pc_num <- min(co1, co2)
      sce <- sce %>% RunTSNE(dims = 1:pc_num) %>% RunUMAP(dims = 1:pc_num)
      sce <- FindNeighbors(sce, dims = 1:pc_num) %>% FindClusters(resolution = c(1:10/10)) 
      
      
      sweep.res.list <- paramSweep_v3(sce, PCs = 1:pc_num, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn_pbmc <- find.pK(sweep.stats)
      pK_bcmvn <- bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)] %>% as.character() %>% as.numeric()
      
      DoubletRate= ncol(sce)*8*1e-6
      
      homotypic.prop <- modelHomotypic(sce$seurat_clusters)
      nExp_poi <- round(DoubletRate*ncol(sce)) 
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
      
      
      sce <- doubletFinder_v3(sce, PCs= 1:pc_num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
      colnames(sce@meta.data)[ncol(sce@meta.data)] <- "DF.classification"
      colnames(sce@meta.data)[ncol(sce@meta.data) - 1] <- "pANN"
      
      scelist_rm_double[[i]] <- sce
      
      print(paste0("______!!!!sample: ",i,"!!!!!______"))
      
      
    }
    
    
    sce <- merge(scelist_rm_double[[1]], scelist_rm_double[2:length(scelist_rm_double)])
    
    Idents(sce)="DF.classification"
    
    sce <- subset(sce,idents = "Singlet")
    
    saveRDS(sce, ".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\1_PCa_Tumor_sce_afterQC.rds")
    dim(sce)
    
    sce_harmony <- sce
    
    sce_harmony <- NormalizeData(sce_harmony) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
    sce_harmony <- RunPCA(sce_harmony, features = VariableFeatures(object = sce_harmony), verbose = F)     
    
    sce_harmony <- RunHarmony(sce_harmony, group.by.vars = "orig.ident",
                              assay.use = "RNA", 
                              max.iter.harmony = 20, 
                              lambda = 1,
                              plot_convergence = TRUE)
    
    
    harmony_embeddings  <-  Embeddings(sce_harmony, 'harmony')
    DimPlot(object = sce_harmony, reduction = "harmony", group.by = "orig.ident")
    VlnPlot(object = sce_harmony, features = "harmony_1", group.by = "orig.ident")
    
    ElbowPlot(sce_harmony, ndims = 50,reduction = "pca")
    
    pct <- sce_harmony [["pca"]]@stdev / sum( sce_harmony [["pca"]]@stdev)* 100
    cumu <- cumsum(pct)
    co1 <- which(cumu > 90 & pct < 5)[1]
    co1  
    co2 <-sort(which((pct[1:length(pct) - 1]- pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    co2 
    pc_num <- min(co1, co2)
    
    sce_harmony <- FindNeighbors(sce_harmony,reduction = "harmony") %>% FindClusters(resolution = c(1:10/10))
    sce_harmony <- sce_harmony %>% RunTSNE(reduction = "harmony", dims = 1:pc_num) %>% RunUMAP(reduction = "harmony", dims = 1:pc_num) 
    
    colnames(sce_harmony@meta.data)
    library(clustree)
    p1 <- clustree(sce_harmony, prefix = 'RNA_snn_res.')+coord_flip()
    p1
    
    Idents(sce_harmony) <- "RNA_snn_res.0.7"
    DimPlot(sce_harmony, label = T, reduction = "tsne")
    DimPlot(sce_harmony, label = T, reduction = "umap")
    
    
    
    saveRDS(sce_harmony, ".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\2_PCa_Tumor_sce_aftercluster.rds")
    
    Idents(sce_harmony) <- "RNA_snn_res.0.7"
    testdata <- GetAssayData(sce_harmony, slot="data")
    clusters <- sce_harmony@meta.data$RNA_snn_res.0.7
    
    ref <- readRDS(".\\1_scRNA_landsacape\\singleR_ref\\singleR_human_ref.rds")
    
    cellpred_main <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                             clusters = clusters, 
                             assay.type.test = "logcounts", assay.type.ref = "logcounts")
    plotScoreHeatmap(cellpred_main,
                     max.labels = 8,
                     clusters = sce_harmony$RNA_snn_res.0.5,
                     order.by = c("labels"),
                     show.labels = F,
                     cluster_cols = T,
                     cluster_rows = T,
                     show_colnames = T)
    clusterAnn_main=as.data.frame(cellpred_main)
    celltype_main = data.frame(ClusterID=rownames(cellpred_main), celltype_main=cellpred_main$labels, stringsAsFactors = FALSE)
    
    sce_harmony@meta.data$singleR_ann ="NA"
    
    for(i in 1:nrow(celltype_main)){
      sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.7 == celltype_main$ClusterID[i]),'singleR_ann'] <- celltype_main$celltype[i]
    }
    
    
    unique(sce_harmony@meta.data[,c("RNA_snn_res.0.7","singleR_ann")])
    DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_ann")
    DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_ann")
    
    cellpred_fine <- SingleR(test = testdata, ref = ref, labels = ref$label.fine,
                             clusters = clusters,
                             assay.type.test = "logcounts", assay.type.ref = "logcounts")
    
    plotScoreHeatmap(cellpred_fine,
                     max.labels = 8,
                     order.by = c("labels"),
                     show.labels = F,
                     cluster_cols = T,
                     cluster_rows = T,
                     show_colnames = T)
    clusterAnn_fine=as.data.frame(cellpred_fine)
    celltype_fine = data.frame(ClusterID=rownames(cellpred_fine), celltype_fine=cellpred_fine$labels, stringsAsFactors = FALSE)
    sce_harmony@meta.data$singleR_fine_ann ="NA"
    for(i in 1:nrow(celltype_fine)){
      sce_harmony@meta.data[which(sce_harmony@meta.data$RNA_snn_res.0.7 == celltype_fine$ClusterID[i]),'singleR_fine_ann'] <- celltype_fine$celltype[i]
    }
    unique(sce_harmony@meta.data[,c("RNA_snn_res.0.7","singleR_fine_ann")])
    
    DimPlot(sce_harmony, label = T, reduction = "tsne",group.by = "singleR_fine_ann",repel = T)
    DimPlot(sce_harmony, label = T, reduction = "umap",group.by = "singleR_fine_ann",repel = T)
    saveRDS(sce_harmony,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\3_PCa_Tumor_sce_after_singleR.rds")
    
    
    Idents(sce_harmony) <- "RNA_snn_res.0.7"
    sce.markers <- FindAllMarkers(object = sce_harmony,
                                  only.pos = T, 
                                  min.pct = 0.25,
                                  logfc.threshold = 0.5)
    sce.markers$pct=sce.markers$pct.1/sce.markers$pct.2
    saveRDS(sce.markers,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_large_up.rds")
    write.table(sce.markers,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_large_up.xls",sep = "\t",quote = F)
    write.table(sce.markers,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_large_up.txt",sep = "\t",quote = F)
    sce.markers_used <- sce.markers %>%
      dplyr::filter(p_val_adj<0.05&pct.1>pct.2&pct>1.5&avg_log2FC>1) %>%
      group_by(cluster) %>%
      slice_max(n = 20, order_by = avg_log2FC)
    
    saveRDS(sce.markers_used,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_fliter_up.rds")
    write.table(sce.markers_used,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_fliter_up.xls",sep = "\t",quote = F)
    write.table(sce.markers_used,".\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce_ann\\sce_PCa_T_markers_fliter_up.txt",sep = "\t",quote = F)
    feature_gene<-list(T_cells=c('CD3E','CD2','CD3D'),
                       B_cells=c('MS4A1','CD79A','JCHAIN'),
                       Plasma=c('IGHA1','IGHG1'),
                       Myeloid=c('CD68','C1QA','LYZ','CD14'),
                       Mast_cells=c('FCER1A', 'SLC18A2', 'TPSB2'),
                       Cycling_cells=c('MKI67','TOP2A'),
                       Endothelium=c('VWF','RAMP2','PECAM1','CLDN5'),
                       Epithelium=c('EPCAM','KRT18','CD24','KRT19'),
                       Mesenchymal_cells=c('ACTA2','PDGFRB'),
                       Fibroblast=c('COL1A1','TAGLN','DCN')
    )
    feature_gene2<-as.vector(unlist(feature_gene))
    DotPlot(object = sce,group.by = 'RNA_snn_res.0.7',
            features =feature_gene2,cols = c("#bababa", "#ca0020"))+
      scale_y_discrete(position = "left") +
      scale_fill_manual( breaks = rev(names(colors)),values = colors)+
      coord_flip()+
      theme_bw()
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(15)]<-"B_cells"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(0,2,17)]<-"T_cells"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(11,12,20)]<-"Fibroblast"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(1,13)]<-"Myeloid"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(4,5,6,7,8,9,10,18)]<-"Epithelium"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(3,16,19)]<-"Endothelium"
    sce@meta.data$cell_type[sce@meta.data$seurat_clusters%in%c(14)]<-"Mast_cells"
    
  }
}

#cell plot
library(Seurat)
library(ggplot2)
library(patchwork)
#ccRCC

#tumor
{
  
  
  sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
  table(sce_ccRCC_tumor@meta.data$cell_type)
  color_ann <- c("Epithelium" = "#E45A5F",
                 "T_cells" = "#41B749",
                 "B_cells" = "#684797",
                 "Myeloid" = "#1781b5",
                 "Fibroblast" = "#F58135",
                 "Endothelium" = "#C0937E",
                 "Mast_cells"="#ce5e8a")
  order_ann <- rev(c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast","Mast_cells"))
  
  DimPlot(sce_ccRCC_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann,raster=FALSE)+ 
    theme(legend.position = "right",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_ccRCC_tumor_umap_cell_type_legend.pdf"
         ,width = 7, height = 5)
  dev.off()
  DimPlot(sce_ccRCC_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann,raster=FALSE)+ 
    theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_ccRCC_tumor_umap_cell_type.pdf"
         ,width = 5, height = 5)
  dev.off()
  
}


#BC

#tumor
{
  sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")

  color_ann <- c("Epithelium" = "#E45A5F",
                 "T_cells" = "#41B749",
                 "B_cells" = "#684797",
                 "Myeloid" = "#1781b5",
                 "Fibroblast" = "#F58135",
                 "Endothelium" = "#C0937E")
  order_ann <- rev(c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast"))
  
  DimPlot(sce_BC_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann)+ 
    theme(legend.position = "right",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_BC_tumor_umap_cell_type_legend.pdf"
         ,width = 7, height = 5)
  dev.off()
  
  
  
  DimPlot(sce_BC_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann)+ 
    theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_BC_tumor_umap_cell_type.pdf"
         ,width = 5, height = 5)
  dev.off()
  
  
}

#PCa
#tumor
{
  sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")
  
  
  table(sce_PCa_tumor@meta.data$cell_type)
  color_ann <- c("Epithelium" = "#E45A5F",
                 "T_cells" = "#41B749",
                 "B_cells" = "#684797",
                 "Myeloid" = "#1781b5",
                 "Fibroblast" = "#F58135",
                 "Endothelium" = "#C0937E",
                 "Mast_cells"="#ce5e8a")
  order_ann <- rev(c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast","Mast_cells"))
  
  DimPlot(sce_PCa_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann)+ 
    theme(legend.position = "right",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_PCa_tumor_umap_cell_type_legend.pdf"
         ,width = 7, height = 5)
  dev.off()
  
  
  
  DimPlot(sce_PCa_tumor, label = F,reduction = "umap",group.by = "cell_type",repel = T,cols = color_ann,order = order_ann)+ 
    theme(legend.position = "none",axis.title.y=element_blank(), axis.text.y=element_blank(),axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.line=element_blank(),axis.ticks=element_blank(),plot.title = element_blank() )
  ggsave(filename = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\1_PCa_tumor_umap_cell_type.pdf"
         ,width = 5, height = 5)
  dev.off()
  
}


{
sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")

sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")


sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")


Idents(sce_ccRCC_tumor) <- "cell_type"

Idents(sce_BC_tumor) <- "cell_type"

Idents(sce_PCa_tumor) <- "cell_type"



Cellratio_ccRCC_tumor <- prop.table(table(Idents(sce_ccRCC_tumor), sce_ccRCC_tumor$orig.ident), margin = 2)

Cellratio_ccRCC_tumor <- as.data.frame(Cellratio_ccRCC_tumor)




Cellratio_BC_tumor <- prop.table(table(Idents(sce_BC_tumor), sce_BC_tumor$orig.ident), margin = 2)

Cellratio_BC_tumor <- as.data.frame(Cellratio_BC_tumor)



Cellratio_PCa_tumor <- prop.table(table(Idents(sce_PCa_tumor), sce_PCa_tumor$orig.ident), margin = 2)

Cellratio_PCa_tumor <- as.data.frame(Cellratio_PCa_tumor)

Cellratio_all <- rbind(Cellratio_ccRCC_tumor,
                       Cellratio_BC_tumor,
                       Cellratio_PCa_tumor)



save(Cellratio_all,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_all.Rdata")

Cellratio_ccRCC_tumor_all <- prop.table(table(Idents(sce_ccRCC_tumor)))
Cellratio_ccRCC_tumor_all <- as.data.frame(Cellratio_ccRCC_tumor_all)
save(Cellratio_ccRCC_tumor_all,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_ccRCC_tumor_all.Rdata")



Cellratio_BC_tumor_all <- prop.table(table(Idents(sce_BC_tumor)))
Cellratio_BC_tumor_all <- as.data.frame(Cellratio_BC_tumor_all)
save(Cellratio_BC_tumor_all,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_BC_tumor_all.Rdata")


Cellratio_PCa_tumor_all <- prop.table(table(Idents(sce_PCa_tumor)))
Cellratio_PCa_tumor_all <- as.data.frame(Cellratio_PCa_tumor_all)
save(Cellratio_PCa_tumor_all,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_PCa_tumor_all.Rdata")


}


{
  load(file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_all.Rdata")
  
  
  
  color_ann <- c("Epithelium" = "#E45A5F",
                 "T_cells" = "#41B749",
                 "B_cells" = "#684797",
                 "Myeloid" = "#1781b5",
                 "Fibroblast" = "#F58135",
                 "Endothelium" = "#C0937E",
                 "Mast_cells"="#ce5e8a")
  
  
  
  
  
  library(ggplot2)
  
  Cellratio_all$Var1 <- factor(Cellratio_all$Var1,
                               levels =rev(c("Epithelium","T_cells","B_cells","Myeloid",
                                             "Endothelium","Fibroblast","Mast_cells")))
  
  Cellratio_all$Var2 <- factor(Cellratio_all$Var2,
                               levels =rev(c(
                                 "ccRCC_tumour1", "ccRCC_tumour2","ccRCC_tumour3","ccRCC_tumour4",
                                 "ccRCC_tumour5","ccRCC_tumour6","ccRCC_tumour7","ccRCC_tumour8", "ccRCC_tumour9","ccRCC_tumour10","ccRCC_tumour11",
                                 "ccRCC_tumour12","ccRCC_tumour13","ccRCC_tumour14","ccRCC_tumour15","ccRCC_tumour16","ccRCC_tumour17",
                                 "BC_tumour1","BC_tumour2","BC_tumour3","BC_tumour4","BC_tumour5","BC_tumour6",
                                 "PCa_tumour1","PCa_tumour2","PCa_tumour3","PCa_tumour4",
                                 "PCa_tumour5","PCa_tumour6","PCa_tumour7","PCa_tumour8","PCa_tumour9","PCa_tumour10","PCa_tumour11","PCa_tumour12",
                                 "PCa_tumour13","PCa_tumour14","PCa_tumour15","PCa_tumour16","PCa_tumour17","PCa_tumour18","PCa_tumour19","PCa_tumour20",
                                 "PCa_tumour21","PCa_tumour22","PCa_tumour23"
                                 
                               )))
  
  p_all <-  ggplot(Cellratio_all) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
    theme_classic() +
    scale_fill_manual(values=color_ann)+
    labs(x='Sample',y = 'Percentage')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+
    theme(legend.position = "none")
  ggsave(p_all,file="D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_all.pdf"
         ,width = 3,height = 10)
  
  
  
}

{
  Cellratio_ccRCC <- Cellratio_ccRCC_tumor
  
  Cellratio_ccRCC$Var1 <- factor(Cellratio_ccRCC$Var1,
                                 levels =rev(c("Epithelium","T_cells","B_cells","Myeloid",
                                               "Endothelium","Fibroblast","Mast_cells")))

  
  Cellratio_ccRCC$Var2 <- factor(Cellratio_ccRCC$Var2,
                                 levels =rev(c("ccRCC_tumour1", "ccRCC_tumour2","ccRCC_tumour3","ccRCC_tumour4",
                                               "ccRCC_tumour5","ccRCC_tumour6","ccRCC_tumour7","ccRCC_tumour8", "ccRCC_tumour9","ccRCC_tumour10","ccRCC_tumour11",
                                               "ccRCC_tumour12","ccRCC_tumour13","ccRCC_tumour14","ccRCC_tumour15","ccRCC_tumour16","ccRCC_tumour17"
                                 )))
  
  
  save(Cellratio_ccRCC,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_ccRCC.Rdata")
  p1 <- ggplot(Cellratio_ccRCC) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
    theme_classic() +
    scale_fill_manual(values=color_ann)+
    labs(x='Sample',y = 'Percentage')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+  
    theme(legend.position="none")
  p1 
  ggsave(p1,file="D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_ccRCC.pdf"
         ,width = 5,height = 0.05*length(Cellratio_ccRCC$Var2))
  
}

{
  
  Cellratio_BC <- Cellratio_BC_tumor
  
  Cellratio_BC$Var1 <- factor(Cellratio_BC$Var1,
                              levels =rev(c("Epithelium","T_cells","B_cells","Myeloid",
                                            "Endothelium","Fibroblast","Mast_cells")))
  Cellratio_BC$Var2 <- factor(Cellratio_BC$Var2,
                              levels =rev(c(
                                "BC_tumour1","BC_tumour2","BC_tumour3","BC_tumour4","BC_tumour5","BC_tumour6"
                              )))
  
  p2 <- ggplot(Cellratio_BC) + 
    geom_bar(x =Cellratio_BC$Var2, y= Cellratio_BC$Freq, fill = Cellratio_BC$Var1,stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
    theme_classic() +
    scale_fill_manual(values=color_ann)+
    labs(x='Sample',y = 'Percentage')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+  
    theme(legend.position="none")
  p2
  
  ggsave(p2,file="D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_BC.pdf"
         ,width  = 5,  height = 0.05*length(Cellratio_BC$Var2))
  
  
}

{
  Cellratio_PCa <- rbind(Cellratio_PCa_normal,Cellratio_PCa_tumor)
  
  Cellratio_PCa$Var1 <- factor(Cellratio_PCa$Var1,
                               levels =rev(c("Epithelium","T_cells","B_cells","Myeloid",
                                             "Endothelium","Fibroblast","Mast_cells")))
  Cellratio_PCa$Var2 <- factor(Cellratio_PCa$Var2,
                               levels =rev(c("PCa_normal1","PCa_normal2","PCa_normal3","PCa_normal4","PCa_normal5","PCa_normal6","PCa_normal7","PCa_normal8",
                                             "PCa_normal9","PCa_normal10","PCa_normal11","PCa_normal12","PCa_tumour1","PCa_tumour2","PCa_tumour3","PCa_tumour4",
                                             "PCa_tumour5","PCa_tumour6","PCa_tumour7","PCa_tumour8","PCa_tumour9","PCa_tumour10","PCa_tumour11","PCa_tumour12",
                                             "PCa_tumour13","PCa_tumour14","PCa_tumour15","PCa_tumour16","PCa_tumour17","PCa_tumour18","PCa_tumour19","PCa_tumour20",
                                             "PCa_tumour21","PCa_tumour22","PCa_tumour23"
                               )))
  
  save(Cellratio_PCa,file = "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\2_Cellratio_PCa.Rdata")
  p3 <- ggplot(Cellratio_PCa) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = 'white')+ 
    theme_classic() +
    scale_fill_manual(values=color_ann)+
    labs(x='Sample',y = 'Percentage')+
    coord_flip()+
    theme(panel.border = element_rect(fill=NA,color="gray", size=1, linetype="solid"))+  
    theme(legend.position="none")
  
  
}


library(ggplot2)
feature_gene<-list(T_cells=c('CD3E','CD2','CD3D'),
                   B_cells=c('MS4A1','CD79A','JCHAIN'),
                   Plasma=c('IGHA1','IGHG1'),
                   Myeloid=c('CD68','C1QA','LYZ','CD14'),
                   Mast_cells=c('FCER1A', 'SLC18A2', 'TPSB2'),
                   Cycling_cells=c('MKI67','TOP2A'),
                   Endothelium=c('VWF','RAMP2','PECAM1','CLDN5'),
                   Epithelium=c('EPCAM','KRT18','CD24','KRT19'),
                   Mesenchymal_cells=c('ACTA2','PDGFRB'),
                   Fibroblast=c('COL1A1','TAGLN','DCN')
)


marker_genes <- rev(c('EPCAM','KRT18','CD24','KRT19', #Epithelium
                      'CD3E','CD2','CD3D',  #T_cells
                      'MS4A1','CD79A','JCHAIN', #B_cells
                      'CD68','C1QA','LYZ','CD14', #Myeloid
                      'VWF','RAMP2','PECAM1','CLDN5', #Endothelium
                      'COL1A1','TAGLN','DCN', #Fibroblast
                      'FCER1A', 'SLC18A2', 'TPSB2' #Mast_cells
))


#ccRCC_tumor
{
  sce_ccRCC_tumor@meta.data$cell_type <- factor(sce_ccRCC_tumor@meta.data$cell_type,levels = c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast","Mast_cells"))
  
  exp_ccRCC_tumor <- DotPlot(sce_ccRCC_tumor,features = marker_genes,group.by = "cell_type")$data
  min(exp_ccRCC_tumor$avg.exp.scaled)
  max(exp_ccRCC_tumor$avg.exp.scaled)
  #min:-0.6094959  max: 2.267787
  
  
  
  p2_ccRCC_tumor <- DotPlot(sce_ccRCC_tumor,features = marker_genes,group.by = "cell_type")+
    coord_flip()+theme_bw()+
    theme(legend.position = "none",panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
    scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-0.7231908,2.267787)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
  
}


#BC_tumor
{
  sce_BC_tumor@meta.data$cell_type <- factor(sce_BC_tumor@meta.data$cell_type,levels = c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast","Mast_cells"))
  
  exp_BC_tumor <- DotPlot(sce_BC_tumor,features = marker_genes,group.by = "cell_type")$data
  min(exp_BC_tumor$avg.exp.scaled)
  max(exp_BC_tumor$avg.exp.scaled)
  #min:-0.6511926  max: 2.041227
  
  
  
  p4_BC_tumor <- DotPlot(sce_BC_tumor,features = marker_genes,group.by = "cell_type")+
    coord_flip()+theme_bw()+
    theme(legend.position = "none",panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
    scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-0.7231908,2.267787)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
  
}


#PCa_tumor
{
  sce_PCa_tumor@meta.data$cell_type <- factor(sce_PCa_tumor@meta.data$cell_type,levels = c("Epithelium","T_cells","B_cells","Myeloid","Endothelium","Fibroblast","Mast_cells"))
  
  exp_PCa_tumor <- DotPlot(sce_PCa_tumor,features = marker_genes,group.by = "cell_type")$data
  min(exp_PCa_tumor$avg.exp.scaled)
  max(exp_PCa_tumor$avg.exp.scaled)
  #min:-0.6268526  max: 2.267786
  
  
  
  p6_PCa_tumor <- DotPlot(sce_PCa_tumor,features = marker_genes,group.by = "cell_type")+
    coord_flip()+theme_bw()+
    theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45,hjust = 1,vjust=1.0))+
    scale_color_gradientn(values = seq(0,1,0.2),limits=(c(-0.7231908,2.267787)),colours = c("#8A959B","#C8CABF","#fa7e23","#ED6766"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
  
}

p_all <- p1_ccRCC_normal+p2_ccRCC_tumor+p3_BC_normal+p4_BC_tumor+p5_PCa_normal+p6_PCa_tumor+plot_layout(nrow = 1, byrow = T)

ggsave(p_all,file="D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\result1_plot\\3_dot_cell_type_marker.pdf",
       width = 15,height = 5)

#Cell disteny-------------
PCa_sce_anno_second <- readRDS("D:/Urinary system tumors/work/4_merge/PCa_sce_anno_second.rds")
PCa_sce_anno_second@meta.data$G_group <- "none"
PCa_sce_anno_second@meta.data[which(PCa_sce_anno_second@meta.data$gleason %in% c("3+4","4+3","3+3")),]$G_group <- "LG"
PCa_sce_anno_second@meta.data[which(PCa_sce_anno_second@meta.data$gleason == "4+5"),]$G_group <- "HG"

data <- as.data.frame(cbind(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings,
                            as.character(PCa_sce_anno_second@meta.data$cell_type) ))
colnames(data) <- c("umap_x","umap_y","celltype") 
data$umap_x <- as.numeric(data$umap_x)
data$umap_y <- as.numeric(data$umap_y)

library(ks)
#Epithelium
Epithelium <- data[which(data$celltype == "Epithelium"),c(1,2)]
colnames(Epithelium) <- c("umap_x","umap_y")
Epithelium$umap_x <- as.numeric(Epithelium$umap_x)
Epithelium$umap_y <- as.numeric(Epithelium$umap_y)
Epithelium <- kde(x=Epithelium)
Epithelium <- with(Epithelium, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#T_cells
T_cells <- data[which(data$celltype == "T_cells"),c(1,2)]
colnames(T_cells) <- c("umap_x","umap_y")
T_cells$umap_x <- as.numeric(T_cells$umap_x)
T_cells$umap_y <- as.numeric(T_cells$umap_y)
T_cells <- kde(x=T_cells)
T_cells <- with(T_cells, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#B_cells
B_cells <- data[which(data$celltype == "B_cells"),c(1,2)]
colnames(B_cells) <- c("umap_x","umap_y")
B_cells$umap_x <- as.numeric(B_cells$umap_x)
B_cells$umap_y <- as.numeric(B_cells$umap_y)
B_cells <- kde(x=B_cells)
B_cells <- with(B_cells, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#Myeloid
Myeloid <- data[which(data$celltype == "Myeloid"),c(1,2)]
colnames(Myeloid) <- c("umap_x","umap_y")
Myeloid$umap_x <- as.numeric(Myeloid$umap_x)
Myeloid$umap_y <- as.numeric(Myeloid$umap_y)
Myeloid <- kde(x=Myeloid)
Myeloid <- with(Myeloid, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#Endothelium
Endothelium <- data[which(data$celltype == "Endothelium"),c(1,2)]
colnames(Endothelium) <- c("umap_x","umap_y")
Endothelium$umap_x <- as.numeric(Endothelium$umap_x)
Endothelium$umap_y <- as.numeric(Endothelium$umap_y)
Endothelium <- kde(x=Endothelium)
Endothelium <- with(Endothelium, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#Fibroblast
Fibroblast <- data[which(data$celltype == "Fibroblast"),c(1,2)]
colnames(Fibroblast) <- c("umap_x","umap_y")
Fibroblast$umap_x <- as.numeric(Fibroblast$umap_x)
Fibroblast$umap_y <- as.numeric(Fibroblast$umap_y)
Fibroblast <- kde(x=Fibroblast)
Fibroblast <- with(Fibroblast, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
#Mast_cells
Mast_cells <- data[which(data$celltype == "Mast_cells"),c(1,2)]
colnames(Mast_cells) <- c("umap_x","umap_y")
Mast_cells$umap_x <- as.numeric(Mast_cells$umap_x)
Mast_cells$umap_y <- as.numeric(Mast_cells$umap_y)
Mast_cells <- kde(x=Mast_cells)
Mast_cells <- with(Mast_cells, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])

#raster polygon
UMAP_EMB <- subset(PCa_sce_anno_second,subset = G_group == "LG")
UMAP_EMB <- UMAP_EMB@reductions[["umap"]]@cell.embeddings
UMAP_EMB <- as.data.frame(UMAP_EMB)
xlims=c(min(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])-1,max(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])+1)
ylims=c(min(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])-1,max(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])+1)
library(viridis)
LG_plot=ggplot(UMAP_EMB, aes(x=UMAP_1, y=UMAP_2) ) + xlim(xlims) + ylim(ylims)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
  #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
  
  scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
  geom_point(aes(x=UMAP_1, y=UMAP_2), col='#FCFDBFFF',size=0.00001,alpha=0.5)+
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 

LG_plot=LG_plot+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

LG_plot <- LG_plot + theme(axis.title.x=element_blank(),
                           # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_line(),
                           axis.title.y=element_blank(),
                           axis.ticks.y=element_blank(),
                           # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                           axis.text.y=element_blank(),
                           axis.ticks.length = unit(0, "cm")
)
con.color <- 'gray40'
con.size=1
cl='white'
LG_plot=LG_plot+geom_path(aes(x, y), data=data.frame(Epithelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(T_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(B_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myeloid),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Endothelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Fibroblast),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Mast_cells),col=cl,linetype = 2,size=con.size)



UMAP_EMB <- subset(PCa_sce_anno_second,subset = G_group == "HG")
UMAP_EMB <- UMAP_EMB@reductions[["umap"]]@cell.embeddings
UMAP_EMB <- as.data.frame(UMAP_EMB)
xlims=c(min(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])-1,max(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])+1)
ylims=c(min(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])-1,max(PCa_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])+1)
library(viridis)
HG_plot=ggplot(UMAP_EMB, aes(x=UMAP_1, y=UMAP_2) ) + xlim(xlims) + ylim(ylims)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
  #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
  
  scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
  geom_point(aes(x=UMAP_1, y=UMAP_2), col='#FCFDBFFF',size=0.00001,alpha=0.5)+
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 

HG_plot=HG_plot+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

HG_plot <- HG_plot + theme(axis.title.x=element_blank(),
                           # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_line(),
                           axis.title.y=element_blank(),
                           axis.ticks.y=element_blank(),
                           # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                           axis.text.y=element_blank(),
                           axis.ticks.length = unit(0, "cm")
)
con.color <- 'gray40'
con.size=1
cl='white'
HG_plot=HG_plot+geom_path(aes(x, y), data=data.frame(Epithelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(T_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(B_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myeloid),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Endothelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Fibroblast),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Mast_cells),col=cl,linetype = 2,size=con.size)

LG_plot+HG_plot
ggsave("PCa_Destiny.png",height = 6,width = 12)
w() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 

T1a_plot=T1a_plot+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

T1a_plot <- T1a_plot + theme(axis.title.x=element_blank(),
                             # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_line(),
                             axis.title.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.y=element_blank(),
                             axis.ticks.length = unit(0, "cm")
)
con.color <- 'gray40'
con.size=1
cl='white'
T1a_plot=T1a_plot+geom_path(aes(x, y), data=data.frame(Epithelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(T_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(B_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myeloid),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Endothelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Fibroblast),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Mast_cells),col=cl,linetype = 2,size=con.size)


UMAP_EMB <- subset(ccRCC_sce_anno_second,subset = TNM == "T1bN0M0")
UMAP_EMB <- UMAP_EMB@reductions[["umap"]]@cell.embeddings
UMAP_EMB <- as.data.frame(UMAP_EMB)
xlims=c(min(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])-1,max(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])+1)
ylims=c(min(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])-1,max(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])+1)
library(viridis)
T1b_plot=ggplot(UMAP_EMB, aes(x=UMAP_1, y=UMAP_2) ) + xlim(xlims) + ylim(ylims)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
  #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
  
  scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
  geom_point(aes(x=UMAP_1, y=UMAP_2), col='#FCFDBFFF',size=0.00001,alpha=0.5)+
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 

T1b_plot=T1b_plot+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

T1b_plot <- T1b_plot + theme(axis.title.x=element_blank(),
                             # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_line(),
                             axis.title.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.y=element_blank(),
                             axis.ticks.length = unit(0, "cm")
)
con.color <- 'gray40'
con.size=1
cl='white'
T1b_plot=T1b_plot+geom_path(aes(x, y), data=data.frame(Epithelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(T_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(B_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myeloid),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Endothelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Fibroblast),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Mast_cells),col=cl,linetype = 2,size=con.size)


UMAP_EMB <- subset(ccRCC_sce_anno_second,subset = TNM == "T2aN0M0")
UMAP_EMB <- UMAP_EMB@reductions[["umap"]]@cell.embeddings
UMAP_EMB <- as.data.frame(UMAP_EMB)
xlims=c(min(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])-1,max(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,1])+1)
ylims=c(min(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])-1,max(ccRCC_sce_anno_second@reductions[["umap"]]@cell.embeddings[,2])+1)
library(viridis)
T2a_plot=ggplot(UMAP_EMB, aes(x=UMAP_1, y=UMAP_2) ) + xlim(xlims) + ylim(ylims)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
  #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
  
  scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
  geom_point(aes(x=UMAP_1, y=UMAP_2), col='#FCFDBFFF',size=0.00001,alpha=0.5)+
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 

T2a_plot=T2a_plot+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))

T2a_plot <- T2a_plot + theme(axis.title.x=element_blank(),
                             # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_line(),
                             axis.title.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                             axis.text.y=element_blank(),
                             axis.ticks.length = unit(0, "cm")
)

con.color <- 'gray40'
con.size=1
cl='white'
T2a_plot=T2a_plot+geom_path(aes(x, y), data=data.frame(Epithelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(T_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(B_cells),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Myeloid),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Endothelium),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Fibroblast),col=cl,linetype = 2,size=con.size)+
  geom_path(aes(x, y), data=data.frame(Mast_cells),col=cl,linetype = 2,size=con.size)

library(patchwork)
(WHOI_plot+T1a_plot)/(WHOII_plot+T1b_plot)/(WHOIII_plot+T2a_plot)
ggsave("ccRCC_Destiny.png",height = 16,width = 12)

#Spatial data ------
library(Seurat)
library(stringr)
library(ggplot2)
sce_PCa_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")

options(future.globals.maxSize=4000000000000000000)
PCa_Spatial <- list()
sce_PCa_tumor <- SCTransform(sce_PCa_tumor, assay = "RNA", return.only.var.genes = T, verbose = FALSE)
saveRDS(sce_PCa_tumor, "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\sp\\PCa_sc_ref.rds")

all_count <- readRDS("D:\\Urinary system tumors\\work\\0_data\\3_PCa\\3_PCa_sp_GSE181294\\slide.seq.raw.counts.rds")

all_ano <- readRDS("D:\\Urinary system tumors\\work\\0_data\\3_PCa\\3_PCa_sp_GSE181294\\slide.seq.ano.rds")
all_ano$sample <- str_extract(row.names(all_ano), ".+?_")
all_ano$sample <- gsub("[_]","",all_ano$sample)
all_ano$barcode <- row.names(all_ano)

all_bead <- readRDS("D:\\Urinary system tumors\\work\\0_data\\3_PCa\\3_PCa_sp_GSE181294\\slide.seq.BeadLocations.rds")

i <- "Tumor01"
bead <- all_bead[[i]]
count <- all_count[[i]]
slide.seq = CreateSeuratObject(count, assay="Spatial")
slide.seq@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = bead
)
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
ElbowPlot(slide.seq ,ndims = 50)
slide.seq <- RunUMAP(slide.seq, dims = 1:20)
slide.seq <- FindNeighbors(slide.seq, dims = 1:20)
slide.seq <- FindClusters(slide.seq, resolution = 0.9, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE) 
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2


anchors <- FindTransferAnchors(reference = sce_PCa_tumor, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = sce_PCa_tumor$cell_type, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Epithelium", "Endothelium"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("Fibroblast", "Myeloid"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("T-cells"), alpha = c(0.1, 1))

PCa_Spatial[[i]] <- slide.seq



i <- "Tumor02"
bead <- all_bead[[i]]
count <- all_count[[i]]
slide.seq = CreateSeuratObject(count, assay="Spatial")
slide.seq@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = bead
)
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
ElbowPlot(slide.seq ,ndims = 50)
slide.seq <- RunUMAP(slide.seq, dims = 1:20)
slide.seq <- FindNeighbors(slide.seq, dims = 1:20)
slide.seq <- FindClusters(slide.seq, resolution = 0.9, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE) 
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2


anchors <- FindTransferAnchors(reference = sce_PCa_tumor, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = sce_PCa_tumor$cell_type, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Epithelium", "Endothelium"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("Fibroblast", "Myeloid"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("T-cells"), alpha = c(0.1, 1))

PCa_Spatial[[i]] <- slide.seq

i <- "Tumor07"
bead <- all_bead[[i]]
count <- all_count[[i]]
slide.seq = CreateSeuratObject(count, assay="Spatial")
slide.seq@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = bead
)
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
ElbowPlot(slide.seq ,ndims = 50)
slide.seq <- RunUMAP(slide.seq, dims = 1:20)
slide.seq <- FindNeighbors(slide.seq, dims = 1:20)
slide.seq <- FindClusters(slide.seq, resolution = 0.9, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE) 
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2


anchors <- FindTransferAnchors(reference = sce_PCa_tumor, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = sce_PCa_tumor$cell_type, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Epithelium", "Endothelium"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("Fibroblast", "Myeloid"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("T-cells"), alpha = c(0.1, 1))

PCa_Spatial[[i]] <- slide.seq

i <- "Tumor08"
bead <- all_bead[[i]]
count <- all_count[[i]]
slide.seq = CreateSeuratObject(count, assay="Spatial")
slide.seq@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = bead
)
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
ElbowPlot(slide.seq ,ndims = 50)
slide.seq <- RunUMAP(slide.seq, dims = 1:20)
slide.seq <- FindNeighbors(slide.seq, dims = 1:20)
slide.seq <- FindClusters(slide.seq, resolution = 0.9, verbose = FALSE)
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE) 
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2


anchors <- FindTransferAnchors(reference = sce_PCa_tumor, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = sce_PCa_tumor$cell_type, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Epithelium", "Endothelium"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("Fibroblast", "Myeloid"), alpha = c(0.1, 1))
SpatialFeaturePlot(slide.seq, features = c("T-cells"), alpha = c(0.1, 1))

PCa_Spatial[[i]] <- slide.seq
saveRDS(PCa_Spatial, "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\sp\\PCa_Spatial.rds")


setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\LG_3")
LG_3 <- Read10X_h5("LG_3_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'LG_3',assay = "Spatial")
coord.df = read.csv("./pos/LG_3_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
LG_3@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
LG_3 <- SCTransform(LG_3, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
LG_3 <- RunPCA(LG_3, assay = "SCT", verbose = FALSE)
LG_3 <- FindNeighbors(LG_3, reduction = "pca", dims = 1:30)
LG_3 <- FindClusters(LG_3, verbose = FALSE)
LG_3 <- RunUMAP(LG_3, reduction = "pca", dims = 1:30)

SpatialDimPlot(LG_3,pt.size.factor = 4)

#HG1.3
setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\HG_1.3")
HG_1.3 <- Read10X_h5("HG_1.3_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'HG_1.3',assay = "Spatial")


coord.df = read.csv("./pos/HG_1.3_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
HG_1.3@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
HG_1.3 <- SCTransform(HG_1.3, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
HG_1.3 <- RunPCA(HG_1.3, assay = "SCT", verbose = FALSE)
HG_1.3 <- FindNeighbors(HG_1.3, reduction = "pca", dims = 1:30)
HG_1.3 <- FindClusters(HG_1.3, verbose = FALSE)
HG_1.3 <- RunUMAP(HG_1.3, reduction = "pca", dims = 1:30)

SpatialDimPlot(HG_1.3,pt.size.factor = 4)




#LG_2
setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\LG_2")
LG_2 <- Read10X_h5("LG_2_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'LG_2',assay = "Spatial")


coord.df = read.csv("./pos/LG_2_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
LG_2@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
LG_2 <- SCTransform(LG_2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
LG_2 <- RunPCA(LG_2, assay = "SCT", verbose = FALSE)
LG_2 <- FindNeighbors(LG_2, reduction = "pca", dims = 1:30)
LG_2 <- FindClusters(LG_2, verbose = FALSE)
LG_2 <- RunUMAP(LG_2, reduction = "pca", dims = 1:30)

SpatialDimPlot(LG_2,pt.size.factor = 4)


#HG_1.2
setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\HG_1.2")
HG_1.2 <- Read10X_h5("HG_1.2_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'HG_1.2',assay = "Spatial")


coord.df = read.csv("./pos/HG_1.2_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
HG_1.2@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
HG_1.2 <- SCTransform(HG_1.2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
HG_1.2 <- RunPCA(HG_1.2, assay = "SCT", verbose = FALSE)
HG_1.2 <- FindNeighbors(HG_1.2, reduction = "pca", dims = 1:30)
HG_1.2 <- FindClusters(HG_1.2, verbose = FALSE)
HG_1.2 <- RunUMAP(HG_1.2, reduction = "pca", dims = 1:30)

SpatialDimPlot(HG_1.2,pt.size.factor = 4)



#HG_2
setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\HG_2")
HG_2 <- Read10X_h5("HG_2_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'HG_2',assay = "Spatial")


coord.df = read.csv("./pos/HG_2_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
HG_2@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
HG_2 <- SCTransform(HG_2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
HG_2 <- RunPCA(HG_2, assay = "SCT", verbose = FALSE)
HG_2 <- FindNeighbors(HG_2, reduction = "pca", dims = 1:30)
HG_2 <- FindClusters(HG_2, verbose = FALSE)
HG_2 <- RunUMAP(HG_2, reduction = "pca", dims = 1:30)

SpatialDimPlot(HG_2,pt.size.factor = 4)

#HG_3.2
setwd("D:\\Urinary system tumors\\work\\0_data\\1_ccRCC\\3_ccRCC_sp\\E-MTAB-12767\\HG_3.2")
HG_3.2 <- Read10X_h5("HG_3.2_filtered_feature_bc_matrix.h5")%>% CreateSeuratObject(project = 'HG_3.2',assay = "Spatial")


coord.df = read.csv("./pos/HG_3.2_tissue_positions_list.csv",header = T,row.names = 1)
coord.df <- coord.df[which(coord.df$X0 != 0),]
coord.df <- coord.df[,c(4:5)]
HG_3.2@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)
HG_3.2 <- SCTransform(HG_3.2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
HG_3.2 <- RunPCA(HG_3.2, assay = "SCT", verbose = FALSE)
HG_3.2 <- FindNeighbors(HG_3.2, reduction = "pca", dims = 1:30)
HG_3.2 <- FindClusters(HG_3.2, verbose = FALSE)
HG_3.2 <- RunUMAP(HG_3.2, reduction = "pca", dims = 1:30)

SpatialDimPlot(HG_3.2,pt.size.factor = 4)
ccRCC <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/1_ccRCC/1_ccRCC_tumor_sce/4_ccRCC_tumor_sce_after_cell_type.rds")
set.seed(123)
selected_cells <- ccRCC@meta.data %>%
  group_by(cell_type)%>%
  sample_frac(0.1)%>%
  pull(cell_names)
ccRCC_sc_subset <- subset(ccRCC,cells=selected_cells)
ccRCC_sc_subset <- SCTransform(ccRCC_sc_subset, assay = "RNA", return.only.var.genes = T, verbose = FALSE)
saveRDS(ccRCC_sc_subset, "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\sp\\ccRCC_sc_subset.rds")
ccRCC_sc_subset2 <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\sp\\ccRCC_sc_subset.rds")
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = LG_3, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = LG_3[["pca"]], dims = 1:30)
LG_3[["predictions"]] <- predictions.assay

DefaultAssay(LG_3) <- "predictions"
#LG
SpatialFeaturePlot(LG_3, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(LG_3, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(LG_3, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)



#HG
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = HG_1.3, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = HG_1.3[["pca"]], dims = 1:30)
HG_1.3[["predictions"]] <- predictions.assay

DefaultAssay(HG_1.3) <- "predictions"
SpatialFeaturePlot(HG_1.3, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_1.3, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_1.3, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)

##LG_2
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = LG_2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = LG_2[["pca"]], dims = 1:30)
LG_2[["predictions"]] <- predictions.assay

DefaultAssay(LG_2) <- "predictions"
SpatialFeaturePlot(LG_2, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(LG_2, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(LG_2, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)


##HG_2
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = HG_2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = HG_2[["pca"]], dims = 1:30)
HG_2[["predictions"]] <- predictions.assay

DefaultAssay(HG_2) <- "predictions"
SpatialFeaturePlot(HG_2, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_2, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_2, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)

##HG_3.2
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = HG_3.2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = HG_3.2[["pca"]], dims = 1:30)
HG_3.2[["predictions"]] <- predictions.assay

DefaultAssay(HG_3.2) <- "predictions"
SpatialFeaturePlot(HG_3.2, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_3.2, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_3.2, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)


##HG_1.2
anchors <- FindTransferAnchors(reference = ccRCC_sc_subset2, query = HG_1.2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = ccRCC_sc_subset2$cell_type, prediction.assay = TRUE,
                                  weight.reduction = HG_1.2[["pca"]], dims = 1:30)
HG_1.2[["predictions"]] <- predictions.assay

DefaultAssay(HG_1.2) <- "predictions"
SpatialFeaturePlot(HG_1.2, features = c("Epithelium", "Endothelium"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_1.2, features = c("Fibroblast", "Myeloid"), pt.size.factor = 4, ncol = 2, crop = TRUE)
SpatialFeaturePlot(HG_1.2, features = c("T-cells"), pt.size.factor = 4, ncol = 2, crop = TRUE)


ccRCC_Spatial <- list()
ccRCC_Spatial[["LG_3"]] <- LG_3
ccRCC_Spatial[["HG_1.2"]] <- HG_1.2
ccRCC_Spatial[["HG_1.3"]] <- HG_1.3
ccRCC_Spatial[["HG_3.2"]] <- HG_3.2
ccRCC_Spatial[["HG_2"]] <- HG_2
saveRDS(ccRCC_Spatial, "D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\sp\\ccRCC_Spatial.rds")

#mistyR---------
# MISTy
library(mistyR)
library(future)
# data manipulation
library(dplyr)
library(purrr)
library(distances)
# plotting
library(ggplot2)
library(Seurat)
setwd("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR")
data("synthetic")
ccRCC_Spatial <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/ccRCC_Spatial.rds")
PCa_Spatial <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/PCa_Spatial.rds")

ccRCC_misty <- list()
for(i in names(ccRCC_Spatial)){
  x <- ccRCC_Spatial[[i]]@images[["image"]]@coordinates
  colnames(x) <- c("row", "col")
  y <- t(ccRCC_Spatial[[i]]@assays[["predictions"]]@data)
  y <- y[row.names(x),]
  colnames(y) <- gsub("\\-", "_", colnames(y))
  x <- cbind(x,y)
  colSums(x)
  x <- x[,-which(colSums(x)==0)]
  x <- x[,-which(colnames(x) %in% c("max"))]
  ccRCC_misty[[i]] <- x
}



names(ccRCC_misty) <- paste("ccRCC",names(ccRCC_misty),sep = "_")

PCa_misty <- list()
for(i in names(PCa_Spatial)){
  x <- PCa_Spatial[[i]]@images[["image"]]@coordinates
  colnames(x) <- c("row", "col")
  y <- t(PCa_Spatial[[i]]@assays[["predictions"]]@data)
  y <- y[row.names(x),]
  colnames(y) <- gsub("\\-", "_", colnames(y))
  x <- cbind(x,y)
  if(0 %in% colSums(x)){
    x <- x[,-which(colSums(x)==0)]}
  x <- x[,-which(colnames(x) %in% c("max"))]
  PCa_misty[[i]] <- x
}
names(PCa_misty) <- paste("PCa",names(PCa_misty),sep = "_")
for(i in 1:5){
  geom_dist <- as.matrix(distances(ccRCC_misty[[i]][,c(1,2)]))
  dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
  paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))
  paraview_radius
  result.folders <- ccRCC_misty[i] %>% imap_chr(function(sample, name) {
    sample.expr <- sample %>% select(-c(row, col))
    sample.pos <- sample %>% select(row, col)
    create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = paraview_radius) %>%
      run_misty(results.folder = paste0("ccRCCresults", .Platform$file.sep, name),seed = 1)
  }) 
}
result.folders <- c("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/ccRCCresults/ccRCC_HG_1.2",
                    "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/ccRCCresults/ccRCC_HG_1.3",
                    "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/ccRCCresults/ccRCC_HG_2",
                    "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/ccRCCresults/ccRCC_HG_3.2",
                    "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/ccRCCresults/ccRCC_LG_3")
names(result.folders) <- c("ccRCC_HG_1.2",
                           "ccRCC_HG_1.3",
                           "ccRCC_HG_2",
                           "ccRCC_HG_3.2",
                           "ccRCC_LG_3")
for(i in 1:5){
  misty.results <- collect_results(result.folders[i])
  misty.results %>%
    plot_improvement_stats("gain.R2") %>%
    plot_improvement_stats("gain.RMSE")
  
  misty.results$improvements %>%
    filter(measure == "p.R2") %>%
    group_by(target) %>%
    summarize(mean.p = mean(value)) %>%arrange(mean.p)
  
  
  misty.results %>% plot_view_contributions()
  
  
  misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.5)
  ggsave(paste0(names(result.folders[i]),"_intra.pdf"),height = 4,width = 4)
  #misty.results %>% plot_interaction_heatmap(view = "para.25", cutoff = 0.5)
  #ggsave(paste0(names(result.folders[i]),"_para.25.pdf"),height = 4,width = 4)
  #misty.results %>% plot_contrast_heatmap("intra", "para.25", cutoff = 0.5)
  
  pdf(paste0(names(result.folders[i]),"_intra_interaction.pdf"),height = 4,width = 4)
  misty.results %>% plot_interaction_communities("intra")
  dev.off()
  #pdf(paste0(names(result.folders[i]),"_para.25_interaction.pdf"),height = 4,width = 4)
  #misty.results %>% plot_interaction_communities("para.25", cutoff = 0.5)
  #dev.off()
}

for(i in 1:5){
  result.folders <- ccRCC_misty %>% imap_chr(function(sample, name) {
    sample.expr <- sample %>% select(-c(row, col))
    sample.pos <- sample %>% select(row, col)
    create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 25) %>%
      run_misty(results.folder = paste0("ccRCCresults_25", .Platform$file.sep, name),seed = 1)
  }) 
}

ccRCC_25_result.folders <- result.folders
misty.results <- collect_results(ccRCC_25_result.folders[2:5])
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>%
  summarize(mean.p = mean(value)) %>%arrange(mean.p)


misty.results %>% plot_view_contributions()


misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.5)
misty.results %>% plot_interaction_heatmap(view = "para.25", cutoff = 0.5)
misty.results %>% plot_contrast_heatmap("intra", "para.25", cutoff = 0.5)


misty.results %>% plot_interaction_communities("intra")
ggsave("ccRCC_HG_intra.pdf",height = 4,width = 4)


misty.results <- collect_results(ccRCC_25_result.folders[1])
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>%
  summarize(mean.p = mean(value)) %>%arrange(mean.p)

misty.results %>% plot_view_contributions()

misty.results %>% plot_interaction_communities("intra")
ggsave("ccRCC_LG_intra.pdf",height = 4,width = 4)

for(i in 1:4){
  geom_dist <- as.matrix(distances(PCa_misty[[i]][,c(1,2)]))
  dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
  paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))
  paraview_radius
  PCa_result.folders <- PCa_misty[i] %>% imap_chr(function(sample, name) {
    sample.expr <- sample %>% select(-c(row, col))
    sample.pos <- sample %>% select(row, col)
    create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = paraview_radius) %>%
      run_misty(results.folder = paste0("PCa_results", .Platform$file.sep, name),seed = 1)
  })
}

for(i in 1:4){
  PCa_result.folders <- PCa_misty %>% imap_chr(function(sample, name) {
    sample.expr <- sample %>% select(-c(row, col))
    sample.pos <- sample %>% select(row, col)
    create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 20) %>%
      run_misty(results.folder = paste0("PCa_results_20", .Platform$file.sep, name),seed = 1)
  })
}
PCa_20_result.folders <- PCa_result.folders

PCa_result.folders <- c("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/PCa_results_20/PCa_Tumor01",
                        "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/PCa_results_20/PCa_Tumor02",
                        "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/PCa_results_20/PCa_Tumor07",
                        "D:/Urinary system tumors/work/1_scRNA_landsacape/sp/mistyR/PCa_results_20/PCa_Tumor08" )
names(PCa_result.folders) <- c("PCa_Tumor01",
                               "PCa_Tumor02",
                               "PCa_Tumor07",
                               "PCa_Tumor08")
for(i in 1:4){
  misty.results <- collect_results(PCa_result.folders[i])
  misty.results %>%
    plot_improvement_stats("gain.R2") %>%
    plot_improvement_stats("gain.RMSE")
  
  misty.results$improvements %>%
    filter(measure == "p.R2") %>%
    group_by(target) %>%
    summarize(mean.p = mean(value)) %>%arrange(mean.p)
  
  misty.results %>% plot_view_contributions()
  
  
  misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.5)
  ggsave(paste0(names(PCa_result.folders[i]),"_intra.pdf"),height = 4,width = 4)
  #misty.results %>% plot_interaction_heatmap(view = "para", cutoff = 0.5)
  #ggsave(paste0(names(PCa_result.folders[i]),"_para.pdf"),height = 4,width = 4)
  #misty.results %>% plot_contrast_heatmap("intra", "para.10", cutoff = 0.5)
  
  pdf(paste0(names(PCa_result.folders[i]),"_intra_interaction.pdf"),height = 4,width = 4)
  misty.results %>% plot_interaction_communities("intra")
  dev.off()
  #pdf(paste0(names(PCa_result.folders[i]),"_para_interaction.pdf"),height = 4,width = 4)
  #misty.results %>% plot_interaction_communities("para", cutoff = 0.5)
  #dev.off()
}

#LG
misty.results <- collect_results(PCa_20_result.folders[3:4])
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>%
  summarize(mean.p = mean(value)) %>%arrange(mean.p)


misty.results %>% plot_view_contributions()


misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.5)
misty.results %>% plot_interaction_heatmap(view = "para.20", cutoff = 0.5)
misty.results %>% plot_contrast_heatmap("intra", "para.20", cutoff = 0.5)

"PCa_LG_intra.pdf"
misty.results %>% plot_interaction_communities("intra")



#HG
misty.results <- collect_results(PCa_20_result.folders[1:2])
misty.results %>%
  plot_improvement_stats("gain.R2") %>%
  plot_improvement_stats("gain.RMSE")

misty.results$improvements %>%
  filter(measure == "p.R2") %>%
  group_by(target) %>%
  summarize(mean.p = mean(value)) %>%arrange(mean.p)


misty.results %>% plot_view_contributions()


misty.results %>% plot_interaction_heatmap(view = "intra", cutoff = 0.5)
misty.results %>% plot_interaction_heatmap(view = "para.20", cutoff = 0.5)
misty.results %>% plot_contrast_heatmap("intra", "para.20", cutoff = 0.5)


misty.results %>% plot_interaction_communities("intra")
"PCa_HG_intra.pdf"

#scATAC analysis ------
##ccRCC-----
# step.1 
# 

subdirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)

# 
sampleNames <- sapply(subdirs, function(subdir) {
  fragments <- list.files(subdir, pattern = "fragments.tsv.gz$", full.names = TRUE)
  if (length(fragments) > 0) {
    return(basename(subdir)) 
  } else {
    return(NA)
  }
})


sampleNames <- sampleNames[!is.na(sampleNames)]

inputFiles <- file.path(subdirs, "fragments.tsv.gz")
# inputFiles <- inputFiles[-c(1,2,3,23)]

inputFiles
names(inputFiles) <- sampleNames  


print(inputFiles)
print(sampleNames)

ArrowFiles <- createArrowFiles(
  inputFiles      = inputFiles,
  sampleNames     = sampleNames, 
  minTSS          = 8,           
  minFrags        = 1000,      
  addTileMat      = TRUE,        
  addGeneScoreMat = TRUE,
  subThreading    = FALSE)        


proj <- ArchRProject(
  ArrowFiles       = ArrowFiles,
  outputDirectory  = outputDirectory,
  copyArrows       = TRUE,
  geneAnnotation   = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo         = TRUE,
  threads          = getArchRThreads())
proj

proj <- addIterativeLSI(ArchRProj        = proj, 
                        useMatrix        = "TileMatrix", 
                        name             = "IterativeLSI", 
                        outlierQuantiles = NULL)
proj <- addClusters(proj,reducedDims     = "IterativeLSI")
proj <- addUMAP(proj,reducedDims         = "IterativeLSI")
proj <- addImputeWeights(proj)

proj <- addHarmony(
  ArchRProj   = proj,
  reducedDims = "IterativeLSI",
  name        = "Harmony",
  groupBy     = "Sample",
  force       = TRUE)

proj <- addClusters(
  input       = proj,
  reducedDims = "Harmony",
  method      = "Seurat",
  name        = "Clusters",
  resolution  = 0.6,
  maxClusters = 12,
  force       = TRUE)

proj <- addUMAP(  ArchRProj = proj,reducedDims = "Harmony",name = "UMAPHarmony",nNeighbors = 30,minDist = 0.5,metric = "cosine",force = TRUE)

library(Seurat)
scRNA    <- readRDS(   'F:/Mingjie/Project_Urinary/01.Data/02.scRNAData/ccRCC_scRNA.rds')
metadata <- read.table('F:/Mingjie/Project_Urinary/01.Data/02.scRNAData/ccRCC_scRNA_anno_second.txt',sep = '\t',header = FALSE, fill = NA, row.names = 1)
colnames <- metadata[1,]
colnames <- c('orig.ident',colnames)
colnames <- colnames[-28]
colnames(metadata) <- colnames
metadata <- metadata[-1,]
scRNA <- AddMetaData(object = scRNA,metadata = metadata,col.name = c('cell_type_second','filter'))

Idents(scRNA) <- scRNA$cell_type
table(Idents(scRNA))

proj <- addGeneIntegrationMatrix(
  ArchRProj   = proj, 
  useMatrix   = "GeneScoreMatrix",
  matrixName  = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA       = scRNA,
  addToArrow  = TRUE,
  groupRNA    = "cell_type",
  nameCell    = "predictedCell_Un",
  nameGroup   = "predictedGroup_Un",
  nameScore   = "predictedScore_Un",
  force       = TRUE
)

### Peak calling
proj <- addGroupCoverages(ArchRProj        = proj, groupBy = "predictedGroup_Un")
proj <- addReproduciblePeakSet(ArchRProj   = proj, groupBy = "predictedGroup_Un", peakMethod = 'Tiles',method = 'p')
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

markersPeaks <- getMarkerFeatures(
  ArchRProj  = proj, 
  useMatrix  = "PeakMatrix", 
  groupBy    = "predictedGroup_Un",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")
saveRDS(markersPeaks,"ccRCC_markersPeaks.rds")
markerList_GR <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet  = "cisbp", name = "Motif")

### peakAnnoEnrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker       = markersPeaks,
  ArchRProj      = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# Peak2GeneLinkages
getAvailableMatrices(proj)
proj <- addCoAccessibility(ArchRProj = proj,reducedDims = "Harmony")
proj <- addPeak2GeneLinks(ArchRProj  = proj,reducedDims = "Harmony")

p1  <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample",            embedding = "UMAPHarmony")
p2  <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters",          embedding = "UMAPHarmony")
p3  <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAPHarmony", pal = pal)
p4  <- plotMarkerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5",transpose = TRUE)
p5  <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE) # ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")  
p6  <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "predictedGroup_Un")

plotPDF(p1,p2,p3, name = "1.UMAP-Harmony.pdf"                    , ArchRProj = proj, addDOC = FALSE, width = 6 , height = 6)
plotPDF(p4,       name = "2.Peak-Marker-Heatmap.pdf"             , ArchRProj = proj, addDOC = FALSE, width = 12, height = 6)
plotPDF(p5,       name = "3.Motifs-Enriched-Marker-Heatmap_5.pdf", ArchRProj = proj, addDOC = FALSE, width = 6 , height = 6)
plotPDF(p6,       name = "4.Peak2GeneHeatmap.pdf"                , ArchRProj = proj, addDOC = FALSE, width = 6 , height = 6)

saveArchRProject(ArchRProj = proj, outputDirectory = outputDirectory)




##PCa -------
#obtain fragment
stty susp undef

（1）prefetch --option-file atac_ara_test.txt  
atac_ara_test.txt
SRR14151148
SRR14151149
（2）
mkdir atac_data

less atac_ara_test.txt |while  read id;
do 
cd /home/lmh/atac/atac_test/$id;
mv * /home/lmh/atac/atac_test/atac_data;
done

ls atac_data/
  fasterq-dump --split-files SRR14151148.sra --include-technical
gzip SRR14151148_1fastq


while read id; do
id_without_extension=${id%.sra}

fasterq-dump --split-files $id --include-technical

gzip ${id_without_extension}_1.fastq
gzip ${id_without_extension}_2.fastq

done < atac_ara_test.txt



cat atac_ara_test.txt | parallel -j 4 'fasterq-dump --split-files /home/lmh/atac/atac_test/atac_data/{} --include-technical'


parallel -j 4 'gzip *.fastq'

parallel -j 4 "    
    fastq-dump --split-3 --gzip {1}    
" ::: $(ls *.sra)     

#2、fasq preprocessing


paste <(zcat SRR14151148_3.fastq.gz | awk 'NR % 4 == 2 {print "CB:Z:" $0}')       <(zcat SRR14151148_1.fastq.gz | awk 'NR % 4 == 1 {print $0}')       > merged.fastq


fastqc -t 4 -o /home/lmh/atac/test_2 SRR14151148.fastq.gz
cd /home/lmh/atac/test_2/fastq


awk 'NR % 4 == 2 {print substr($0, 1, 18)}' SRR14151148.fastq.gz > barcodes.txt

less /home/lmh/atac/atac_test/SRA/atac_ara_test.txt |while  read id;
do 
fastp -i ${id}_1.fastq.gz -o ../${id}.fastq.gz  -f 16 -t 2 -L
done

fastqc -t 4 -o /home/lmh/atac/atac_test/atac_data/fastp/fastp_fastqc/ /home/lmh/atac/atac_test/atac_data/fastp/*.gz
cd /home/lmh/atac/atac_test/atac_data/fastp/
  
  multiqc .

trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o /mnt/d/ATAC/trim/  $fq1 $fq2 &
  
  trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --cores 4 -o /home/lmh/atac/test_2 SRR14151148.fastq.gz &
  
  
  
  source venv/bin/activate


$ which bwa
/usr/bin/bwa 
$ snaptools index-genome  \
--input-fasta=gencode.v47.transcripts.fa  \
--output-prefix=hg38  \
--aligner=bwa  \
--path-to-aligner=/usr/bin/  \
--num-threads=4

bowtie2-build -f /home/lmh/atac/test_2/GRCh38.fa --threads 4 GRCh38


snaptools align-single-end  \   
--input-reference=gencode.v47.transcripts.fa  \
--input-fastq=SRR14151148.fastq.gz  \
--output-bam=SRR14151148.bam  \
--aligner=bwa  \
--path-to-aligner=/usr/bin/  \
--read-fastq-command=zcat  \
--min-cov=0  \
--num-threads=4  \
--if-sort=True  \
--tmp-folder=./  \
--overwrite=TRUE 

bowtie2_index
hg38.1.bt2
hg38.2.bt2
hg38.3.bt2
hg38.4.bt2
hg38.rev.1.bt2
hg38.rev.2.bt2


mkdir -p /home/lmh/atac/test_2/alignment
bowtie2_index=/home/lmh/atac/test_2/GRCh38
align_dir=/home/lmh/atac/atac_test/SRA/alignment


bowtie2  -p 6 -x  $bowtie2_index --very-sensitive -X 2000 -U  /home/lmh/atac/test_2/alignment/SRR14151148_fastp_fastqc.zip \
2>/home/lmh/atac/test_2/alignment/SRR14151148.summary \
-S /home/lmh/atac/test_2/alignment/SRR14151148.sam 

while read id; do
bowtie2 -p 7 -x "$bowtie2_index" --very-sensitive -X 2000 -U "${id}.fastq.gz" \
2>"$align_dir/${id}.summary" \
-S "$align_dir/${id}.sam"
done < /home/lmh/atac/atac_test/SRA/atac_ara_test.txt



cd /home/lmh/atac/atac_test/atac_data/fastp/alignment
parallel -k -j 6 "
  samtools sort  {1}.sam > {1}.sort.bam    
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

samtools sort matched.sam > matched.sort.bam

parallel -j 6 "
  java -jar /home/lmh/picard/picard.jar MarkDuplicates \
         INPUT=${sample}.sort.bam \
	 OUTPUT=../rmdup/${sample}.rmdup.bam \
	 METRICS_FILE=../rmdup/${sample}.log \
	 REMOVE_DUPLICATES true \
	 VALIDATION_STRINGENCY =LENIENT 
  samtools index -@ 7 ../rmdup/{1}.rmdup.bam 
  samtools flagstat -@ 7  ../rmdup/{1}.rmdup.bam > ../rmdup/{1}.rmdup.stat 
" ::: $( ls *.sort.bam)

parallel -j 4 "
  java -jar /home/lmh/picard/picard.jar MarkDuplicates \
         INPUT={} \
         OUTPUT=../rmdup/{/.}.rmdup.bam \
         METRICS_FILE=../rmdup/{/.}.log \
         REMOVE_DUPLICATES=true \
         VALIDATION_STRINGENCY=LENIENT && \
  samtools index -@ 7 ../rmdup/{/.}.rmdup.bam && \
  samtools flagstat -@ 7 ../rmdup/{/.}.rmdup.bam > ../rmdup/{/.}.rmdup.stat
" ::: *.sort.bam


java -jar /home/lmh/picard/picard.jar MarkDuplicates \
INPUT=SRR14151155_matched.sort.bam \
OUTPUT=../rmdup/SRR14151155_matched.rmdup.bam \
METRICS_FILE=../rmdup/SRR14151155_matched.log \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT


java -jar picard.jar MarkDuplicates \
I=input.bam \
O=marked_duplicates.bam \
M=marked_dup_metrics.txt



/mnt/share/F


python add_barcode_to_sam.py SRR14151150.sam SRR14151150_3.fastq SRR14151150_matched.sam


while read id; do
python add_barcode_to_sam.py "${id}.sam" "${id}_3.fastq" "${id}_matched.sam"
done < /home/lmh/atac/atac_test/SRA/atac_ara_test.txt


less /home/lmh/atac/atac_test/atac_data/atac_ara_test.txt |while  read id;
do 
python add_barcode_to_sam.py ${id}.sam ${id}_3.fastq ${id}_matched.sam
done


python calculate_read_support.py matched.sort.bam fragment_output.txt

while read id; do
python calculate_read_support.py "${id}_matched.sort.rmdup.bam" "${id}_fragment.tsv" 
done < /home/lmh/atac/atac_test/atac_data/atac_ara_test.txt

library(ArchR)
library(Seurat)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
setwd("D:/Urinary system tumors/work/6_scATAC")
subdirs <- list.dirs(path, full.names = TRUE, recursive = FALSE)

# 
sampleNames <- sapply(subdirs, function(subdir) {
  fragments <- list.files(subdir, pattern = "fragments.tsv.gz$", full.names = TRUE)
  if (length(fragments) > 0) {
    return(basename(subdir)) 
  } else {
    return(NA)
  }
})


sampleNames <- sampleNames[!is.na(sampleNames)]

inputFiles <- file.path(subdirs, "fragments.tsv.gz")
# inputFiles <- inputFiles[-c(1,2,3,23)]

inputFiles
names(inputFiles) <- sampleNames  


print(inputFiles)
print(sampleNames)

ArrowFiles <- createArrowFiles(
  inputFiles      = inputFiles,
  sampleNames     = sampleNames, 
  minTSS          = 4,           
  minFrags        = 1000,      
  addTileMat      = TRUE,        
  addGeneScoreMat = TRUE,
  subThreading    = FALSE)    
PCa_ArchR_all_proj <- ArchRProject(
  ArrowFiles       = ArrowFiles,
  outputDirectory  = outputDirectory,
  copyArrows       = TRUE,
  geneAnnotation   = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  showLogo         = TRUE,
  threads          = getArchRThreads())
PCa_ArchR_all_proj

PCa_ArchR_all_proj <- addIterativeLSI(ArchRProj        = PCa_ArchR_all_proj, 
                                      useMatrix        = "TileMatrix", 
                                      name             = "IterativeLSI", 
                                      outlierQuantiles = NULL)
PCa_ArchR_all_proj <- addClusters(PCa_ArchR_all_proj,reducedDims     = "IterativeLSI")
PCa_ArchR_all_proj <- addUMAP(PCa_ArchR_all_proj,reducedDims         = "IterativeLSI")
PCa_ArchR_all_proj <- addImputeWeights(PCa_ArchR_all_proj)

PCa_ArchR_all_proj <- addHarmony(
  ArchRProj   = PCa_ArchR_all_proj,
  reducedDims = "IterativeLSI",
  name        = "Harmony",
  groupBy     = "Sample",
  force       = TRUE)

PCa_ArchR_all_proj <- addClusters(
  input       = PCa_ArchR_all_proj,
  reducedDims = "Harmony",
  method      = "Seurat",
  name        = "Clusters",
  resolution  = 0.6,
  maxClusters = 12,
  force       = TRUE)
scRNA <- readRDS('D:/Urinary system tumors/work/4_merge/PCa_sce_anno_second.rds')
PCa_ArchR_all_proj@projectMetadata$outputDirectory <- "D:/Urinary system tumors/work/6_scATAC/PCa"
temp <- str_split_fixed(PCa_ArchR_all_proj@sampleColData@listData$ArrowFiles,"/", n = 9)
temp[,9]
PCa_ArchR_all_proj@sampleColData@listData$ArrowFiles <- paste0("D:/Urinary system tumors/work/6_scATAC/PCa/ArrowFiles/",
                                                               temp[,9])
PCa_ArchR_all_proj <- addGeneIntegrationMatrix(
  ArchRProj   = PCa_ArchR_all_proj,
  useMatrix   = "GeneScoreMatrix",
  matrixName  = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA       = scRNA,
  groupRNA    = "cell_type",
  nameCell    = "predictedCell_Un",
  nameGroup   = "predictedGroup_Un",
  nameScore   = "predictedScore_Un"
)
groupList <- SimpleList(
  cluster = SimpleList(
    ATAC = PCa_ArchR_all_proj$cellNames,
    RNA = colnames(scRNA)
  )
)

PCa_ArchR_all_proj <- addGeneIntegrationMatrix(
  ArchRProj = PCa_ArchR_all_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scRNA,
  groupList = groupList,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
)
plotEmbedding(ArchRProj = PCa_ArchR_all_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotEmbedding(ArchRProj = PCa_ArchR_all_proj, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")

cell <- proj@cellColData@rownames[which(grepl("Epithelium",PCa_ArchR_all_proj@cellColData@listData[["predictedGroup_Co"]]))]

plotEmbedding(ArchRProj = proj_epi, colorBy = "cellColData", name = "predictedGroup_Co", embedding = "UMAP")
PCa_ArchR_all_proj@projectMetadata$outputDirectory <- "D:/Urinary system tumors/work/6_scATAC/ALL_arrow"
temp <- str_split_fixed(PCa_ArchR_all_proj@sampleColData@listData$ArrowFiles,"/", n = 9)
temp[,9]
PCa_ArchR_all_proj@sampleColData@listData$ArrowFiles <- paste0("D:/Urinary system tumors/work/6_scATAC/ALL_arrow/",
                                                               temp[,7])

PCa_ArchR_all_proj  <- addGroupCoverages(ArchRProj = PCa_ArchR_all_proj,groupBy = "predictedGroup_Co",force = T) 
PCa_ArchR_all_proj  <- addReproduciblePeakSet(ArchRProj= PCa_ArchR_all_proj ,groupBy = "predictedGroup_Co", 
                                              peakMethod = 'Tiles', method = 'p',force = T)
PCa_ArchR_all_proj  <- addPeakMatrix(PCa_ArchR_all_proj,force = T)
saveRDS(PCa_ArchR_all_proj,"PCa_ArchR_all_proj.rds")

if("Motif" %ni% names(PCa_ArchR_all_proj@peakAnnotation)){
  PCa_ArchR_all_proj <- addMotifAnnotations(ArchRProj = PCa_ArchR_all_proj, motifSet = "cisbp", name = "Motif",force = T)
}
PCa_ArchR_all_proj <- addBgdPeaks(PCa_ArchR_all_proj,force = TRUE)
PCa_ArchR_all_proj <- addDeviationsMatrix(
  ArchRProj = PCa_ArchR_all_proj, 
  peakAnnotation = "Motif",
  matrixName = "DeviationsMatrix2",
  force = T
)
getAvailableMatrices(PCa_ArchR_all_proj)
DeviationsMatrix <- ArchR::getMatrixFromProject(ArchRProj = PCa_ArchR_all_proj,useMatrix = "DeviationsMatrix2")
DeviationsMatrix2 <- as.data.frame(DeviationsMatrix@assays@data$z)

cell_types_df <- data.frame(
  Barcode = rownames(PCa_ArchR_all_proj@cellColData),
  CellType = as.character(PCa_ArchR_all_proj$predictedGroup_Co),
  stringsAsFactors = FALSE
)
unique(cell_types_df$CellType)

DeviationsMatrix2_avg <- data.frame(
  "Epithelium" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="Epithelium")]),
  "Myeloid" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="Myeloid")]),
  "Endothelium" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="Endothelium")]),
  "T_cells" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="T_cells")]),
  "Fibroblast" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="Fibroblast")]),
  "Mast_cells" = rowMeans(DeviationsMatrix2[,which(cell_types_df$CellType=="Mast_cells")]),check.names = F)

max_celltype <- apply(DeviationsMatrix2_avg, 1, function(row) {
  max_value <- max(row)
  clusters <- names(row)[which(row == max_value)]
  return(c(paste(clusters, collapse = ","), max_value))
})
max_clusters_df <- as.data.frame(t(max_celltype))
colnames(max_clusters_df) <- c("max_celltype", "max_value")

max_clusters_df$module <- rownames(DeviationsMatrix2_avg)

unique_TF <- max_clusters_df %>%
  group_by(module) %>%
  filter(max_value == max(max_value)) %>%
  ungroup()
unique_TF$max_celltype <- factor(unique_TF$max_celltype, levels = colnames(DeviationsMatrix2_avg))

unique_TF <- unique_TF[order(unique_TF$max_celltype), ]

DeviationsMatrix2_cancer <- DeviationsMatrix2_avg[unique_TF$module,]

plotVarDev <- getVarDeviations(PCa_ArchR_all_proj, name = "DeviationsMatrix2", plot = T)
plotVarDev
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = PCa_ArchR_all_proj, addDOC = FALSE)
seGroupMotif <- getGroupSE(ArchRProj = PCa_ArchR_all_proj, useMatrix = "DeviationsMatrix2", groupBy = "predictedGroup_Co")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(
  ArchRProj = PCa_ArchR_all_proj,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "DeviationsMatrix2",
  reducedDims = "IterativeLSI"
)
corGIM_MM <- correlateMatrices(
  ArchRProj = PCa_ArchR_all_proj,
  useMatrix1 = "GeneIntegrationMatrix_all",
  useMatrix2 = "DeviationsMatrix2",
  reducedDims = "IterativeLSI"
)
# For each correlation analyses, we annotate each motif with the maximum delta observed between clusters
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$DeviationsMatrix2_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$DeviationsMatrix2_name, rowData(seZ)$name), "maxDelta"]
# "we consider positive regulators as those TFs whose correlation between motif and gene score 
# (or gene expression) is greater than 0.5 with an adjusted p-value less than 0.01 and a maximum 
# inter-cluster difference in deviation z-score that is in the top quartile (Max TF Motif Delta)."
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"DeviationsMatrix2_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
#corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

df <- data.frame(corGSM_MM)
rownames(df) <- df$GeneScoreMatrix_name
gene_plot <- df[which(df$TFRegulator=="YES"), ]
library(ggplot2)
library(ggrepel)
p <- ggplot(df, aes(x=cor, y=maxDelta)) +
  geom_hline(aes(yintercept=0), color = "#999999", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0), color = "#999999", linetype="dashed", size=1) + 
  geom_point(data = df[(df$cor >0),], stroke = 0.5, size=2, shape=16, color="#E58267",alpha =0.6) + 
  geom_point(data = df[(df$cor < 0),], stroke = 0.5, size=2, shape=16,color="#6AACD0",alpha =0.6) + 
  geom_point(data = df[rownames(df) %in% rownames(gene_plot),], stroke = 0.5, size=4, shape=16,color="olivedrab3") + 
  labs(x = "Correlation To Gene Score",y = "Max TF Motif Delta", title = "") + 
  theme_bw() + #下面就是ggplot主题修饰
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black")) +
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 8, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 16)) +
  theme(legend.position = "none") + 
  geom_text_repel(data=gene_plot, aes(label=GeneScoreMatrix_name), color="black", size=4, fontface="italic",
                  arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
pdf(paste0("corGSM_MM_posTFregulators.pdf"), width=8, height=8)
p
dev.off()
# Same thing for RNA:
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"DeviationsMatrix2_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
#corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"

sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])



df <- data.frame(corGIM_MM)
rownames(df) <- df$GeneIntegrationMatrix_name
gene_plot <- df[which(df$TFRegulator=="YES"), ]
library(ggplot2)
library(ggrepel)
p <- ggplot(df, aes(x=cor, y=maxDelta)) +
  geom_hline(aes(yintercept=0), color = "#999999", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0), color = "#999999", linetype="dashed", size=1) + 
  geom_point(data = df[(df$cor >0),], stroke = 0.5, size=2, shape=16, color="#E58267",alpha =0.6) +
  geom_point(data = df[(df$cor < 0),], stroke = 0.5, size=2, shape=16,color="#6AACD0",alpha =0.6) + 
  geom_point(data = df[rownames(df) %in% rownames(gene_plot),], stroke = 0.5, size=4, shape=16,color="#117733") + 
  labs(x = "Correlation To Gene Expression",y = "Max TF Motif Delta", title = "") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black")) +
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 8, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 16)) +
  theme(legend.position = "none") + 
  geom_text_repel(data=gene_plot, aes(label=GeneIntegrationMatrix_all_name), color="black", size=4, fontface="italic", 
                  arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
pdf(paste0("corGIM_MM_posTFregulators.pdf"), width=8, height=8)
p
dev.off()


PCa_ArchR_all_proj <- addCoAccessibility(
  ArchRProj = PCa_ArchR_all_proj,
  dimsToUse = 1:30,
  reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
  ArchRProj = PCa_ArchR_all_proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA_loop <- getCoAccessibility(
  ArchRProj = PCa_ArchR_all_proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = T
)

cA_data <- as.data.frame(cA)
cA_loop_data <- as.data.frame(cA_loop)
PCa_ArchR_all_proj <- addPeak2GeneLinks(
  ArchRProj = PCa_ArchR_all_proj,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix_all"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = PCa_ArchR_all_proj,
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
color <- c("Epithelium" = "#E45A5F",
           "T_cells" = "#41B749",
           "B_cells" = "#684797",
           "Myeloid" = "#1781B5",
           "Endothelium" = "#C0937E",
           "Fibroblast" = "#F58135",
           "Mast_cells" = "#CE5E8A")
plotEmbedding(ArchRProj = PCa_ArchR_all_proj, colorBy = "cellColData", name = "predictedGroup_Co", 
              embedding = "UMAP",pal = color, size = 0.5,labelAsFactors = F)


markerPeaks <- getMarkerFeatures(
  ArchRProj = PCa_ArchR_all_proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_Co",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")


enrichMotifs <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = PCa_ArchR_all_proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
# Multiply by 50 since this is a super small test sample
heatmapEM <- plotEnrichHeatmap(enrichMotifs,cutOff = -log10(0.05),n=8,transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "PCa_motif", width = 10, height = 5, ArchRProj = PCa_ArchR_all_proj, addDOC = FALSE)


heatmapEM <- plotEnrichHeatmap(ccRCC_enrichMotifs ,cutOff = -log10(0.05),n=8,transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "ccRCC_motif", width = 10, height = 5, ArchRProj = PCa_ArchR_all_proj, addDOC = FALSE)

#decoupleR -------
library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
plan("multisession", workers = 4)
net <- get_collectri(organism='human', split_complexes=FALSE)
PCa <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/3_PCa/three_datasets/1_PCa_tumor_sce/4_PCa_tumor_sce_after_cell_type.rds")
mat <- as.matrix(PCa@assays$RNA@data)
acts <- run_ulm(mat, net,minsize = 5,
                .source='source', .target='target',.mor='mor')

PCa[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = PCa) <- "tfsulm"

library(future)
plan("multicore") 
options(future.globals.maxSize = 10 * 1024^3)  
PCa <- ScaleData(PCa)
PCa@assays$tfsulm@data <- PCa@assays$tfsulm@scale.data

n_tfs <- 25

df_PCa <- t(as.matrix(PCa@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = PCa@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))

y <- c()
for(i in unique(df$source)){
  x <- df[which(df$source == i),]
  if(max(x$mean) > x$mean[order(x$mean, decreasing = TRUE)][2]+x$mean[order(x$mean, decreasing = TRUE)][2] ){
    y2 <- c(x$cluster[which(x$mean == max(x$mean))],x$source[which(x$mean == max(x$mean))])
    y <- rbind(y, y2)
  }
}



tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)


top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()


palette_length = 100
my_color = colorRampPalette(c("Darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))


pheatmap(top_acts_mat, 
         border_color = NA,
         color=my_color, 
         breaks = my_breaks) 
saveRDS(PCa,"PCa_TF.rds")
rm("PCa")
rm("mat")
BC <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/2_BC/1_BC_tumor_sce/BC_sce_anno_second.rds")

mat <- as.matrix(BC@assays$RNA@data)

acts <- run_ulm(mat, net,minsize = 5,
                .source='source', .target='target',.mor='mor')


BC[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = BC) <- "tfsulm"

library(future)
plan("multicore") 
options(future.globals.maxSize = 10 * 1024^3)  
BC <- ScaleData(BC)
BC@assays$tfsulm@data <- BC@assays$tfsulm@scale.data


n_tfs <- 25

df <- t(as.matrix(BC@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = BC@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))


tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)


top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()


palette_length = 100
my_color = colorRampPalette(c("Darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))


pheatmap(top_acts_mat, 
         border_color = NA,
         color=my_color, 
         breaks = my_breaks) 
saveRDS(BC,"BC_TF.rds")
rm("BC")
rm("mat")

ccRCC <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/1_ccRCC/1_ccRCC_tumor_sce/ccRCC_sce_anno_second.rds")
mat <- as.matrix(ccRCC@assays$RNA@data)

acts <- run_ulm(mat, net,minsize = 5,
                .source='source', .target='target',.mor='mor')


ccRCC[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = ccRCC) <- "tfsulm"

library(future)
plan("multicore") 
options(future.globals.maxSize = 10 * 1024^3)  
ccRCC <- ScaleData(ccRCC)
ccRCC@assays$tfsulm@data <- ccRCC@assays$tfsulm@scale.data


n_tfs <- 25

df <- t(as.matrix(ccRCC@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = ccRCC@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))
tfs <- c("SPIB","SPIC","SPI1","BCL11A","BCL11B","CEBPB")

tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)


top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()


palette_length = 100
my_color = colorRampPalette(c("Darkblue","white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))


pheatmap(top_acts_mat, 
         border_color = NA,
         color=my_color, 
         breaks = my_breaks) 
saveRDS(ccRCC,"ccRCC_TF.rds")


df_PCa <- t(as.matrix(PCa@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = PCa@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))
write.table(df_PCa,"PCa_TF.txt",quote = F,sep = "\t")
y <- c()
for(i in unique(df_PCa$source)){
  x <- df_PCa[which(df_PCa$source == i),]
  if(max(x$mean) > x$mean[order(x$mean, decreasing = TRUE)][2]*3){
    y2 <- c(x$cluster[which(x$mean == max(x$mean))],x$source[which(x$mean == max(x$mean))])
    y <- rbind(y, y2)
  }
}
PCa_TFy <- y
df_PCa <- df_PCa[which(df_PCa$source %in% y[,2]),]

df_BC <- t(as.matrix(BC@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = BC@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))
write.table(df_BC,"BC_TF.txt",quote = F,sep = "\t")
y <- c()
for(i in unique(df_BC$source)){
  x <- df_BC[which(df_BC$source == i),]
  if(max(x$mean) > x$mean[order(x$mean, decreasing = TRUE)][2]*3 ){
    y2 <- c(x$cluster[which(x$mean == max(x$mean))],x$source[which(x$mean == max(x$mean))])
    y <- rbind(y, y2)
  }
}
BC_TFy <- y
df_BC <- df_BC[which(df_BC$source %in% y[,2]),]
df_ccRCC <- t(as.matrix(ccRCC@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = ccRCC@meta.data$cell_type) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source", 
               values_to = "score") %>% 
  group_by(cluster, source) %>% 
  summarise(mean = mean(score))
write.table(df_ccRCC,"ccRCC_TF.txt",quote = F,sep = "\t")
y <- c()
for(i in unique(df_ccRCC$source)){
  x <- df_ccRCC[which(df_ccRCC$source == i),]
  if(max(x$mean) > x$mean[order(x$mean, decreasing = TRUE)][2]*3 ){
    y2 <- c(x$cluster[which(x$mean == max(x$mean))],x$source[which(x$mean == max(x$mean))])
    y <- rbind(y, y2)
  }
}
ccRCC_TFy <- y
df_ccRCC <- df_ccRCC[which(df_ccRCC$source %in% y[,2]),]

#TCGA xCell -----
library(xCell)
library(data.table)
PCa_exp <- fread("TCGA-PRAD.star_tpm.tsv", data.table = F )
library(biomaRt)
library(ggplot2)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

PCa_exp$Ensembl_ID <- gsub("\\..*", "", PCa_exp$Ensembl_ID)

genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "ensembl_gene_id",
               values = PCa_exp$Ensembl_ID,
               mart = mart)

genes <- genes[genes$hgnc_symbol != "", ]

PCa_exp$ensembl_gene_id <- PCa_exp$Ensembl_ID
PCa_exp_mapped <- merge(genes, PCa_exp, by = "ensembl_gene_id")
PCa_exp_mapped <- PCa_exp_mapped[,-c(1,3)]
library(dplyr)
PCa_exp_symbol <- PCa_exp_mapped  %>%
  group_by(hgnc_symbol) %>%
  summarise_all(mean)

PCa_exp_symbol <- as.data.frame(PCa_exp_symbol)
rownames(PCa_exp_symbol) <- PCa_exp_symbol$hgnc_symbol
PCa_exp_symbol <- PCa_exp_symbol[, -1]
pd <- fread("TCGA-PRAD.clinical.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(PCa_exp_symbol),]
pd$sample <- "Tumor"
pd$rowname <- row.names(pd)
pd$sample[grep("-11",pd$rowname)] <- "Normal"
temp <- pd[,c("rowname","sample","secondary_gleason_grade.diagnoses","primary_gleason_grade.diagnoses")]
temp <- temp[,-1]
temp$gleason <- paste0(temp$primary_gleason_grade.diagnoses,"+",temp$secondary_gleason_grade.diagnoses)
temp$gleason
temp$group <- "none"
temp[which(temp$gleason %in% c("Pattern 2+Pattern 4",
                               "Pattern 3+Pattern 3",
                               "Pattern 3+Pattern 4",
                               "Pattern 3+Pattern 5",
                               "Pattern 4+Pattern 3")),]$group <- "LG"

temp[which(temp$gleason %in% c("Pattern 4+Pattern 4",
                               "Pattern 4+Pattern 5",
                               "Pattern 5+Pattern 5",
                               "Pattern 5+Pattern 3",
                               "Pattern 5+Pattern 4")),]$group <- "HG"
temp[which(temp$sample == "Normal"),"group"] <- "Normal"
PCa_temp <- temp
library(devtools)
library(MCPcounter)
PCa_tscores2 <- xCellAnalysis(PCa_exp_symbol,rnaseq = TRUE,scale = TRUE, alpha = 0.5, 
                              parallel.sz = 4, parallel.type = "SOCK")
plot <- data.frame(t(PCa_tscores2))
plot <- cbind(plot[row.names(PCa_temp),],PCa_temp)
plot$ImmuneScore
ylim1 <- boxplot.stats(plot$ImmuneScore)$stats[c(1, 5)]
ggplot(plot, aes(x = gleason, y = ImmuneScore, fill = sample)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 12) + 
  theme(panel.grid=element_blank()) +
  geom_signif(comparisons = list(c("Normal","LG"),c("Normal","HG")),
              test = t.test, 
              y_position = c(0.2,0.3),
              tip_length = c(c(0.05,0.05),c(0.05,0.05)),
              size=0.8,color="black")+ coord_cartesian(ylim = ylim1*1.05)


ccRCC_exp <- fread("D:/Urinary system tumors/work/TCGA/TCGA-KIRC.star_tpm.tsv", data.table = F )

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ccRCC_exp$Ensembl_ID <- gsub("\\..*", "", ccRCC_exp$Ensembl_ID)

genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "ensembl_gene_id",
               values = ccRCC_exp$Ensembl_ID,
               mart = mart)

genes <- genes[genes$hgnc_symbol != "", ]

ccRCC_exp$ensembl_gene_id <- ccRCC_exp$Ensembl_ID
ccRCC_exp_mapped <- merge(genes, ccRCC_exp, by = "ensembl_gene_id")
ccRCC_exp_mapped <- ccRCC_exp_mapped[,-c(1,3)]
library(dplyr)
ccRCC_exp_symbol <- ccRCC_exp_mapped  %>%
  group_by(hgnc_symbol) %>%
  summarise_all(mean)

ccRCC_exp_symbol <- as.data.frame(ccRCC_exp_symbol)
rownames(ccRCC_exp_symbol) <- ccRCC_exp_symbol$hgnc_symbol
ccRCC_exp_symbol <- ccRCC_exp_symbol[, -1]
pd <- fread("D:/Urinary system tumors/work/TCGA/TCGA-KIRC.clinical.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(ccRCC_exp_symbol),]
pd$sample <- "Tumor"
pd$rowname <- row.names(pd)
pd$sample[grep("-11",pd$rowname)] <- "Normal"
temp <- pd[,c("rowname","sample","gender.demographic","ajcc_pathologic_stage.diagnoses")]
temp <- temp[which(temp$gender.demographic == "male"),]
ccRCC_temp <- temp
ccRCC_temp[which(ccRCC_temp$sample == "Normal"),"ajcc_pathologic_stage.diagnoses"] <- "Normal"
ccRCC_temp <- ccRCC_temp[which(ccRCC_temp$ajcc_pathologic_stage.diagnoses != ""),]
ccRCC_tscores2 <- xCellAnalysis(ccRCC_exp_symbol,rnaseq = TRUE,scale = TRUE, alpha = 0.5, 
                                parallel.sz = 4, parallel.type = "SOCK")
plot <- data.frame(t(ccRCC_tscores2))
plot <- cbind(plot[row.names(ccRCC_temp),],ccRCC_temp)
plot$ImmuneScore
ylim1 <- boxplot.stats(plot$ImmuneScore)$stats[c(1, 5)]
ggplot(plot, aes(x = ajcc_pathologic_stage.diagnoses, y = ImmuneScore, fill = sample)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw(base_size = 12) + 
  theme(panel.grid=element_blank()) +
  geom_signif(comparisons = list(c("Normal","Stage I"),c("Normal","Stage II"),c("Normal","Stage III"),c("Normal","Stage IV")),
              test = t.test, 
              y_position = c(0.21,0.26,0.3,0.34),
              tip_length = c(c(0.02,0.02),c(0.02,0.02),c(0.02,0.02),c(0.02,0.02),c(0.02,0.02)),
              size=0.8,color="black")+ coord_cartesian(ylim = ylim1*1.05)


BC_exp <- fread("TCGA-BLCA.star_tpm.tsv", data.table = F )

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

BC_exp$Ensembl_ID <- gsub("\\..*", "", BC_exp$Ensembl_ID)

genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "ensembl_gene_id",
               values = BC_exp$Ensembl_ID,
               mart = mart)

genes <- genes[genes$hgnc_symbol != "", ]

BC_exp$ensembl_gene_id <- BC_exp$Ensembl_ID
BC_exp_mapped <- merge(genes, BC_exp, by = "ensembl_gene_id")
BC_exp_mapped <- BC_exp_mapped[,-c(1,3)]
library(dplyr)
BC_exp_symbol <- BC_exp_mapped  %>%
  group_by(hgnc_symbol) %>%
  summarise_all(mean)

BC_exp_symbol <- as.data.frame(BC_exp_symbol)
rownames(BC_exp_symbol) <- BC_exp_symbol$hgnc_symbol
BC_exp_symbol <- BC_exp_symbol[, -1]
pd <- fread("D:/Urinary system tumors/work/TCGA/TCGA-BLCA.clinical.tsv", data.table = F )
rownames(pd)=pd[,1]
pd <- pd[colnames(BC_exp_symbol),]
pd$sample <- "Tumor"
pd$rowname <- row.names(pd)
pd$sample[grep("-11",pd$rowname)] <- "Normal"
temp <- pd[,c("rowname","sample","gender.demographic")]
temp <- temp[which(temp$gender.demographic == "female"),]
BC_temp <- temp
annotation_col <- temp[which(temp$sample == "Tumor"),]
annotation_col2 <- temp[which(temp$sample == "Normal"),]
BC_tscores2 <- xCellAnalysis(BC_exp_symbol,rnaseq = TRUE,scale = TRUE, alpha = 0.5, 
                             parallel.sz = 4, parallel.type = "SOCK")
plot <- data.frame(t(BC_tscores2))
plot <- cbind(plot[row.names(BC_temp),],BC_temp)
plot$T.cells
ylim1 <- boxplot.stats(plot$ImmuneScore)$stats[c(1, 5)]

ggplot(plot, aes(x = sample, y = ImmuneScore, fill = sample)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw(base_size = 12) + 
  theme(panel.grid=element_blank()) +
  geom_signif(comparisons = list(c("Normal","Tumor")),
              test = t.test, 
              y_position = c(0.2),
              tip_length = c(c(0.02,0.02)),
              size=0.8,color="black")