#inferCNV -------
library(Seurat)
library(dplyr)
library(ggplot2)
library(jjAnno)
library(scRNAtoolVis)
library(rjags)
library(infercnv) 
library(AnnoProbe) 
library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
library(reticulate)
library(circlize)
library(RColorBrewer)

{
  
  sce_ccRCC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
  dim(sce_ccRCC_tumor)
  
  sce_ccRCC_cnv <-  subset(sce_ccRCC_tumor, cell_type %in% c("T_cells","Epithelium"))
  dim(sce_ccRCC_cnv)
  
  saveRDS(sce_ccRCC_cnv,"D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\sce_ccRCC_cnv_run.rds")
  
  
  sce_BC_tumor <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")
  dim(sce_BC_tumor)
  sce_BC_cnv <-  subset(sce_BC_tumor, cell_type %in% c("T_cells","Epithelium"))
  dim(sce_BC_cnv)
  saveRDS(sce_BC_cnv,"D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\sce_BC_cnv_run.rds")
  
  
  sce_BC_cnv <- readRDS("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\sce_BC_cnv_run.rds")
  group_file <- as.matrix(sce_BC_cnv@meta.data$infercnv_group)
  rownames(group_file) <- rownames(sce_BC_cnv@meta.data)
  
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(sce_BC_cnv,"counts"),
                                       annotations_file = group_file,
                                       delim = "\t",
                                       gene_order_file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\0_run_infercnv_data\\gene_pos_gencode.v35.annotation.txt",
                                       ref_group_names = "BC_normal_adjacent_Epithelium")
  
  infercnv_obj_run <- infercnv::run(infercnv_obj,
                                    out_dir = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\out_result_cnv_BC",
                                    cutoff = 0.1,
                                    analysis_mode="subclusters",
                                    cluster_by_groups = T,
                                    HMM = T,
                                    denoise = T,
                                    write_expr_matrix=TRUE,
                                    #up_to_step = 15,
                                    num_threads = 10,
                                    leiden_resolution = 0.0001,
                                    output_format = "pdf")
  
  
  
  plot_cnv(infercnv_obj_run,
           out_dir = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\plot_cnv_BC",
           title="infercnv figure",
           obs_title = "Melanoma cells",
           ref_title = "Melanocytes",
           cluster_by_groups = TRUE,
           cluster_references = TRUE,
           color_safe_pal = TRUE,
           output_filename ="infercnv_scaled_to_chr",
           output_format =  'pdf')
  
  
  plot_per_group(infercnv_obj_run,
                 on_references = FALSE,
                 on_observations = TRUE,
                 sample = FALSE,
                 out_dir = "/home/data/t220416/Melanoma/inferCNV/out_result_acral_skin/Fig_sample",
                 output_format =  'pdf')
  
  
  sce_PCa <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")
  dim(sce_PCa)
  sce_PCa_cnv <-  subset(sce_PCa, cell_type %in% c("T_cells","Epithelium"))
  dim(sce_PCa_cnv)
  saveRDS(sce_PCa_cnv,"D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\sce_PCa_cnv_run.rds")
  
  
  sce_PCa_cnv <- readRDS("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\sce_PCa_cnv_run.rds")
  group_file <- as.matrix(sce_PCa_cnv@meta.data$infercnv_group)
  rownames(group_file) <- rownames(sce_PCa_cnv@meta.data)
  
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(sce_PCa_cnv,"counts"),
                                       annotations_file = group_file,
                                       delim = "\t",
                                       gene_order_file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\0_run_infercnv_data\\gene_pos_gencode.v35.annotation.txt",
                                       ref_group_names = "PCa_normal_adjacent_Epithelium")
  
  infercnv_obj_run <- infercnv::run(infercnv_obj,
                                    out_dir = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\out_result_cnv_PCa",
                                    cutoff = 0.1,
                                    analysis_mode="subclusters",
                                    cluster_by_groups = T,
                                    HMM = T,
                                    denoise = T,
                                    write_expr_matrix=TRUE,
                                    #up_to_step = 15,
                                    leiden_resolution = 0.0001,
                                    output_format = "pdf")
  
  save(infercnv_obj,file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\infercnv_obj.RData")
  save(infercnv_obj_run,file ="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\infercnv_obj_run.RData")
  
  
  plot_cnv(infercnv_obj_run,
           out_dir = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\plot_cnv_PCa",
           title="infercnv figure",
           obs_title = "PCa_Epithelium",
           ref_title = "PCa_normal_adjacent_Epithelium",
           cluster_by_groups = TRUE,
           cluster_references = TRUE,
           color_safe_pal = TRUE,
           output_filename ="infercnv_scaled_to_chr",
           output_format =  'pdf')
  
  plot_per_group(infercnv_obj_run,
                 on_references = FALSE,
                 on_observations = TRUE,
                 sample = FALSE,
                 out_dir = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\plot_cnv_PCa\\plot_samples",
                 output_format =  'pdf')
  
}

{
  
  infercnv_obj<-readRDS('D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\out_result_cnv_ccRCC\\run.final.infercnv_obj')
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$T_cells
  
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc <- test_loc$Epithelium
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("T_cells",length(normal_loc)),rep("Epithelium",length(test_loc)))
  )
  
  
  head(anno.df)
  gn <- rownames(expr)
  geneFile <- read.table("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\0_run_infercnv_data\\gene_pos_gencode.v35.annotation.txt",header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  set.seed(42)
  kmeans.result <- kmeans(t(expr), 7)
  
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
  kmeans_df_s=arrange(kmeans_df,kmeans_class) 
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
  head(kmeans_df_s)
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(10, "Dark2")[1:7]
  names(color_v)=as.character(1:7)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("T_cells"="#41B749","Epithelium" = "#E45A5F"),kmeans_class=color_v))
  
  pdf("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\plot_cnv_ccRCC\\2_inferCNV_kmean_hot.pdf",width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
               col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
               column_gap = unit(2, "mm"),
               
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
               
               top_annotation = top_anno,left_annotation = left_anno,
               row_title = NULL,column_title = NULL)
  draw(ht,heatmap_legend_side = "right")
  dev.off()
  
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  
  p <- CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_fill_manual(values = color_v)+theme_bw()
  pdf()
  ggsave(p,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\plot_cnv_ccRCC\\2_kcluster_CNV_score.pdf"
         ,width =10,height =6, units = "cm")
  dev.off()
  
  save(expr,ht,kmeans_df,kmeans_df_s,file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\plot_cnv_ccRCC\\2_inferCNV_kmean.Rdata")
  
  
  
  scRNA_Ductal <- readRDS("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\sce_ccRCC_cnv_run.rds")
  class(scRNA_Ductal@meta.data)
  kmeans_df_s_order <- kmeans_df_s[rownames(scRNA_Ductal@meta.data),]
  identical(rownames(kmeans_df_s_order),rownames(scRNA_Ductal@meta.data))
  
  scRNA_Ductal@meta.data$kmeans_class <- kmeans_df_s_order$kmeans_class
  
  
  table(kmeans_df_s$kmeans_class,kmeans_df_s$class)
  #kmeans_class 1 3 = Epithelium_normal kmeans_class  2 4 5 6 7 = Epithelium_tumor   
  
  scRNA_Ductal@meta.data$Epithelium_type <- recode(scRNA_Ductal@meta.data$kmeans_class,
                                                   "1"="Epithelium_normal",
                                                   "3"="Epithelium_normal",
                                                   "2"="Epithelium_tumor",
                                                   "4"="Epithelium_tumor",
                                                   "5"="Epithelium_tumor",
                                                   "6"="Epithelium_tumor",
                                                   "7"="Epithelium_tumor")
  
  
  table(scRNA_Ductal@meta.data$Epithelium_type)
  
  col_Epithelium_type <- c("Epithelium_tumor" = "#DC0000",
                           "Epithelium_normal" = "#5050FF")
  DimPlot(scRNA_Ductal,reduction = "umap",group.by = "Epithelium_type",cols =col_Epithelium_type )
  
  saveRDS(scRNA_Ductal,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\1_ccRCC\\2_ccRCC_after_ident_tumor.rds")
}

#
{
  
  infercnv_obj<-readRDS('D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\out_result_cnv_BC\\run.final.infercnv_obj')
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$T_cells
  
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc <- test_loc$Epithelium
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("T_cells",length(normal_loc)),rep("Epithelium",length(test_loc)))
  )
  
  
  head(anno.df)
  gn <- rownames(expr)
  geneFile <- read.table("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\0_run_infercnv_data\\gene_pos_gencode.v35.annotation.txt",header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  set.seed(42)
  kmeans.result <- kmeans(t(expr), 7)
  
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
  kmeans_df_s=arrange(kmeans_df,kmeans_class) 
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
  head(kmeans_df_s)
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(10, "Dark2")[1:7] 
  names(color_v)=as.character(1:7)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("T_cells"="#41B749","Epithelium" = "#E45A5F"),kmeans_class=color_v))
  
  pdf("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\plot_cnv_BC\\2_inferCNV_kmean_hot.pdf",width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
               col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
               column_gap = unit(2, "mm"),
               
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
               
               top_annotation = top_anno,left_annotation = left_anno, 
               row_title = NULL,column_title = NULL)
  draw(ht,heatmap_legend_side = "right")
  dev.off()
  
  
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  
  p <- CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_fill_manual(values = color_v)+theme_bw()
  pdf()
  ggsave(p,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\plot_cnv_BC\\2_kcluster_CNV_score.pdf"
         ,width =10,height =6, units = "cm")
  dev.off()
  
  save(expr,ht,kmeans_df,kmeans_df_s,file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\plot_cnv_BC\\2_inferCNV_kmean.Rdata")
  
  
  
  scRNA_Ductal <- readRDS("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\sce_BC_cnv_run.rds")
  class(scRNA_Ductal@meta.data)
  kmeans_df_s_order <- kmeans_df_s[rownames(scRNA_Ductal@meta.data),]
  identical(rownames(kmeans_df_s_order),rownames(scRNA_Ductal@meta.data))
  
  scRNA_Ductal@meta.data$kmeans_class <- kmeans_df_s_order$kmeans_class
  
  #kmeans_class 4 7= Epithelium_normal kmeans_class 1 2 3 5 6 = Epithelium_tumor   
  
  scRNA_Ductal@meta.data$Epithelium_type <- recode(scRNA_Ductal@meta.data$kmeans_class,
                                                   "4"="Epithelium_normal",
                                                   "7"="Epithelium_normal",
                                                   "1"="Epithelium_tumor",
                                                   "2"="Epithelium_tumor",
                                                   "3"="Epithelium_tumor",
                                                   "5"="Epithelium_tumor",
                                                   "6"="Epithelium_tumor")
  
  
  table(scRNA_Ductal@meta.data$Epithelium_type)
  
  col_Epithelium_type <- c("Epithelium_tumor" = "#DC0000",
                           "Epithelium_normal" = "#5050FF")
  DimPlot(scRNA_Ductal,reduction = "umap",group.by = "Epithelium_type",cols =col_Epithelium_type )
  
  saveRDS(scRNA_Ductal,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_BC\\2_BC_after_ident_tumor.rds")
}



#PCa
{
  
  
  
  
  infercnv_obj<-readRDS('D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\out_result_cnv_PCa\\run.final.infercnv_obj')
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$T_cells
  
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc <- test_loc$Epithelium
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("T_cells",length(normal_loc)),rep("Epithelium",length(test_loc)))
  )
  
  
  head(anno.df)
  gn <- rownames(expr)
  geneFile <- read.table("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\0_run_infercnv_data\\gene_pos_gencode.v35.annotation.txt",header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
  expr=expr[intersect(gn,geneFile$V1),]
  head(sub_geneFile,4)
  expr[1:4,1:4]
  set.seed(42)
  kmeans.result <- kmeans(t(expr), 7)
  
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
  kmeans_df_s=arrange(kmeans_df,kmeans_class) 
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
  head(kmeans_df_s)
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(10, "Dark2")[1:7] 
  names(color_v)=as.character(1:7)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("T_cells"="#41B749","Epithelium" = "#E45A5F"),kmeans_class=color_v))
  
  pdf("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\plot_cnv_PCa\\2_inferCNV_kmean_hot.pdf",width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
               col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
               column_gap = unit(2, "mm"),
               
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
               
               top_annotation = top_anno,left_annotation = left_anno, 
               row_title = NULL,column_title = NULL)
  draw(ht,heatmap_legend_side = "right")
  dev.off()
  
  
  
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  
  p <- CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_fill_manual(values = color_v)+theme_bw()+ylim(0,0.0075)
  pdf()
  ggsave(p,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\plot_cnv_PCa\\2_kcluster_CNV_score.pdf"
         ,width =10,height =6, units = "cm")
  dev.off()
  
  save(expr,ht,kmeans_df,kmeans_df_s,file = "D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\plot_cnv_PCa\\2_inferCNV_kmean.Rdata")
  
  
  
  scRNA_Ductal <- readRDS("D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\3_PCa\\sce_PCa_cnv_run.rds")
  class(scRNA_Ductal@meta.data)
  kmeans_df_s_order <- kmeans_df_s[rownames(scRNA_Ductal@meta.data),]
  identical(rownames(kmeans_df_s_order),rownames(scRNA_Ductal@meta.data))
  
  
  scRNA_Ductal@meta.data$kmeans_class <- kmeans_df_s_order$kmeans_class
  
  #kmeans_class 1 3= Epithelium_normal kmeans_class 2 4 5 6 7= Epithelium_tumor   
  
  scRNA_Ductal@meta.data$Epithelium_type <- recode(scRNA_Ductal@meta.data$kmeans_class,
                                                   "1"="Epithelium_normal",
                                                   "3"="Epithelium_normal",
                                                   "2"="Epithelium_tumor",
                                                   "4"="Epithelium_tumor",
                                                   "5"="Epithelium_tumor",
                                                   "6"="Epithelium_tumor",
                                                   "7"="Epithelium_tumor")
  
  
  table(scRNA_Ductal@meta.data$Epithelium_type)
  
  col_Epithelium_type <- c("Epithelium_tumor" = "#DC0000",
                           "Epithelium_normal" = "#5050FF")
  DimPlot(scRNA_Ductal,reduction = "umap",group.by = "Epithelium_type",cols =col_Epithelium_type )
  saveRDS(scRNA_Ductal,file="D:\\Urinary system tumors\\work\\2_InferCNV_tumor\\2_PCa_after_ident_tumor.rds")
}
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

sce_PCa <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\3_PCa\\three_datasets\\1_PCa_tumor_sce\\4_PCa_tumor_sce_after_cell_type.rds")
sce_ccRCC <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\1_ccRCC\\1_ccRCC_tumor_sce\\4_ccRCC_tumor_sce_after_cell_type.rds")
sce_BC <- readRDS("D:\\Urinary system tumors\\work\\1_scRNA_landsacape\\2_BC\\1_BC_tumor_sce\\4_BC_tumor_sce_after_cell_type.rds")


sce_PCa@meta.data[row.names(PCa@meta.data),"cell_type"] <- PCa@meta.data$cell_type_second
sce_PCa <-  subset(sce_PCa, cell_type %in% c("T_cells","Tumor cells"))
group_file <- as.matrix(sce_PCa@meta.data$cell_type)
rownames(group_file) <- rownames(sce_PCa@meta.data)
##out_result_acral_skin：leiden算法，leiden_resolution为0.0001
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(sce_PCa,"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T_cells")

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 40,
                                  leiden_resolution = 0.0001,
                                  output_format = "pdf")

sce_ccRCC@meta.data[row.names(ccRCC@meta.data),"cell_type"] <- ccRCC@meta.data$cell_type_second
sce_ccRCC <-  subset(sce_ccRCC, cell_type %in% c("T_cells","Tumor cells"))
group_file <- as.matrix(sce_ccRCC@meta.data$cell_type)
rownames(group_file) <- rownames(sce_ccRCC@meta.data)
##out_result_acral_skin：leiden算法，leiden_resolution为0.0001
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(sce_ccRCC,"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T_cells")

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 40,
                                  leiden_resolution = 0.0001,
                                  output_format = "pdf")
sce_BC@meta.data[row.names(BC@meta.data),"cell_type"] <- BC@meta.data$cell_type_second
sce_BC <-  subset(sce_BC, cell_type %in% c("T_cells","Tumor cells"))
group_file <- as.matrix(sce_BC@meta.data$cell_type)
rownames(group_file) <- rownames(sce_BC@meta.data)
##out_result_acral_skin：leiden算法，leiden_resolution为0.0001
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(sce_BC,"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T_cells")

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 40,
                                  leiden_resolution = 0.0001,
                                  output_format = "pdf")
#Tumor subclone analysis ------
library(Seurat)
#ccRCC
ccRCC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor.rds")
ccRCC <-  subset(ccRCC, cell_type %in% c("Epithelium"))
ccRCC@meta.data$cell_type_second <- "Normal epithelium"
ccRCC@meta.data[which(ccRCC@meta.data$cell_type == "Epithelium"&ccRCC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
ccRCC_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC/ccRCC_infercnv.observation_groupings.txt",header = T)
ccRCC@meta.data$CNV_type <- "Normal"
ccRCC@meta.data[row.names(ccRCC_CNV), ]$CNV_type <- ccRCC_CNV$Dendrogram.Group

ccRCC_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                             sep = "\t",header = T)
ccRCC_CNV_gene<- ccRCC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                             "1"="Complete loss",
                                                             "2"="One copy_loss",
                                                             "3"="Neutral",
                                                             "4"="One copy_gain",
                                                             "5"="Two copies_gain",
                                                             "6"=">Two copies_gain"))
ccRCC_CNV_gene <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$CNV != "Neutral"),]

ccRCC_CNV_DEG <-list()
for( i in unique(ccRCC@meta.data$CNV_type)){
  z <- paste0("Tumor cells.",i)
  temp <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(ccRCC,
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  ccRCC_CNV_DEG[[i]] <- CM_diff
}

#BC
BC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/2_BC/2_BC_after_ident_tumor.rds")
BC <-  subset(BC, cell_type %in% c("Epithelium"))
BC@meta.data$cell_type_second <- "Normal epithelium"
BC@meta.data[which(BC@meta.data$cell_type == "Epithelium"&BC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
BC_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC/BC_infercnv.observation_groupings.txt",header = T)
BC@meta.data$CNV_type <- "Normal"
BC@meta.data[row.names(BC_CNV), ]$CNV_type <- BC_CNV$Dendrogram.Group

BC_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                          sep = "\t",header = T)
BC_CNV_gene<- BC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                       "1"="Complete loss",
                                                       "2"="One copy_loss",
                                                       "3"="Neutral",
                                                       "4"="One copy_gain",
                                                       "5"="Two copies_gain",
                                                       "6"=">Two copies_gain"))
BC_CNV_gene <- BC_CNV_gene[which(BC_CNV_gene$CNV != "Neutral"),]

BC_CNV_DEG <-list()
for( i in unique(BC@meta.data$CNV_type)){
  z <- paste0("Tumor cells.",i)
  temp <- BC_CNV_gene[which(BC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(BC,
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  BC_CNV_DEG[[i]] <- CM_diff
}
#PCa
PCa <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/3_PCa/2_PCa_after_ident_tumor.rds")
PCa <-  subset(PCa, cell_type %in% c("Epithelium"))
PCa@meta.data$cell_type_second <- "Normal epithelium"
PCa@meta.data[which(PCa@meta.data$cell_type == "Epithelium"&PCa@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
PCa_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa/PCa_infercnv.observation_groupings.txt",header = T)
PCa@meta.data$CNV_type <- "Normal"
PCa@meta.data[row.names(PCa_CNV), ]$CNV_type <- PCa_CNV$Dendrogram.Group

PCa_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                           sep = "\t",header = T)
PCa_CNV_gene<- PCa_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                         "1"="Complete loss",
                                                         "2"="One copy_loss",
                                                         "3"="Neutral",
                                                         "4"="One copy_gain",
                                                         "5"="Two copies_gain",
                                                         "6"=">Two copies_gain"))
PCa_CNV_gene <- PCa_CNV_gene[which(PCa_CNV_gene$CNV != "Neutral"),]

PCa_CNV_DEG <-list()
for( i in unique(PCa@meta.data$CNV_type)){
  z <- paste0("Tumor cells.",i)
  temp <- PCa_CNV_gene[which(PCa_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(PCa,
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  PCa_CNV_DEG[[i]] <- CM_diff
}

library(fgsea)         
library(data.table)    
library(ggplot2)       
library(dplyr)          
library(msigdb)         
library(GSEABase)       
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)


for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )]){
  BC_CNV_DEG[[i]]$group <- "different"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC > 0 & BC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  BC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC < 0 & BC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  BC_CNV_DEG[[i]][pos,]$group <- "LOH"
}


BC_CNV_function_group <- list()
for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )] ){
  temp <- BC_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        set.seed(123)
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  BC_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}


for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )]){
  ccRCC_CNV_DEG[[i]]$group <- "different"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC > 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC < 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "LOH"
}


ccRCC_CNV_function_group <- list()
for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )] ){
  temp <- ccRCC_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        set.seed(123)
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  ccRCC_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}


for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )]){
  PCa_CNV_DEG[[i]]$group <- "different"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC > 0 & PCa_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  PCa_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC < 0 & PCa_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  PCa_CNV_DEG[[i]][pos,]$group <- "LOH"
}

PCa_CNV_function_group <- list()
for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )] ){
  temp <- PCa_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        set.seed(123)
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  PCa_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}
save(BC_CNV_function_group,ccRCC_CNV_function_group,PCa_CNV_function_group,
     file = "CNV_chr_group_function.Rdata")

intersect( intersect(BC_CNV_function[["Tumor cells_s2"]][["chr6"]]$Description, BC_CNV_function[["Tumor cells_s1"]][["chr6"]]$Description),
           BC_CNV_function[["Tumor cells_s3"]][["chr6"]]$Description)

BC_CNV_region <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
                            header = T,sep = "\t")
BC_CNV_region <- BC_CNV_region[grepl("Tumor",BC_CNV_region$cell_group_name),]
BC_CNV_region<- BC_CNV_region %>% dplyr::mutate(CNV=recode(state,
                                                           "1"="Complete loss",
                                                           "2"="One copy_loss",
                                                           "3"="Neutral",
                                                           "4"="One copy_gain",
                                                           "5"="Two copies_gain",
                                                           "6"=">Two copies_gain"))
BC_CNV_region <- BC_CNV_region[which(BC_CNV_region$CNV != "Neutral"),]


intersect( intersect(BC_CNV_function[["Tumor cells_s2"]][["chr6"]]$Description, BC_CNV_function[["Tumor cells_s1"]][["chr6"]]$Description),
           BC_CNV_function[["Tumor cells_s3"]][["chr6"]]$Description)

BC@meta.data$tumor <- paste0("BC ",BC@meta.data$cell_type_second)
ccRCC@meta.data$tumor <- paste0("ccRCC ",ccRCC@meta.data$cell_type_second)
PCa@meta.data$tumor <- paste0("PCa ",PCa@meta.data$cell_type_second)
tumor <- merge(BC, ccRCC)
tumor <- merge(tumor, PCa)
saveRDS(tumor,"tumor.rds")
library(Seurat)
library(dplyr)
tumor <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor_function/tumor.rds")
tumor <- NormalizeData(tumor) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("orig.ident","percent.mt"))
#ScaleData()
tumor <- RunPCA(tumor, features = VariableFeatures(object = tumor), verbose = F)

DimPlot(tumor, reduction = "pca")
DimHeatmap(tumor, dims = 1:5, cells = 500, balanced = TRUE)

ElbowPlot(tumor, ndims = 50)

pct <- tumor [["pca"]]@stdev / sum( tumor [["pca"]]@stdev)* 100
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
##FindClusters
pc_num= 1:14
tumor <- FindNeighbors(tumor, dims = pc_num) %>% FindClusters(resolution = c(1:10/10)) 

Idents(tumor) <- "tumor"
tumor <- tumor %>% RunTSNE(dims = pc_num) %>% RunUMAP(dims = pc_num)
tumor@meta.data$group <- 0
tumor@meta.data[grepl("BC ",tumor@meta.data$tumor),]$group <- "BC"
tumor@meta.data[grepl("ccRCC ",tumor@meta.data$tumor),]$group <- "ccRCC"
tumor@meta.data[grepl("PCa ",tumor@meta.data$tumor),]$group <- "PCa"

tumor_color <- c("BC" = "#55AC48", 
                 "ccRCC" = "#0F68A8", 
                 "PCa" = "#EB686C")
DimPlot(tumor,group.by = "group")+scale_color_manual(values=c("BC" = "#55AC48", 
                                                              "ccRCC" = "#0F68A8", 
                                                              "PCa" = "#EB686C"))+ guides(color = FALSE)
p <- DimPlot(tumor,group.by = "tumor")+scale_color_manual(values=c("BC Tumor cells" = "#2C5B77",
                                                                   "BC Normal epithelium" = "#69BFF4",
                                                                   "ccRCC Tumor cells" = "#AA212B", 
                                                                   "ccRCC Normal epithelium" = "#FCA7AD", 
                                                                   "PCa Tumor cells" = "#278417",
                                                                   "PCa Normal epithelium" = "#B0F7A4"))+ guides(color = FALSE)

ggsave("Tumor_cluster.pdf",height = 6,width = 6)

metadata <- as.data.frame(table(tumor@meta.data[,c("tumor","group")]))
color2 <- c("BC Tumor cells" = "#2C5B77",
            "BC Normal epithelium" = "#69BFF4",
            "ccRCC Tumor cells" = "#AA212B", 
            "ccRCC Normal epithelium" = "#FCA7AD", 
            "PCa Tumor cells" = "#278417",
            "PCa Normal epithelium" = "#B0F7A4")
ggplot(metadata,aes(group,Freq,fill=tumor))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,colour="black",size=11))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =color2) +
  scale_color_manual(values = color2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Tumor_percent.pdf",height = 6,width = 6)


metadata <- as.data.frame(table(tumor@meta.data[,c("tumor","gleason")]))
color2 <- c("BC Tumor cells" = "#2C5B77",
            "BC Normal epithelium" = "#69BFF4",
            "ccRCC Tumor cells" = "#AA212B", 
            "ccRCC Normal epithelium" = "#FCA7AD", 
            "PCa Tumor cells" = "#278417",
            "PCa Normal epithelium" = "#B0F7A4")
ggplot(metadata,aes(Gleason,Freq,fill=tumor))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,colour="black",size=11))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =color2) +
  scale_color_manual(values = color2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Tumor_percent_sample.pdf",height = 6,width = 6)
library(viridis)
library(ks)
data <- cbind(tumor@reductions[["umap"]]@cell.embeddings,tumor@meta.data[,c("tumor","group")])
colnames(data) <- c("umap_x","umap_y","cell","group")

BC <- data[which(data$group == "BC"),c(1,2)]
BC$umap_x <- as.numeric(BC$umap_x)
BC$umap_y <- as.numeric(BC$umap_y)
BC <- kde(x=BC)
BC1 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
BC10 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
BC25 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
BC50 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
BC75 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])

ccRCC <- data[which(data$group == "ccRCC"),c(1,2)]
colnames(ccRCC) <- c("umap_x","umap_y")
ccRCC$umap_x <- as.numeric(ccRCC$umap_x)
ccRCC$umap_y <- as.numeric(ccRCC$umap_y)
ccRCC <- kde(x=ccRCC)
ccRCC1 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
ccRCC10 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
ccRCC25 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
ccRCC50 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
ccRCC75 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])


PCa  <- data[which(data$group == "PCa"),c(1,2)]
colnames(PCa ) <- c("umap_x","umap_y")
PCa $umap_x <- as.numeric(PCa $umap_x)
PCa $umap_y <- as.numeric(PCa $umap_y)
PCa  <- kde(x=PCa )
PCa1  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
PCa10  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
PCa25  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
PCa50  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
PCa75  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])


con.color <- 'gray40'
con.size=0.8
cl='gray40'
p+geom_path(aes(x, y), data=data.frame(BC1),col="#053851",linetype = 1,size=con.size)+
  geom_path(aes(x, y), data=data.frame(ccRCC1),col="#6B0710",linetype = 1,size=con.size)+
  geom_path(aes(x, y), data=data.frame(PCa1),col="#095B15",linetype = 1,size=con.size)
ggsave("Tumor_cluster.pdf",height = 6,width = 6)

DimPlot(tumor,group.by = "CNV_type")

load("D:/Urinary system tumors/work/2_InferCNV_tumor/3_PCa/plot_cnv_PCa/2_inferCNV_kmean.Rdata")
expr2=expr-1
expr2=expr2 ^ 2
PCa_CNV_score=as.data.frame(colMeans(expr2))
colnames(PCa_CNV_score)="CNV_score"

tumor@meta.data$CNV_score <- "none"
tumor@meta.data[row.names(BC_CNV_score),]$CNV_score <- BC_CNV_score$CNV_score
tumor@meta.data[row.names(ccRCC_CNV_score),]$CNV_score <- ccRCC_CNV_score$CNV_score
tumor@meta.data[row.names(PCa_CNV_score),]$CNV_score <- PCa_CNV_score$CNV_score
tumor@meta.data$CNV_score <- as.numeric(tumor@meta.data$CNV_score)

BC_tumor = subset(tumor,group %in% c("BC"))
FeaturePlot(BC_tumor, features = "CNV_score")+ guides(color = FALSE)
ggsave("BC_tumor_CNV_score.pdf",height = 4,width = 4)
ccRCC_tumor = subset(tumor,group %in% c("ccRCC"))
FeaturePlot(ccRCC_tumor, features = "CNV_score")+ guides(color = FALSE)
ggsave("ccRCC_tumor_CNV_score.pdf",height = 4,width = 4)
PCa_tumor = subset(tumor,group %in% c("PCa"))
FeaturePlot(PCa_tumor, features = "CNV_score")+ guides(color = FALSE)
ggsave("PCa_tumor_CNV_score.pdf",height = 4,width = 4)




metadata <- as.data.frame(table(tumor@meta.data[,c("CNV_type","group")]))
color2 <- c("Normal" = "#FCA7AD",
            "Tumor cells_s1" = "#C8D3C3",
            "Tumor cells_s2" = "#EFB36E",
            "Tumor cells_s3" = "#BBDC76",
            "Tumor cells_s4" = "#D9A0AB",
            "Tumor cells_s5" = "#FFED6F",
            "Tumor cells_s6" = "#F0D1E1",
            "Tumor cells_s7" = "#8DD3C7",
            "Tumor cells_s8" = "#C8A7C9",
            "Tumor cells_s9" = "#A9A0B2",
            "Tumor cells_s10" = "#F0EFBB")
ggplot(metadata,aes(group,Freq,fill=CNV_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,colour="black",size=11))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =color2) +
  scale_color_manual(values = color2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("ccRCC_percent.pdf",height = 6,width = 6)


color2 <- c("Normal" = "#69BFF4",
            "Tumor cells_s1" = "#8DD3C7",
            "Tumor cells_s2" = "#FFED6F",
            "Tumor cells_s3" = "#D8C965")
ggplot(metadata,aes(group,Freq,fill=CNV_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,colour="black",size=11))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =color2) +
  scale_color_manual(values = color2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("BC_percent.pdf",height = 6,width = 6)

color2 <- c("Normal" = "#B0F7A4",
            "Tumor cells_s1" = "#FFED6F",
            "Tumor cells_s2" = "#8DD3C7",
            "Tumor cells_s3" = "#D8C965")
ggplot(metadata,aes(group,Freq,fill=CNV_type))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_classic()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45,hjust = 1,colour="black",size=11))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =color2) +
  scale_color_manual(values = color2)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("PCa_percent.pdf",height = 6,width = 6)
saveRDS(tumor,"D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor2.rds")


BC_CNV_function_subtype <- list()
for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )] ){
  temp <- BC_CNV_DEG[[i]]
  temp3 <- temp[which(temp$group %in% c("AMP","LOH") ),]
  x = temp3[order(temp3$avg_log2FC,decreasing = T),]
  genelist <- structure(x$avg_log2FC,names=x$gene)
  set.seed(123)
  res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
  res_CM_2 <- data.frame(res_CM)
  res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
  res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
  res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
  res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
  BC_CNV_function_subtype[[i]] <- res_CM_2
}

ccRCC_CNV_function_subtype <- list()
for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )] ){
  temp <- ccRCC_CNV_DEG[[i]]
  temp3 <- temp[which(temp$group %in% c("AMP","LOH") ),]
  x = temp3[order(temp3$avg_log2FC,decreasing = T),]
  genelist <- structure(x$avg_log2FC,names=x$gene)
  set.seed(123)
  res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
  res_CM_2 <- data.frame(res_CM)
  res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
  res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
  res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
  res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
  ccRCC_CNV_function_subtype[[i]] <- res_CM_2
}
PCa_CNV_function_subtype <- list()
for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )] ){
  temp <- PCa_CNV_DEG[[i]]
  temp3 <- temp[which(temp$group %in% c("AMP","LOH") ),]
  x = temp3[order(temp3$avg_log2FC,decreasing = T),]
  genelist <- structure(x$avg_log2FC,names=x$gene)
  set.seed(123)
  res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
  res_CM_2 <- data.frame(res_CM)
  res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
  res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
  res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
  res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
  PCa_CNV_function_subtype[[i]] <- res_CM_2
}

save(BC_CNV_function_subtype,ccRCC_CNV_function_subtype,PCa_CNV_function_subtype,
     file = "CNV_subtype_function.Rdata")

library(enrichplot)

gseaplot2(ccRCC_CNV_function_group[["Tumor cells_s2"]][["chr3"]][["LOH"]], geneSetID = c(1,2,4))
gseaplot2(ccRCC_CNV_function_group[["Tumor cells_s2"]][["chr7"]][["AMP"]], geneSetID = c(7,11,15))

gseaplot2(BC_CNV_function_group[["Tumor cells_s1"]][["chr6"]][["LOH"]], geneSetID = c(1,2,3,4))

gseaplot2(PCa_CNV_function_group[["Tumor cells_s3"]][["chr7"]][["AMP"]], geneSetID = c(1))
gseaplot2(PCa_CNV_function_group[["Tumor cells_s3"]][["chr19"]][["AMP"]], geneSetID = c(1,2))

x <- PCa_CNV_function_group
for(i in names(x) ){
  for(z in names(x[[i]])){
    for(n in names(x[[i]][[z]] )){
      if(dim(x[[i]][[z]][[n]]@result)[1] != 0){
        write.table(x[[i]][[z]][[n]]@result, paste0("PCa", " ",i, " ",z," ",n,".txt"), sep = '\t',quote = F, row.names = F)
      }
    }
  }
}

x <- BC_CNV_function_group
for(i in names(x) ){
  for(z in names(x[[i]])){
    for(n in names(x[[i]][[z]] )){
      if(dim(x[[i]][[z]][[n]]@result)[1] != 0){
        write.table(x[[i]][[z]][[n]]@result, paste0("BC", " ",i, " ",z," ",n,".txt"), sep = '\t',quote = F, row.names = F)
      }
    }
  }
}

x <- ccRCC_CNV_function_group
for(i in names(x) ){
  for(z in names(x[[i]])){
    for(n in names(x[[i]][[z]] )){
      if(dim(x[[i]][[z]][[n]]@result)[1] != 0){
        write.table(x[[i]][[z]][[n]]@result, paste0("ccRCC", " ",i, " ",z," ",n,".txt"), sep = '\t',quote = F, row.names = F)
      }
    }
  }
}
x <- ccRCC_CNV_function_subtype
for(i in names(x) ){
  write.table(x[[i]], paste0("ccRCC_GSEA_res", " ",i,".txt"), sep = '\t',quote = F, row.names = F)
}
x <- BC_CNV_function_subtype
for(i in names(x) ){
  write.table(x[[i]], paste0("BC_GSEA_res", " ",i,".txt"), sep = '\t',quote = F, row.names = F)
}
x <- PCa_CNV_function_subtype
for(i in names(x) ){
  write.table(x[[i]], paste0("PCa_GSEA_res", " ",i,".txt"), sep = '\t',quote = F, row.names = F)
}


data <- BC_CNV_function_subtype[["Tumor cells_s1"]]
data <- data[1:20,]

data$type_order=factor(rev(1:length(rev(data$Description))),labels=rev(data$Description))

ggplot(data=data, aes(x=type_order,y = -log10(p.adjust), fill=-log10(p.adjust) )) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() + 
  theme_bw() + scale_fill_gradient(
    low = "#20518B",   
    high = "#F6F6E2")
ggsave("BC_Tumor cells_s1_GSEA.pdf", height = 6,width = 18)



data <- ccRCC_CNV_function_subtype[["Tumor cells_s2"]]
data <- data[1:20,]
data$type_order=factor(rev(1:length(rev(data$Description))),labels=rev(data$Description))

ggplot(data=data, aes(x=type_order,y = -log10(p.adjust), fill=-log10(p.adjust) )) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() + 
  theme_bw() + scale_fill_gradient(
    low = "#4E081B",   
    high = "#C55557")
ggsave("ccRCC_Tumor cells_s2_GSEA.pdf", height = 6,width = 12)



data <- PCa_CNV_function_subtype[["Tumor cells_s3"]]
data <- data[1:20,]

data$type_order=factor(rev(1:length(rev(data$Description))),labels=rev(data$Description))

ggplot(data=data, aes(x=type_order,y = -log10(p.adjust), fill=-log10(p.adjust) )) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() + 
  theme_bw() + scale_fill_gradient(
    low = "#085F5A", 
    high = "#90D0D0")
ggsave("PCa_Tumor cells_s3_GSEA.pdf", height = 6,width = 12)

tumor <-readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor2.rds")
subclone_color <- c("ccRCC Normal" = "#FCA7AD",
                    "ccRCC Tumor cells_s1" = "#C8D3C3",
                    "ccRCC Tumor cells_s2" = "#EFB36E",
                    "ccRCC Tumor cells_s3" = "#BBDC76",
                    "ccRCC Tumor cells_s4" = "#D9A0AB",
                    "ccRCC Tumor cells_s5" = "#FC56DC",#FFED6F
                    "ccRCC Tumor cells_s6" = "#F0D1E1",
                    "ccRCC Tumor cells_s7" = "#FF9F06",#8DD3C7
                    "ccRCC Tumor cells_s8" = "#C8A7C9",
                    "ccRCC Tumor cells_s9" = "#A9A0B2",
                    "ccRCC Tumor cells_s10" = "#F0EFBB",
                    "BC Normal" = "#69BFF4",
                    "BC Tumor cells_s1" = "#3972FF",#FFED6F
                    "BC Tumor cells_s2" = "#A7A7EA",#8DD3C7
                    "BC Tumor cells_s3" = "#B0FFFB",#D8C965
                    "PCa Normal" = "#B0F7A4",
                    "PCa Tumor cells_s1" = "#FFED6F",
                    "PCa Tumor cells_s2" = "#8DD3C7",
                    "PCa Tumor cells_s3" = "#D8C965")
tumor@meta.data$subclone <- paste0(tumor@meta.data$group," ",tumor@meta.data$CNV_type)
p <- DimPlot(tumor,group.by = "subclone")+scale_color_manual(values=subclone_color)+ guides(color = FALSE)
library(viridis)
library(ks)
data <- cbind(tumor@reductions[["umap"]]@cell.embeddings,tumor@meta.data[,c("tumor","group")])
colnames(data) <- c("umap_x","umap_y","cell","group")

BC <- data[which(data$group == "BC"),c(1,2)]
BC$umap_x <- as.numeric(BC$umap_x)
BC$umap_y <- as.numeric(BC$umap_y)
BC <- kde(x=BC)
BC1 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
BC10 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
BC25 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
BC50 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
BC75 <- with(BC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])

ccRCC <- data[which(data$group == "ccRCC"),c(1,2)]
colnames(ccRCC) <- c("umap_x","umap_y")
ccRCC$umap_x <- as.numeric(ccRCC$umap_x)
ccRCC$umap_y <- as.numeric(ccRCC$umap_y)
ccRCC <- kde(x=ccRCC)
ccRCC1 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
ccRCC10 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
ccRCC25 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
ccRCC50 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
ccRCC75 <- with(ccRCC, contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])

PCa  <- data[which(data$group == "PCa"),c(1,2)]
colnames(PCa ) <- c("umap_x","umap_y")
PCa $umap_x <- as.numeric(PCa $umap_x)
PCa $umap_y <- as.numeric(PCa $umap_y)
PCa  <- kde(x=PCa )
PCa1  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["1%"])[[1]])
PCa10  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["10%"])[[1]])
PCa25  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["25%"])[[1]])
PCa50  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["50%"])[[1]])
PCa75  <- with(PCa , contourLines(x=eval.points[[1]], y=eval.points[[2]],z=estimate, levels=cont["75%"])[[1]])

con.color <- 'gray40'
con.size=0.8
cl='gray40'
p+geom_path(aes(x, y), data=data.frame(BC1),col="#053851",linetype = 1,size=con.size)+
  geom_path(aes(x, y), data=data.frame(ccRCC1),col="#6B0710",linetype = 1,size=con.size)+
  geom_path(aes(x, y), data=data.frame(PCa1),col="#095B15",linetype = 1,size=con.size)
ggsave("Tumor_cluster.pdf",height = 6,width = 6)

ccRCC_Spatial <- readRDS("D:/Urinary system tumors/work/1_scRNA_landsacape/sp/ccRCC_Spatial.rds")
LG_3 <- ccRCC_Spatial[["LG_3"]]
DefaultAssay(LG_3) <- "SCT"
SpatialFeaturePlot(LG_3,features = "PDK4",pt.size.factor = 4)
ggsave("ccRCC_LG_PDK4.pdf",width = 5,height = 5)
HG_2 <- ccRCC_Spatial[["HG_2"]]
DefaultAssay(HG_2) <- "SCT"
SpatialFeaturePlot(HG_2,features = "PDK4",pt.size.factor = 4)
ggsave("ccRCC_HG_PDK4.pdf",width = 5,height = 5)


SpatialFeaturePlot(LG_3,features = "HSPB1",pt.size.factor = 4)
ggsave("ccRCC_LG_HSPB1.pdf",width = 5,height = 5)
SpatialFeaturePlot(HG_2,features = "HSPB1",pt.size.factor = 4)
ggsave("ccRCC_HG_HSPB1.pdf",width = 5,height = 5)

LG <- PCa_spaital[["Tumor08"]]
DefaultAssay(LG) <- "SCT"
LG <- AddModuleScore(LG,features = list,name = names(list))
SpatialFeaturePlot(LG,features = "KLK3",pt.size.factor = 2)
ggsave("PCa_LG_KLK3.pdf",width = 5,height = 5)
SpatialFeaturePlot(LG,features = "AZGP1",pt.size.factor = 2)
ggsave("PCa_LG_AZGP1.pdf",width = 5,height = 5)

HG <- PCa_spaital[["Tumor02"]]
DefaultAssay(HG) <- "SCT"
HG <- AddModuleScore(HG,features = list,name = names(list))
SpatialFeaturePlot(HG,features = "KLK3",pt.size.factor = 2)
ggsave("PCa_HG_KLK3.pdf",width = 5,height = 5)

SpatialFeaturePlot(HG,features = "AZGP1",pt.size.factor = 2)
ggsave("PCa_HG_AZGP1.pdf",width = 5,height = 5)

##Tumor cells vs Normal epithelium ----
library(Seurat)
#ccRCC
ccRCC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/1_ccRCC/2_ccRCC_after_ident_tumor.rds")
ccRCC <-  subset(ccRCC, cell_type %in% c("Epithelium"))
ccRCC@meta.data$cell_type_second <- "Normal epithelium"
ccRCC@meta.data[which(ccRCC@meta.data$cell_type == "Epithelium"&ccRCC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
ccRCC_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC/ccRCC_infercnv.observation_groupings.txt",header = T)
ccRCC@meta.data$CNV_type <- "Normal"
ccRCC@meta.data[row.names(ccRCC_CNV), ]$CNV_type <- ccRCC_CNV$Dendrogram.Group

ccRCC_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                             sep = "\t",header = T)
ccRCC_CNV_gene<- ccRCC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                             "1"="Complete loss",
                                                             "2"="One copy_loss",
                                                             "3"="Neutral",
                                                             "4"="One copy_gain",
                                                             "5"="Two copies_gain",
                                                             "6"=">Two copies_gain"))
ccRCC_CNV_gene <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$CNV != "Neutral"),]
temp <- ccRCC_CNV_gene
CM_diff <- FindMarkers(ccRCC,
                       only.pos = FALSE,
                       logfc.threshold = 0,
                       min.pct = 0,
                       ident.1 = "Tumor cells",
                       ident.2 = "Normal epithelium",
                       group.by = "cell_type_second")
ccRCC_CM_diff2 <- CM_diff
CM_diff$gene <- row.names(CM_diff)
CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
ccRCC_CNV_DEG_all <- CM_diff

PCa <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/3_PCa/2_PCa_after_ident_tumor.rds")
PCa <-  subset(PCa, cell_type %in% c("Epithelium"))
PCa@meta.data$cell_type_second <- "Normal epithelium"
PCa@meta.data[which(PCa@meta.data$cell_type == "Epithelium"&PCa@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
PCa_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa/PCa_infercnv.observation_groupings.txt",header = T)
PCa@meta.data$CNV_type <- "Normal"
PCa@meta.data[row.names(PCa_CNV), ]$CNV_type <- PCa_CNV$Dendrogram.Group

PCa_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                           sep = "\t",header = T)
PCa_CNV_gene<- PCa_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                         "1"="Complete loss",
                                                         "2"="One copy_loss",
                                                         "3"="Neutral",
                                                         "4"="One copy_gain",
                                                         "5"="Two copies_gain",
                                                         "6"=">Two copies_gain"))
PCa_CNV_gene <- PCa_CNV_gene[which(PCa_CNV_gene$CNV != "Neutral"),]
temp <- PCa_CNV_gene
CM_diff <- FindMarkers(PCa,
                       only.pos = FALSE,
                       logfc.threshold = 0,
                       min.pct = 0,
                       ident.1 = "Tumor cells",
                       ident.2 = "Normal epithelium",
                       group.by = "cell_type_second")
PCa_CM_diff2 <- CM_diff
CM_diff$gene <- row.names(CM_diff)
CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
PCa_CNV_DEG_all <- CM_diff

BC <- readRDS("D:/Urinary system tumors/work/2_InferCNV_tumor/2_BC/2_BC_after_ident_tumor.rds")
BC <-  subset(BC, cell_type %in% c("Epithelium"))
BC@meta.data$cell_type_second <- "Normal epithelium"
BC@meta.data[which(BC@meta.data$cell_type == "Epithelium"&BC@meta.data$Epithelium_type == "Epithelium_tumor"),]$cell_type_second <- "Tumor cells"
BC_CNV <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC/BC_infercnv.observation_groupings.txt",header = T)
BC@meta.data$CNV_type <- "Normal"
BC@meta.data[row.names(BC_CNV), ]$CNV_type <- BC_CNV$Dendrogram.Group

BC_CNV_gene <- read.table("D:/Urinary system tumors/work/2_InferCNV_tumor_function/BC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                          sep = "\t",header = T)
BC_CNV_gene<- BC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                       "1"="Complete loss",
                                                       "2"="One copy_loss",
                                                       "3"="Neutral",
                                                       "4"="One copy_gain",
                                                       "5"="Two copies_gain",
                                                       "6"=">Two copies_gain"))
BC_CNV_gene <- BC_CNV_gene[which(BC_CNV_gene$CNV != "Neutral"),]
temp <- BC_CNV_gene
CM_diff <- FindMarkers(BC,
                       only.pos = FALSE,
                       logfc.threshold = 0,
                       min.pct = 0,
                       ident.1 = "Tumor cells",
                       ident.2 = "Normal epithelium",
                       group.by = "cell_type_second")
BC_CM_diff2 <- CM_diff
CM_diff$gene <- row.names(CM_diff)
CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
BC_CNV_DEG_all <- CM_diff

library(fgsea)         
library(data.table)     
library(ggplot2)      
library(dplyr)        
library(msigdb)      
library(GSEABase)       
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)

temp <- ccRCC_CNV_DEG_all
temp <- temp[order(temp$avg_log2FC,decreasing = T),]
temp <- distinct(temp[,-c(7,11)])
genelist <- structure(temp$avg_log2FC,names=temp$gene)
set.seed(123)
res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)
res_CM_2 <- data.frame(res_CM)
res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
ccRCC_res_CM_2 <- res_CM_2

temp <- BC_CNV_DEG_all
temp <- temp[order(temp$avg_log2FC,decreasing = T),]
temp <- distinct(temp[,-c(7,11)])
genelist <- structure(temp$avg_log2FC,names=temp$gene)
set.seed(123)
res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)

res_CM_2 <- data.frame(res_CM)
res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
BC_res_CM_2 <- res_CM_2




temp <- PCa_CNV_DEG_all
temp <- temp[order(temp$avg_log2FC,decreasing = T),]
temp <- distinct(temp[,-c(7,11)])
genelist <- structure(temp$avg_log2FC,names=temp$gene)
set.seed(123)
res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0,seed = T)

res_CM_2 <- data.frame(res_CM)
res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
PCa_res_CM_2 <- res_CM_2


library(ggplot2)
library(ggnewscale)
library(scales)
squash_axis <- function(from, to, factor) { 
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

pathway <- c("GOBP_CELL_CELL_ADHESION",
             "HALLMARK_GLYCOLYSIS",
             "GOBP_INFLAMMATORY_RESPONSE",
             "GOBP_EPITHELIUM_DEVELOPMENT",
             "GOBP_CELL_MIGRATION",
             "GOBP_CELL_CELL_ADHESION",
             "HALLMARK_ANDROGEN_RESPONSE",
             "GOBP_HORMONE_METABOLIC_PROCESS",
             "HALLMARK_ESTROGEN_RESPONSE_LATE",
             "HALLMARK_ESTROGEN_RESPONSE_EARLY",
             "HALLMARK_INFLAMMATORY_RESPONSE",
             "GOBP_INFLAMMATORY_RESPONSE",
             "HALLMARK_IL6_JAK_STAT3_SIGNALING",
             "HALLMARK_IL2_STAT5_SIGNALING")
ggplot(PCa_res_CM_2, aes(x =NES, y=log10adjp)) +
  geom_point(aes(size=PCa_res_CM_2$count_ratio,color=PCa_res_CM_2$NES,alpha=abs(PCa_res_CM_2$NES))) + 
  scale_alpha_continuous(range = c(0.3, 0.7))+
  scale_size(range=c(1,14))+
  scale_color_gradient2(low = "#7398CF", mid = "white", high = "#F3A29B", midpoint = 0,space = "Lab") +
  xlim(c(-3, 3))+ylim(c(0, 10)) + 
  geom_vline(xintercept=0,lty=1,col="black",lwd=0.8) + 
  geom_hline(yintercept =0, lty=1,col="black",lwd=0.8) + 
  labs(x="NES", y="-log10 adjust p value") + 
  ggtitle("PCa") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  coord_trans(x = squash_axis(-1,1,5))+      
  theme(legend.position = 'none')+
  geom_text(data = subset(PCa_res_CM_2, Description %in% pathway), aes(label = Description),size = 3)+
  geom_point(
    data = subset(PCa_res_CM_2, Description %in% pathway),
    color = "red", size = 5)
ggsave("PCa_all_gsea.pdf", width = 8,height = 6)



ggplot(ccRCC_res_CM_2, aes(x =NES, y=log10adjp)) +
  geom_point(aes(size=ccRCC_res_CM_2$count_ratio,color=ccRCC_res_CM_2$NES,alpha=abs(ccRCC_res_CM_2$NES))) + 
  scale_alpha_continuous(range = c(0.3, 0.7))+
  scale_size(range=c(1,14))+
  scale_color_gradient2(low = "#7398CF", mid = "white", high = "#F3A29B", midpoint = 0,space = "Lab") +
  xlim(c(-3, 3))+ylim(c(0, 10)) + 
  geom_vline(xintercept=0,lty=1,col="black",lwd=0.8) + 
  geom_hline(yintercept =0, lty=1,col="black",lwd=0.8) + 
  labs(x="NES", y="-log10 adjust p value") + 
  ggtitle("ccRCC") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  coord_trans(x = squash_axis(-1,1,5))+   
  theme(legend.position = 'none')+
  geom_text(data = subset(ccRCC_res_CM_2, Description %in% pathway), aes(label = Description),size = 3)+
  geom_point(
    data = subset(ccRCC_res_CM_2, Description %in% pathway),
    color = "red", size = 5)
ggsave("ccRCC_all_gsea.pdf", width = 8,height = 6)

ggplot(BC_res_CM_2, aes(x =NES, y=log10adjp)) +
  geom_point(aes(size=BC_res_CM_2$count_ratio,color=BC_res_CM_2$NES,alpha=abs(BC_res_CM_2$NES))) + 
  scale_alpha_continuous(range = c(0.3, 0.7))+
  scale_size(range=c(1,14))+
  scale_color_gradient2(low = "#7398CF", mid = "white", high = "#F3A29B", midpoint = 0,space = "Lab") +
  xlim(c(-3, 3))+ylim(c(0, 10)) + 
  geom_vline(xintercept=0,lty=1,col="black",lwd=0.8) + 
  geom_hline(yintercept =0, lty=1,col="black",lwd=0.8) + 
  labs(x="NES", y="-log10 adjust p value") + 
  ggtitle("BC") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  coord_trans(x = squash_axis(-1,1,5))+      
  theme(legend.position = 'none')+
  geom_text(data = subset(BC_res_CM_2, Description %in% pathway), aes(label = Description),size = 3)+
  geom_point(
    data = subset(BC_res_CM_2, Description %in% pathway),
    color = "red", size = 5)
ggsave("BC_all_gsea.pdf", width = 8,height = 6)
write.table(PCa_res_CM_2,("PCa_res_CM_2.txt"), sep = '\t',quote = F, row.names = F)
write.table(PCa_res_CM_2,("BC_res_CM_2.txt"), sep = '\t',quote = F, row.names = F)
write.table(PCa_res_CM_2,("ccRCC_res_CM_2.txt"), sep = '\t',quote = F, row.names = F)


##inferCNV chrx and chrY-----
BC_scRNA <- readRDS("./ALL_data/BC_sce_anno_second.rds")
BC_scRNA <-  subset(BC_scRNA, cell_type %in% c("T_cells") | cell_type_second %in% c("Tumor cells"))
BC_scRNA$infercnv_cluster <- ifelse(BC_scRNA$cell_type == "T_cells", "T cells",
                                    ifelse(BC_scRNA$cell_type == "Epithelium", "Tumor cells", "unknown"))
BC_group <- read.table("./ALL_data/BC_infercnv.observation_groupings.txt",header = T)
BC_scRNA@meta.data[row.names(BC_group),"infercnv_cluster"] <- BC_group$Dendrogram.Group
group_file <- as.matrix(BC_scRNA@meta.data$infercnv_cluster)
rownames(group_file) <- rownames(BC_scRNA@meta.data)
matrix <- as.matrix(GetAssayData(BC_scRNA[["RNA"]],"counts"))

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(BC_scRNA[["RNA"]],"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./0_run_infercnv_data/gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T cells",
                                     chr_exclude = c("chr1", "chr2", "chr3","chr4", "chr5", "chr6",
                                                     "chr7", "chr8","chr9", "chr10", "chr11", "chr12",
                                                     "chr13", "chr14","chr15", "chr16", "chr17", "chr18",
                                                     "chr19", "chr20","chr21", "chr22", "chrM"))

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "./inferCNV_XY/out_result_cnv_BC",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 10,
                                  plot_steps = FALSE,
                                  leiden_resolution = 0.0001,
                                  #no_plot = TRUE,
                                  output_format = "pdf"
                                  # resume_mode = FALSE
)
save(infercnv_obj,infercnv_obj_run,file ="./inferCNV_XY/BC_infercnv_obj_run.RData")




ccRCC_scRNA <- readRDS("./ALL_data/ccRCC_sce_anno_second2.rds")
ccRCC_scRNA <-  subset(ccRCC_scRNA, cell_type %in% c("T_cells") | cell_type_second %in% c("Tumor cells"))
ccRCC_scRNA$infercnv_cluster <- ifelse(ccRCC_scRNA$cell_type == "T_cells", "T cells",
                                       ifelse(ccRCC_scRNA$cell_type == "Epithelium", "Tumor cells", "unknown"))
ccRCC_group <- read.table("./ALL_data/ccRCC_infercnv.observation_groupings.txt",header = T)
ccRCC_scRNA@meta.data[row.names(ccRCC_group),"infercnv_cluster"] <- ccRCC_group$Dendrogram.Group
group_file <- as.matrix(ccRCC_scRNA@meta.data$infercnv_cluster)
rownames(group_file) <- rownames(ccRCC_scRNA@meta.data)
matrix <- as.matrix(GetAssayData(ccRCC_scRNA[["RNA"]],"counts"))

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(ccRCC_scRNA[["RNA"]],"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./0_run_infercnv_data/gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T cells",
                                     chr_exclude = c("chr1", "chr2", "chr3","chr4", "chr5", "chr6",
                                                     "chr7", "chr8","chr9", "chr10", "chr11", "chr12",
                                                     "chr13", "chr14","chr15", "chr16", "chr17", "chr18",
                                                     "chr19", "chr20","chr21", "chr22", "chrM"))

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "./inferCNV_XY/out_result_cnv_ccRCC",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 10,
                                  plot_steps = FALSE,
                                  leiden_resolution = 0.0001,
                                  #no_plot = TRUE,
                                  output_format = "pdf"
                                  #,resume_mode = FALSE
)
saveRDS(infercnv_obj_run,"./inferCNV_XY/ccRCC_infercnv_obj_run.rds")

PCa_scRNA <- readRDS("./ALL_data/PCa_sce_anno_second.rds")
PCa_scRNA <-  subset(PCa_scRNA, cell_type %in% c("T_cells") | cell_type_second %in% c("Tumor cells"))
PCa_scRNA$infercnv_cluster <- ifelse(PCa_scRNA$cell_type == "T_cells", "T cells",
                                     ifelse(PCa_scRNA$cell_type == "Epithelium", "Tumor cells", "unknown"))
PCa_group <- read.table("./ALL_data/PCa_infercnv.observation_groupings.txt",header = T)
PCa_scRNA@meta.data[row.names(PCa_group),"infercnv_cluster"] <- PCa_group$Dendrogram.Group
group_file <- as.matrix(PCa_scRNA@meta.data$infercnv_cluster)
rownames(group_file) <- rownames(PCa_scRNA@meta.data)
matrix <- as.matrix(GetAssayData(PCa_scRNA[["RNA"]],"counts"))

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = GetAssayData(PCa_scRNA[["RNA"]],"counts"),
                                     annotations_file = group_file,
                                     delim = "\t",
                                     gene_order_file = "./0_run_infercnv_data/gene_pos_gencode.v35.annotation.txt",
                                     ref_group_names = "T cells",
                                     chr_exclude = c("chr1", "chr2", "chr3","chr4", "chr5", "chr6",
                                                     "chr7", "chr8","chr9", "chr10", "chr11", "chr12",
                                                     "chr13", "chr14","chr15", "chr16", "chr17", "chr18",
                                                     "chr19", "chr20","chr21", "chr22", "chrM"))

infercnv_obj_run <- infercnv::run(infercnv_obj,
                                  out_dir = "./inferCNV_XY/out_result_cnv_PCa",
                                  cutoff = 0.1,
                                  analysis_mode="subclusters",
                                  cluster_by_groups = T,
                                  HMM = T,
                                  denoise = T,
                                  write_expr_matrix=TRUE,
                                  #up_to_step = 15,
                                  num_threads = 10,
                                  plot_steps = FALSE,
                                  leiden_resolution = 0.0001,
                                  no_plot = TRUE,
                                  output_format = NA
                                  #,resume_mode = FALSE
)
saveRDS(infercnv_obj_run,"./inferCNV_XY/PCa_infercnv_obj_run.rds")
#plot
#ccRCC
ccRCC_CNV_gene <- read.table("D:/Urinary system tumors/work/inferCNV_XY/out_result_cnv_ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                             sep = "\t",header = T)
ccRCC_CNV_gene<- ccRCC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                             "1"="Complete loss",
                                                             "2"="One copy_loss",
                                                             "3"="Neutral",
                                                             "4"="One copy_gain",
                                                             "5"="Two copies_gain",
                                                             "6"=">Two copies_gain"))
ccRCC_CNV_gene <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$CNV != "Neutral"),]

ccRCC_CNV_DEG <-list()
for( i in unique(ccRCC@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(ccRCC,
                         features = unique(ccRCC_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  ccRCC_CNV_DEG[[i]] <- CM_diff
}

#BC
BC_CNV_gene <- read.table("D:/Urinary system tumors/work/inferCNV_XY/out_result_cnv_BC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                          sep = "\t",header = T)
BC_CNV_gene<- BC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                       "1"="Complete loss",
                                                       "2"="One copy_loss",
                                                       "3"="Neutral",
                                                       "4"="One copy_gain",
                                                       "5"="Two copies_gain",
                                                       "6"=">Two copies_gain"))
BC_CNV_gene <- BC_CNV_gene[which(BC_CNV_gene$CNV != "Neutral"),]

BC_CNV_DEG <-list()
for( i in unique(BC@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- BC_CNV_gene[which(BC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(BC,
                         features = unique(BC_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  BC_CNV_DEG[[i]] <- CM_diff
}

#PCa
PCa_CNV_gene <- read.table("D:/Urinary system tumors/work/inferCNV_XY/out_result_cnv_PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                           sep = "\t",header = T)
PCa_CNV_gene<- PCa_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                         "1"="Complete loss",
                                                         "2"="One copy_loss",
                                                         "3"="Neutral",
                                                         "4"="One copy_gain",
                                                         "5"="Two copies_gain",
                                                         "6"=">Two copies_gain"))
PCa_CNV_gene <- PCa_CNV_gene[which(PCa_CNV_gene$CNV != "Neutral"),]

PCa_CNV_DEG <-list()
for( i in unique(PCa@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- PCa_CNV_gene[which(PCa_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(PCa,
                         features = unique(PCa_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  PCa_CNV_DEG[[i]] <- CM_diff
}
library(fgsea)          
library(data.table)  
library(ggplot2)       
library(dplyr)          
library(msigdb)       
library(GSEABase)     
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)


for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )]){
  BC_CNV_DEG[[i]]$group <- "different"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC > 0 & BC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  BC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC < 0 & BC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  BC_CNV_DEG[[i]][pos,]$group <- "LOH"
}
BC_CNV_function_group <- list()
for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )] ){
  temp <- BC_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  BC_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}

for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )]){
  ccRCC_CNV_DEG[[i]]$group <- "different"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC > 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC < 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "LOH"
}

ccRCC_CNV_function_group <- list()
for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )] ){
  temp <- ccRCC_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  ccRCC_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}
for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )]){
  PCa_CNV_DEG[[i]]$group <- "different"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC > 0 & PCa_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  PCa_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC < 0 & PCa_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  PCa_CNV_DEG[[i]][pos,]$group <- "LOH"
}
PCa_CNV_function_group <- list()
for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )] ){
  temp <- PCa_CNV_DEG[[i]]
  temp_list <- list()
  for(z in unique(temp$chr)){
    temp2 <- temp[which(temp$chr == z),]
    cnv_list <- list()
    for(c in c("AMP","LOH")){
      temp3 <- temp2[which(temp2$group == c),]
      x = temp3[order(temp3$avg_log2FC,decreasing = T),]
      genelist <- structure(x$avg_log2FC,names=x$gene)
      if(length(genelist)){
        res_CM <- GSEA(genelist,TERM2GENE = genesets,eps = 0)
        res_CM_2 <- data.frame(res_CM)
        res_CM_2$log10pvalue <- -log10(res_CM_2$pvalue)
        res_CM_2$log10adjp <- -log10(res_CM_2$p.adjust)
        res_CM_2$count <- apply(res_CM_2, 1, function(x) length(strsplit(as.character(x["core_enrichment"]), "/")[[1]]))
        res_CM_2$count_ratio <- res_CM_2$count/res_CM_2$setSize
        cnv_list[[c]] <- res_CM
      }
    }
    temp_list[[z]] <- cnv_list
  }
  PCa_CNV_function_group[[i]] <- temp_list
  print(paste0(i,z,c))
}

PCa_CNV_DEG2 <- PCa_CNV_DEG
BC_CNV_DEG2 <- BC_CNV_DEG
ccRCC_CNV_DEG2 <- ccRCC_CNV_DEG

for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )] ){
  PCa_CNV_DEG[[i]] <- PCa_CNV_DEG[[i]][which(PCa_CNV_DEG[[i]]$group != "different"),]
  write.table(PCa_CNV_DEG[[i]],paste0("PCa_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )] ){
  BC_CNV_DEG[[i]] <- BC_CNV_DEG[[i]][which(BC_CNV_DEG[[i]]$group != "different"),]
  write.table(BC_CNV_DEG[[i]],paste0("BC_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )] ){
  ccRCC_CNV_DEG[[i]] <- ccRCC_CNV_DEG[[i]][which(ccRCC_CNV_DEG[[i]]$group != "different"),]
  write.table(ccRCC_CNV_DEG[[i]],paste0("ccRCC_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

BC_all_gene_AMP <- unique(c(BC_CNV_DEG[["Tumor cells_s1"]][which(BC_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s2"]][which(BC_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s3"]][which(BC_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"]))

BC_all_gene_LOH <- unique(c(BC_CNV_DEG[["Tumor cells_s1"]][which(BC_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s2"]][which(BC_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s3"]][which(BC_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"]))


PCa_all_gene_AMP <- unique(c(PCa_CNV_DEG[["Tumor cells_s1"]][which(PCa_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s2"]][which(PCa_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s3"]][which(PCa_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"]))

PCa_all_gene_LOH <- unique(c(PCa_CNV_DEG[["Tumor cells_s1"]][which(PCa_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s2"]][which(PCa_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s3"]][which(PCa_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"]))



ccRCC_all_gene_AMP <- unique(c(ccRCC_CNV_DEG[["Tumor cells_s1"]][which(ccRCC_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s2"]][which(ccRCC_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s3"]][which(ccRCC_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s4"]][which(ccRCC_CNV_DEG[["Tumor cells_s4"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s5"]][which(ccRCC_CNV_DEG[["Tumor cells_s5"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s6"]][which(ccRCC_CNV_DEG[["Tumor cells_s6"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s7"]][which(ccRCC_CNV_DEG[["Tumor cells_s7"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s8"]][which(ccRCC_CNV_DEG[["Tumor cells_s8"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s9"]][which(ccRCC_CNV_DEG[["Tumor cells_s9"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s10"]][which(ccRCC_CNV_DEG[["Tumor cells_s10"]]$group == "AMP"),"gene"]))

ccRCC_all_gene_LOH <- unique(c(ccRCC_CNV_DEG[["Tumor cells_s1"]][which(ccRCC_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s2"]][which(ccRCC_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s3"]][which(ccRCC_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s4"]][which(ccRCC_CNV_DEG[["Tumor cells_s4"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s5"]][which(ccRCC_CNV_DEG[["Tumor cells_s5"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s6"]][which(ccRCC_CNV_DEG[["Tumor cells_s6"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s7"]][which(ccRCC_CNV_DEG[["Tumor cells_s7"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s8"]][which(ccRCC_CNV_DEG[["Tumor cells_s8"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s9"]][which(ccRCC_CNV_DEG[["Tumor cells_s9"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s10"]][which(ccRCC_CNV_DEG[["Tumor cells_s10"]]$group == "LOH"),"gene"]))

gene_entrez <- bitr(
  intersect(intersect(BC_all_gene_AMP,PCa_all_gene_AMP), ccRCC_all_gene_AMP),
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
library(enrichplot)

dotplot(ego, showCategory = 10)

barplot(ego, showCategory = 10)

emapplot(pairwise_termsim(ego))

cnetplot(ego, showCategory = 10)

gene_entrez <- bitr(
  intersect(intersect(BC_all_gene_LOH,PCa_all_gene_LOH), ccRCC_all_gene_LOH),
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
library(enrichplot)
dotplot(ego, showCategory = 10)
barplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego))
cnetplot(ego, showCategory = 10)

#TCGA CNV -----
library(dplyr)
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Copy Number Variation", 
                  data.type = "Masked Copy Number Segment")
GDCdownload(query, method = "api", files.per.chunk = 100)
segment_dat <- GDCprepare(query = query)
data_CNV <- segment_dat
dim(data_CNV)
ov_cnv <- data_CNV[,-1]
ov_cnv <- ov_cnv[,c(6,1:5)]
tumor_seg <- ov_cnv[substr(ov_cnv$Sample,14,15)=="01",]
names(tumor_seg) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
dim(tumor_seg) ##336810      6
head(tumor_seg)

extracted_ids <- sub("^([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", tumor_seg$Sample)
load("D:/Urinary system tumors/work/2_SNP/ccRCC_TCGA_clinical.RData")

which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"])
segment_dat2 <- tumor_seg[which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"]),]

names(segment_dat2) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
write.table(segment_dat2,"MaskedCopyNumberSegment.txt",sep="\t",
            quote = F,col.names = F,row.names = F)


hg_marker_file <- read.delim("snp6.na35.remap.hg38.subset.txt.gz")
dim(hg_marker_file) ## 1837507       7
view(head(hg_marker_file))
Hg_marker_file <- hg_marker_file[hg_marker_file$freqcnv=="FALSE",]
head(Hg_marker_file)       
Hg_marker_file <- Hg_marker_file[,1:3]
names(Hg_marker_file) <- c("Marker Name","Chromosome","Marker Position")
head(Hg_marker_file)
write.table(Hg_marker_file,file = "Hg_marker.txt",sep = "\t",col.names = T,row.names = F)


library(dplyr)
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Copy Number Variation", 
                  data.type = "Masked Copy Number Segment")
GDCdownload(query, method = "api", files.per.chunk = 100)
segment_dat <- GDCprepare(query = query)
data_CNV <- segment_dat
dim(data_CNV)
ov_cnv <- data_CNV[,-1]
ov_cnv <- ov_cnv[,c(6,1:5)]
tumor_seg <- ov_cnv[substr(ov_cnv$Sample,14,15)=="01",]
names(tumor_seg) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
dim(tumor_seg) ##336810      6
head(tumor_seg)

extracted_ids <- sub("^([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", tumor_seg$Sample)
load("D:/Urinary system tumors/work/2_SNP/PCa_TCGA_clinical.RData")

which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"])
segment_dat2 <- tumor_seg[which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"]),]

names(segment_dat2) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
write.table(segment_dat2,"MaskedCopyNumberSegment.txt",sep="\t",
            quote = F,col.names = F,row.names = F)


query <- GDCquery(project = "TCGA-BLCA", 
                  data.category = "Copy Number Variation", 
                  data.type = "Masked Copy Number Segment")
GDCdownload(query, method = "api", files.per.chunk = 100)
segment_dat <- GDCprepare(query = query)
data_CNV <- segment_dat
dim(data_CNV)
ov_cnv <- data_CNV[,-1]
ov_cnv <- ov_cnv[,c(6,1:5)]
tumor_seg <- ov_cnv[substr(ov_cnv$Sample,14,15)=="01",]
names(tumor_seg) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
dim(tumor_seg) ##336810      6
head(tumor_seg)

extracted_ids <- sub("^([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", tumor_seg$Sample)
load("D:/Urinary system tumors/work/2_SNP/BC_TCGA_clinical.RData")

which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"])
segment_dat2 <- tumor_seg[which(extracted_ids %in% clinical[which(clinical$gender == "MALE"),"bcr_patient_barcode"]),]

names(segment_dat2) <- c("Sample","Chromosome","Start Position","End Position","Num markers","Seg.CN")
write.table(segment_dat2,"MaskedCopyNumberSegment.txt",sep="\t",
            quote = F,col.names = F,row.names = F)
library(maftools)
all.lesions <- "./PRAD/GISTIC_2.0/all_lesions.conf_90.txt"
amp.genes <- "./PRAD/GISTIC_2.0/amp_genes.conf_90.txt"
del.genes <- "./PRAD/GISTIC_2.0/del_genes.conf_90.txt"
scores.gis <- "./PRAD/GISTIC_2.0/scores.gistic"
coad.gistic = readGistic(gisticAllLesionsFile = all.lesions, 
                         gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, 
                         gisticScoresFile = scores.gis, isTCGA = TRUE)
pdf("PRAD_GISTIC_2.0_maftools.pdf",height = 6,width = 12)
gisticChromPlot(gistic = coad.gistic
                ,markBands = "all"
                ,ref.build = "hg38"
)
dev.off()
#GISTIC_2.0 plot
library(BSgenome.Hsapiens.UCSC.hg38) 
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), 
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
)
df$chromNum <- 1:length(df$chromName) 
df <- df[1:22,] 
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength))

df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])

tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle
scores <- read.table("./PRAD/GISTIC_2.0/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
chromID <- scores$Chromosome

scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
library(ggplot2)
library(ggsci)

df$ypos <- rep(c(0.2,0.25),11)

ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
ggsave("PRAD_GISTIC_2.0_scores.pdf",height = 4, width = 12)

#inferCNV chrX and chr Y ----
#ccRCC
ccRCC_CNV_gene <- read.table("./inferCNV_XY/out_result_cnv_ccRCC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                             sep = "\t",header = T)
ccRCC_CNV_gene<- ccRCC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                             "1"="Complete loss",
                                                             "2"="One copy_loss",
                                                             "3"="Neutral",
                                                             "4"="One copy_gain",
                                                             "5"="Two copies_gain",
                                                             "6"=">Two copies_gain"))
ccRCC_CNV_gene <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$CNV != "Neutral"),]

ccRCC_CNV_DEG <-list()
for( i in unique(ccRCC@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- ccRCC_CNV_gene[which(ccRCC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(ccRCC,
                         features = unique(ccRCC_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  ccRCC_CNV_DEG[[i]] <- CM_diff
}

#BC
BC_CNV_gene <- read.table("./inferCNV_XY/out_result_cnv_BC/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                          sep = "\t",header = T)
BC_CNV_gene<- BC_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                       "1"="Complete loss",
                                                       "2"="One copy_loss",
                                                       "3"="Neutral",
                                                       "4"="One copy_gain",
                                                       "5"="Two copies_gain",
                                                       "6"=">Two copies_gain"))
BC_CNV_gene <- BC_CNV_gene[which(BC_CNV_gene$CNV != "Neutral"),]

BC_CNV_DEG <-list()
for( i in unique(BC@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- BC_CNV_gene[which(BC_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(BC,
                         features = unique(BC_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  BC_CNV_DEG[[i]] <- CM_diff
}

#PCa
PCa_CNV_gene <- read.table("./inferCNV_XY/out_result_cnv_PCa/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat",
                           sep = "\t",header = T)
PCa_CNV_gene<- PCa_CNV_gene %>% dplyr::mutate(CNV=recode(state,
                                                         "1"="Complete loss",
                                                         "2"="One copy_loss",
                                                         "3"="Neutral",
                                                         "4"="One copy_gain",
                                                         "5"="Two copies_gain",
                                                         "6"=">Two copies_gain"))
PCa_CNV_gene <- PCa_CNV_gene[which(PCa_CNV_gene$CNV != "Neutral"),]

PCa_CNV_DEG <-list()
for( i in unique(PCa@meta.data$CNV_type)){
  z <- paste0(i,".",i,"_s1")
  temp <- PCa_CNV_gene[which(PCa_CNV_gene$cell_group_name == z),]
  CM_diff <- FindMarkers(PCa,
                         features = unique(PCa_CNV_gene$gene),
                         only.pos = FALSE,
                         logfc.threshold = 0,
                         min.pct = 0,
                         ident.1 = i,
                         ident.2 = "Normal",
                         group.by = "CNV_type")
  CM_diff$gene <- row.names(CM_diff)
  CM_diff <- CM_diff[which(CM_diff$gene %in% temp$gene),]
  CM_diff <- merge(CM_diff,temp[,c(1,4:8)], by="gene",all.x=T)
  CM_diff <- CM_diff[which(CM_diff$p_val_adj < 0.05),]
  PCa_CNV_DEG[[i]] <- CM_diff
}

library(fgsea)      
library(data.table)    
library(ggplot2)       
library(dplyr)      
library(msigdb)         
library(GSEABase)   
library(msigdbr)
library(GSVA)
library(clusterProfiler)
genesets_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
genesets_hallmark <- subset(genesets_hallmark,select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_GO <- msigdbr(species = "Homo sapiens", category = "C5")
genesets_GO <- subset(genesets_GO,gs_subcat%in%c("GO:BP","GO:CC","GO:MF"),select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets_KEGG <- msigdbr(species = "Homo sapiens", category = "C2")
genesets_KEGG <- subset(genesets_KEGG,gs_subcat=="CP:KEGG",select = c("gs_name","gene_symbol")) %>% as.data.frame()

genesets <- rbind(genesets_GO,genesets_hallmark)
genesets <- rbind(genesets,genesets_KEGG)


for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )]){
  BC_CNV_DEG[[i]]$group <- "different"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC > 0 & BC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  BC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(BC_CNV_DEG[[i]]$avg_log2FC < 0 & BC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  BC_CNV_DEG[[i]][pos,]$group <- "LOH"
}


for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )]){
  ccRCC_CNV_DEG[[i]]$group <- "different"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC > 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(ccRCC_CNV_DEG[[i]]$avg_log2FC < 0 & ccRCC_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  ccRCC_CNV_DEG[[i]][pos,]$group <- "LOH"
}

for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )]){
  PCa_CNV_DEG[[i]]$group <- "different"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC > 0 & PCa_CNV_DEG[[i]]$CNV %in% c("One copy_gain","Two copies_gain",">Two copies_gain") )
  PCa_CNV_DEG[[i]][pos,]$group <- "AMP"
  pos <- which(PCa_CNV_DEG[[i]]$avg_log2FC < 0 & PCa_CNV_DEG[[i]]$CNV %in% c("Complete loss","One copy_loss") )
  PCa_CNV_DEG[[i]][pos,]$group <- "LOH"
}

PCa_CNV_DEG2 <- PCa_CNV_DEG
BC_CNV_DEG2 <- BC_CNV_DEG
ccRCC_CNV_DEG2 <- ccRCC_CNV_DEG

for(i in names(PCa_CNV_DEG)[-which( names(PCa_CNV_DEG) == "Normal" )] ){
  PCa_CNV_DEG[[i]] <- PCa_CNV_DEG[[i]][which(PCa_CNV_DEG[[i]]$group != "different"),]
  write.table(PCa_CNV_DEG[[i]],paste0("PCa_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

for(i in names(BC_CNV_DEG)[-which( names(BC_CNV_DEG) == "Normal" )] ){
  BC_CNV_DEG[[i]] <- BC_CNV_DEG[[i]][which(BC_CNV_DEG[[i]]$group != "different"),]
  write.table(BC_CNV_DEG[[i]],paste0("BC_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

for(i in names(ccRCC_CNV_DEG)[-which( names(ccRCC_CNV_DEG) == "Normal" )] ){
  ccRCC_CNV_DEG[[i]] <- ccRCC_CNV_DEG[[i]][which(ccRCC_CNV_DEG[[i]]$group != "different"),]
  write.table(ccRCC_CNV_DEG[[i]],paste0("ccRCC_",i,"_CNV_gene.txt"),sep = "\t",quote = F)
}

BC_all_gene <- unique(c(BC_CNV_DEG[["Tumor cells_s1"]]$gene,
                        BC_CNV_DEG[["Tumor cells_s2"]]$gene,
                        BC_CNV_DEG[["Tumor cells_s3"]]$gene))

PCa_all_gene <- unique(c(PCa_CNV_DEG[["Tumor cells_s1"]]$gene,
                         PCa_CNV_DEG[["Tumor cells_s2"]]$gene,
                         PCa_CNV_DEG[["Tumor cells_s3"]]$gene))

ccRCC_all_gene <- unique(c(ccRCC_CNV_DEG[["Tumor cells_s1"]]$gene,ccRCC_CNV_DEG[["Tumor cells_s2"]]$gene,
                           ccRCC_CNV_DEG[["Tumor cells_s3"]]$gene,ccRCC_CNV_DEG[["Tumor cells_s4"]]$gene,
                           ccRCC_CNV_DEG[["Tumor cells_s5"]]$gene,ccRCC_CNV_DEG[["Tumor cells_s6"]]$gene,
                           ccRCC_CNV_DEG[["Tumor cells_s7"]]$gene,ccRCC_CNV_DEG[["Tumor cells_s8"]]$gene,
                           ccRCC_CNV_DEG[["Tumor cells_s9"]]$gene,ccRCC_CNV_DEG[["Tumor cells_s10"]]$gene))
intersect(intersect(BC_all_gene,PCa_all_gene), ccRCC_all_gene)
library(clusterProfiler)
library(org.Hs.eg.db)

chrX <- BC_all_gene
gene_entrez <- bitr(
  intersect(intersect(BC_all_gene,PCa_all_gene), ccRCC_all_gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",     # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
library(enrichplot)

dotplot(ego, showCategory = 10)

barplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego))
cnetplot(ego, showCategory = 10)


BC_all_gene_AMP <- unique(c(BC_CNV_DEG[["Tumor cells_s1"]][which(BC_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s2"]][which(BC_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s3"]][which(BC_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"]))

BC_all_gene_LOH <- unique(c(BC_CNV_DEG[["Tumor cells_s1"]][which(BC_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s2"]][which(BC_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                            BC_CNV_DEG[["Tumor cells_s3"]][which(BC_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"]))


PCa_all_gene_AMP <- unique(c(PCa_CNV_DEG[["Tumor cells_s1"]][which(PCa_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s2"]][which(PCa_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s3"]][which(PCa_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"]))

PCa_all_gene_LOH <- unique(c(PCa_CNV_DEG[["Tumor cells_s1"]][which(PCa_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s2"]][which(PCa_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                             PCa_CNV_DEG[["Tumor cells_s3"]][which(PCa_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"]))

ccRCC_all_gene_AMP <- unique(c(ccRCC_CNV_DEG[["Tumor cells_s1"]][which(ccRCC_CNV_DEG[["Tumor cells_s1"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s2"]][which(ccRCC_CNV_DEG[["Tumor cells_s2"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s3"]][which(ccRCC_CNV_DEG[["Tumor cells_s3"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s4"]][which(ccRCC_CNV_DEG[["Tumor cells_s4"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s5"]][which(ccRCC_CNV_DEG[["Tumor cells_s5"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s6"]][which(ccRCC_CNV_DEG[["Tumor cells_s6"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s7"]][which(ccRCC_CNV_DEG[["Tumor cells_s7"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s8"]][which(ccRCC_CNV_DEG[["Tumor cells_s8"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s9"]][which(ccRCC_CNV_DEG[["Tumor cells_s9"]]$group == "AMP"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s10"]][which(ccRCC_CNV_DEG[["Tumor cells_s10"]]$group == "AMP"),"gene"]))

ccRCC_all_gene_LOH <- unique(c(ccRCC_CNV_DEG[["Tumor cells_s1"]][which(ccRCC_CNV_DEG[["Tumor cells_s1"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s2"]][which(ccRCC_CNV_DEG[["Tumor cells_s2"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s3"]][which(ccRCC_CNV_DEG[["Tumor cells_s3"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s4"]][which(ccRCC_CNV_DEG[["Tumor cells_s4"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s5"]][which(ccRCC_CNV_DEG[["Tumor cells_s5"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s6"]][which(ccRCC_CNV_DEG[["Tumor cells_s6"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s7"]][which(ccRCC_CNV_DEG[["Tumor cells_s7"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s8"]][which(ccRCC_CNV_DEG[["Tumor cells_s8"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s9"]][which(ccRCC_CNV_DEG[["Tumor cells_s9"]]$group == "LOH"),"gene"],
                               ccRCC_CNV_DEG[["Tumor cells_s10"]][which(ccRCC_CNV_DEG[["Tumor cells_s10"]]$group == "LOH"),"gene"]))

gene <- c(intersect(intersect(BC_all_gene_AMP,PCa_all_gene_AMP), ccRCC_all_gene_AMP),
          intersect(intersect(BC_all_gene_LOH,PCa_all_gene_LOH), ccRCC_all_gene_LOH))


gene_entrez <- bitr(
  gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",     # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
dotplot(ego, showCategory = 10)
barplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego))
cnetplot(ego, showCategory = 10)


gene_entrez <- bitr(
  intersect(intersect(BC_all_gene_AMP,PCa_all_gene_AMP), ccRCC_all_gene_AMP),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",     # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
library(enrichplot)
dotplot(ego, showCategory = 10)
barplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego))
cnetplot(ego, showCategory = 10)

gene_entrez <- bitr(
  intersect(intersect(BC_all_gene_LOH,PCa_all_gene_LOH), ccRCC_all_gene_LOH),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
ego <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",     # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE      
)
library(enrichplot)
dotplot(ego, showCategory = 10)
barplot(ego, showCategory = 10)
emapplot(pairwise_termsim(ego))
cnetplot(ego, showCategory = 10)

CNV_gene <- read.delim("clipboard", header = T) #table s2 PCa_CNV_DEG_chrX_Y
CNV_gene <- na.omit(CNV_gene)
Time_gene <- read.table("D:/Urinary system tumors/work/3_Pseutime/1_tumor_epi/PCa/2_PCa_tumor_modulated_genes.txt", header = T,sep = "\t")

Time_CNV_gene <- Time_gene[unique(CNV_gene$gene),]

TF_gene <- read.csv("D:/PCa_ATAC/PCa_motif_peak_gene.csv", header = T)
TF_gene <-distinct(TF_gene[,c(7,10)])
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
                     TF_gene$geneName)),"PCa_PPI.txt",quote = F,sep = "\t",row.names = F)

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
write.table(TF_gene,"PCa_TF_gene.txt",quote = F,sep = "\t",row.names = F)
table <- cbind(unique(TF_gene$Target),"gene")
table <- rbind(table,cbind(unique(TF_gene$Source),"TF"))
write.table(table,"PCa_table.txt",quote = F,sep = "\t",row.names = F)


CNV_gene <- read.delim("clipboard", header = T)#table s2 ccRCC_CNV_DEG_chrX_Y
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


write.table(unique(c(TF_gene$Target,TF_gene$Source)),"ccRCC_PPI.txt",quote = F,sep = "\t",row.names = F)
colnames(TF_gene) <- c("Target","Source")
TF_gene$Direction <- "Direction"
write.table(TF_gene,"ccRCC_TF_gene.txt",quote = F,sep = "\t",row.names = F)
table <- cbind(unique(TF_gene$Target),"gene")
table <- rbind(table,cbind(unique(TF_gene$Source),"TF"))
write.table(table,"ccRCC_table.txt",quote = F,sep = "\t",row.names = F)