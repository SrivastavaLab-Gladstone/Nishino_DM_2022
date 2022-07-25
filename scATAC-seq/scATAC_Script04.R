#### Script04; Subset analyses

## Load libraries
library(ArchR)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
addArchRThreads(threads = 28) 
set.seed(1)

## Re-load ArchR project!
E105_DM_Harmony_trial　<- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony_trial_Latest")

### Subset Meso+CM scATAC-seq data and merge it with Meso+CM subset scRNA-seq data 
## Subset Meso+CM scATAC-seq data
idxSample <- BiocGenerics::which(E105_DM_Harmony_trial$Clusters %in% c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C24", "C25"))
cellsSample <- E105_DM_Harmony_trial$cellNames[idxSample]
Meso <- E105_DM_Harmony_trial[cellsSample, ]

## Preprocess Meso+CM scRNA-seq subset data
scRNA <- readRDS("./scRNA_data/04_labeled_Meso_Subset_07-02-2021.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 4000) 
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)

## Integration (Unconstrained, w/o imputation) 
Meso <- addGeneIntegrationMatrix(
  ArchRProj = Meso, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Meso",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE, 
  groupRNA = "clusters",
  nameCell = "predictedCell_Meso",
  nameGroup = "predictedGroup_Meso",
  nameScore = "predictedScore_Meso",
  nGenes = 4000,
  useImputation = T,
  force = T
)
pal <- paletteDiscrete(values = seRNA$clusters)
p1 <- plotEmbedding(
  Meso, 
  colorBy = "cellColData",
  name = "predictedGroup_Meso",
  embedding = "UMAP",
  pal = pal)
p1
plotPDF(p1, name = "Plot-UMAP-RNA-Harmony_Meso-Integration_imputation-T.pdf", ArchRProj = Meso, addDOC = FALSE, width = 5, height = 5)
Meso <- addGeneIntegrationMatrix(
  ArchRProj = Meso,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Meso",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell_Meso",
  nameGroup = "predictedGroup_Meso",
  nameScore = "predictedScore_Meso",
  nGenes = 4000,
  useImputation = F,
  force = T
)
pal <- paletteDiscrete(values = seRNA$clusters)
p1 <- plotEmbedding(
  Meso,
  colorBy = "cellColData",
  name = "predictedGroup_Meso",
  embedding = "UMAP",
  pal = pal)
p1
plotPDF(p1, name = "Plot-UMAP-RNA-Harmony_Meso-Integration_imputation-F.pdf", ArchRProj = Meso, addDOC = FALSE, width = 5, height = 5)
# Save
saveArchRProject(ArchRProj = Meso, outputDirectory = "/path/to/ArchrOutputDir/Meso", load = TRUE)
Meso 

## My manual labeling --> make sure this is a separate metadata column from existing clusters
temp_cluster_names <- as.character(revalue(Meso$Clusters, c(
  "C10"="AHF1",          
  "C4"="Non_Meso",            
  "C24"="EpiC",         
  "C18"="Non_Meso",            
  "C19"="Non_Meso",
  "C14"="ParaxialMeso1",          
  "C11"="pSHF",            
  "C22"="Non_Meso",     
  "C25"="Cardiomyocyte",           
  "C20"="Non_Meso",            
  "C5"="Non_Meso",        
  "C15"="AHF2",            
  "C23"="Non_Meso",    
  "C3"="Non_Meso",           
  "C8"="Non_Meso",            
  "C16"="ParaxialMeso2",          
  "C12"="PharyngealMeso",           
  "C13"="BrM",     
  "C9"="Non_Meso" ,     
  "C1"="Non_Meso",          
  "C6"="Non_Meso" ,     
  "C17"="Non_Meso",            
  "C21"="Non_Meso",     
  "C7"="Non_Meso",           
  "C2"="Non_Meso"            
)))
length(unique(temp_cluster_names))
Meso$Clusters3 <- temp_cluster_names
Meso@cellColData@listData[["Clusters3"]] <- temp_cluster_names
table(Meso$Clusters3,Meso$Sample)

p1 <- plotEmbedding(Meso, colorBy = "cellColData", name = "Clusters3")
p1
plotPDF(p1,
        name = "Plot-UMAP-Harmony_Meso_Manual-Annotations_celltype.pdf",
        ArchRProj = Meso,
        addDOC = FALSE, width = , height = 5)

## Calculate jaccard index score; subset Meso+CM from scATAC-entire data and integrated scRNA-seq Meso+CM data
Meso_ji <- Reduce(bind_rows,lapply(unique(Meso$Clusters), function(TN_label) {
  df <- Reduce(bind_rows,lapply(unique(Meso$predictedGroup_Meso), function(CCA_label) {
    ji <- sum(Meso$Clusters == TN_label & Meso$predictedGroup_Meso == CCA_label)/sum(Meso$Clusters == TN_label | Meso$predictedGroup_Meso ==CCA_label)
    df <- data.frame(CCA_cell_type=CCA_label, JI=ji, stringsAsFactors = F)
    return(df)
  }))
  df$TN_label_cell_type <- TN_label
  return(df)
}))
Meso_ji_mat <- spread(Meso_ji, key = TN_label_cell_type, value = JI)
row.names(Meso_ji_mat) <- Meso_ji_mat$CCA_cell_type
Meso_ji_mat <- Meso_ji_mat[,2:ncol(Meso_ji_mat)]
summary(apply(Meso_ji_mat, 1, max))
summary(apply(Meso_ji_mat, 2, max))
pheatmap::pheatmap(Meso_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100))
pheatmap::pheatmap(Meso_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")
i <- apply(Meso_ji_mat, 1, max) >= 0.15
j <- apply(Meso_ji_mat, 2, max) >= 0.15
pheatmap::pheatmap(Meso_ji_mat[i,j], clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")


### NC compartment analysis 
## Load libraries
library(ArchR)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(tidyr)
library(RColorBrewer)
addArchRThreads(threads = 28) 
set.seed(1)

## Re-load ArchR project!
E105_DM_Harmony_trial　<- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony_trial_Latest")

## Subset NC scATAC-seq data
idxSample <- BiocGenerics::which(E105_DM_Harmony_trial$Clusters %in% c("C8", "C17", "C18", "C19", "C20"))
cellsSample <- E105_DM_Harmony_trial$cellNames[idxSample]
NC <- E105_DM_Harmony_trial[cellsSample, ]

## Preprocess NC subset scRNA-seq data 
scRNA <- readRDS("./scRNA_data/NC_labeled_07-17-2021.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 4000) 
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)
seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)

## Integration (Unconstrained, w/o imputation) 
NC <- addGeneIntegrationMatrix(
  ArchRProj = NC,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_NC",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell_NC",
  nameGroup = "predictedGroup_NC",
  nameScore = "predictedScore_NC",
  nGenes = 4000,
  useImputation = F,
  force = T
)
pal <- paletteDiscrete(values = seRNA$clusters)
p1 <- plotEmbedding(
  NC,
  colorBy = "cellColData",
  name = "predictedGroup_NC",
  embedding = "UMAP",
  pal = pal)
p1
plotPDF(p1, name = "Plot-UMAP-RNA-Harmony_NC-Integration_imputation-F.pdf", ArchRProj = NC, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = NC, outputDirectory = "/path/to/ArchrOutputDir/NC", load = TRUE)
NC 

## Calculate jaccard index score; subset NC from scATAC-entire data and integrated scRNA-seq NC data
NC_ji <- Reduce(bind_rows,lapply(unique(NC$Clusters), function(TN_label) {
  df <- Reduce(bind_rows,lapply(unique(NC$predictedGroup_NC), function(CCA_label) {
    ji <- sum(NC$Clusters == TN_label & NC$predictedGroup_NC == CCA_label)/sum(NC$Clusters == TN_label | NC$predictedGroup_NC ==CCA_label)
    df <- data.frame(CCA_cell_type=CCA_label, JI=ji, stringsAsFactors = F)
    return(df)
  }))
  df$TN_label_cell_type <- TN_label
  return(df)
}))
NC_ji_mat <- spread(NC_ji, key = TN_label_cell_type, value = JI)
row.names(NC_ji_mat) <- NC_ji_mat$CCA_cell_type
NC_ji_mat <- NC_ji_mat[,2:ncol(NC_ji_mat)]
summary(apply(NC_ji_mat, 1, max))
summary(apply(NC_ji_mat, 2, max))
pheatmap::pheatmap(NC_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100))
pheatmap::pheatmap(NC_ji_mat, clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")
i <- apply(NC_ji_mat, 1, max) >= 0.15
j <- apply(NC_ji_mat, 2, max) >= 0.15
pheatmap::pheatmap(NC_ji_mat[i,j], clustering_method = "ward.D2", color = colorRampPalette((brewer.pal(n = 9, name ="Reds")))(100), display_numbers = T, number_color = "White")

## original UMAP for NC
p1 <- plotEmbedding(
  NC,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP",
  pal = pal)
p1
plotPDF(p1, name = "Plot-UMAP-RNA-Harmony_NC-originalClusters.pdf", ArchRProj = NC, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = NC, outputDirectory = "/path/to/ArchrOutputDir/NC", load = TRUE)
NC 


## Session Information
sessionInfo()
# R version 4.0.4 (2021-02-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] gridExtra_2.3                      BSgenome.Mmusculus.UCSC.mm10_1.4.0 BSgenome_1.58.0                    rtracklayer_1.50.0                
# [5] Biostrings_2.58.0                  XVector_0.30.0                     circlize_0.4.13.1001               ComplexHeatmap_2.6.2              
# [9] nabor_0.5.0                        gtable_0.3.0                       SeuratWrappers_0.3.0               dplyr_1.0.5                       
# [13] ArchR_1.0.1                        magrittr_2.0.1                     rhdf5_2.34.0                       Matrix_1.3-2                      
# [17] data.table_1.14.0                  SummarizedExperiment_1.20.0        Biobase_2.50.0                     GenomicRanges_1.42.0              
# [21] GenomeInfoDb_1.26.4                IRanges_2.24.1                     S4Vectors_0.28.1                   BiocGenerics_0.36.0               
# [25] MatrixGenerics_1.2.1               matrixStats_0.58.0                 ggplot2_3.3.3                     
# 
# loaded via a namespace (and not attached):
#   [1] plyr_1.8.6               igraph_1.2.6             lazyeval_0.2.2           splines_4.0.4            BiocParallel_1.24.1      listenv_0.8.0           
# [7] scattermore_0.7          digest_0.6.27            htmltools_0.5.1.1        magick_2.7.0             fansi_0.4.2              tensor_1.5              
# [13] cluster_2.1.1            ROCR_1.0-11              remotes_2.2.0            globals_0.14.0           spatstat.sparse_1.2-1    colorspace_2.0-0        
# [19] ggrepel_0.9.1            xfun_0.22                crayon_1.4.1             RCurl_1.98-1.2           jsonlite_1.7.2           spatstat.data_2.0-0     
# [25] survival_3.2-7           zoo_1.8-9                glue_1.4.2               polyclip_1.10-0          zlibbioc_1.36.0          leiden_0.3.7            
# [31] GetoptLong_1.0.5         DelayedArray_0.16.2      Rhdf5lib_1.12.1          future.apply_1.7.0       shape_1.4.5              abind_1.4-5             
# [37] scales_1.1.1             pheatmap_1.0.12          DBI_1.1.1                miniUI_0.1.1.1           Rcpp_1.0.6               viridisLite_0.3.0       
# [43] xtable_1.8-4             clue_0.3-58              reticulate_1.18          spatstat.core_1.65-5     rsvd_1.0.3               htmlwidgets_1.5.3       
# [49] httr_1.4.2               RColorBrewer_1.1-2       ellipsis_0.3.1           Seurat_4.0.2             ica_1.0-2                XML_3.99-0.5            
# [55] pkgconfig_2.0.3          farver_2.1.0             uwot_0.1.10              deldir_0.2-10            utf8_1.2.1               tidyselect_1.1.0        
# [61] labeling_0.4.2           rlang_0.4.10             reshape2_1.4.4           later_1.1.0.1            munsell_0.5.0            tools_4.0.4             
# [67] generics_0.1.0           ggridges_0.5.3           stringr_1.4.0            fastmap_1.1.0            goftest_1.2-2            fitdistrplus_1.1-3      
# [73] purrr_0.3.4              RANN_2.6.1               pbapply_1.4-3            future_1.21.0            nlme_3.1-152             mime_0.10               
# [79] rstudioapi_0.13          compiler_4.0.4           plotly_4.9.3             png_0.1-7                spatstat.utils_2.0-0     tibble_3.1.0            
# [85] stringi_1.5.3            lattice_0.20-41          vctrs_0.3.6              pillar_1.5.1             lifecycle_1.0.0          rhdf5filters_1.2.0      
# [91] BiocManager_1.30.10      spatstat.geom_1.65-5     lmtest_0.9-38            GlobalOptions_0.1.2      RcppAnnoy_0.0.18         cowplot_1.1.1           
# [97] bitops_1.0-6             irlba_2.3.3              httpuv_1.5.5             patchwork_1.1.1          R6_2.5.0                 promises_1.2.0.1        
# [103] KernSmooth_2.23-18       parallelly_1.24.0        codetools_0.2-18         MASS_7.3-53.1            gtools_3.8.2             assertthat_0.2.1        
# [109] rjson_0.2.20             withr_2.4.1              SeuratObject_4.0.0       GenomicAlignments_1.26.0 Rsamtools_2.6.0          sctransform_0.3.2       
# [115] GenomeInfoDbData_1.2.4   mgcv_1.8-34              rpart_4.1-15             tidyr_1.1.3              Cairo_1.5-12.2           Rtsne_0.15              
# [121] shiny_1.6.0              tinytex_0.30  