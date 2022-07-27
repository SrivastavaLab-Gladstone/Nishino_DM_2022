#### Script03; Integration of scRNA-seq data of corresponding samples to scATAC
### Load libraries and data ####
library(ArchR)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(plyr)
library(dplyr)
addArchRThreads(threads = 27) 
set.seed(1)
E105_DM_Harmony <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony")
E105_DM_Harmony_trial <- E105_DM_Harmony
rm(E105_DM_Harmony)


### Preprocess of scRNA-seq data before integration ####
scRNA <- readRDS("./scRNA_data/02_entire_broad-labeled_08-04-2021.RDS")
scRNA <- FindVariableFeatures(scRNA, assay = "SCT", nfeatures = 4000) 
Idents(scRNA) <- "celltype01" 
scRNA$clusters <- scRNA@active.ident
table(scRNA@active.ident)

seRNA <- as.SingleCellExperiment(scRNA, assay = "SCT")
colnames(colData(seRNA))
table(colData(seRNA)$clusters)
##  integration (Unconstrained Integration, without imputation) 
E105_DM_Harmony_trial <- addGeneIntegrationMatrix(
  ArchRProj = E105_DM_Harmony_trial,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupRNA = "clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  nGenes = 4000,
  useImputation = F,
  force = T
)
pal <- paletteDiscrete(values = seRNA$clusters)
p1 <- plotEmbedding(
  E105_DM_Harmony_trial,
  colorBy = "cellColData",
  name = "predictedGroup",
  embedding = "UMAP",
  pal = pal)
p1
plotPDF(p1, name = "Plot-UMAP-RNA-Harmony_Low_Latest-Integration_imputation-F.pdf", ArchRProj = E105_DM_Harmony_trial, addDOC = FALSE, width = 5, height = 5)
## save
saveArchRProject(ArchRProj = E105_DM_Harmony_trial, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony_trial_Latest", load = TRUE)
## Re-load ArchR project!
E105_DM_Harmony_trial <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony_trial_Latest")


### annotate scATAC-seq clusters based on the scRNA-seq intergration ####
## labeling with cell-type in "Clusters2"
cM <- confusionMatrix(E105_DM_Harmony_trial$Clusters, E105_DM_Harmony_trial$predictedGroup)
cM
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
cbind(labelOld, labelNew)

## manual labeling --> make sure this is a separate metadata column from existing clusters
temp_cluster_names <- as.character(revalue(E105_DM_Harmony_trial$Clusters, c(
  "C10"="Mesoderm_CardiacProgenitor",          
  "C4"="Endothelium",            
  "C24"="Epicardium",         
  "C18"="NeuralCrest",            
  "C19"="NeuralCrest",
  "C14"="Mesoderm_CardiacProgenitor",          
  "C11"="Mesoderm_CardiacProgenitor",            
  "C22"="Endoderm",     
  "C25"="Cardiomyocyte",           
  "C20"="NeuralCrest",            
  "C5"="Endothelium",        
  "C15"="Mesoderm_CardiacProgenitor",            
  "C23"="Endoderm",    
  "C3"="Blood",           
  "C8"="NeuralCrest",            
  "C16"="Mesoderm_CardiacProgenitor",          
  "C12"="Mesoderm_CardiacProgenitor",           
  "C13"="Mesoderm_CardiacProgenitor",     
  "C9"="Endoderm" ,     
  "C1"="Blood",          
  "C6"="Endothelium" ,     
  "C17"="NeuralCrest",            
  "C21"="Endoderm",     
  "C7"="Unknown2",           
  "C2"="Unknown1"            
)))
length(unique(temp_cluster_names))
E105_DM_Harmony_trial$Clusters2 <- temp_cluster_names
E105_DM_Harmony_trial@cellColData@listData[["Clusters2"]] <- temp_cluster_names
table(E105_DM_Harmony_trial$Clusters2,E105_DM_Harmony_trial$Sample)

p1 <- plotEmbedding(E105_DM_Harmony_trial, colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1,
        name = "Plot-UMAP-Harmony_Low_Latest-Manual-Annotations_celltype.pdf",
        ArchRProj = E105_DM_Harmony_trial,
        addDOC = FALSE, width = , height = 5)


### sessionfo ####
#                            STZ01 STZ02 STZ03 Veh01 Veh02 Veh03
# Blood                        155   256   240   368   330  1468
# Cardiomyocyte                400   499   427   976  1078  1355
# Endoderm                     943  1077   930   820  1157  1754
# Endothelium                  274   380   327   674   678   545
# Epicardium                    43   123   100   320   360   578
# Mesoderm_CardiacProgenitor  1996  2912  1804  1346  2081  2663
# NeuralCrest                  934  1211  1030  1849  1799  3078
# Unknown1                       2     4     3     3     1     4
# Unknown2                      35    57    61    35    27    40

# ## 2) labeling with cell-type + condition in "Clusters3"
# temp_cluster_names <- paste0(E105_DM_Harmony_trial$Clusters2,"_",substr(E105_DM_Harmony_trial$Sample, 1,3))
# length(unique(temp_cluster_names))
# E105_DM_Harmony_trial$Clusters3 <- temp_cluster_names
# E105_DM_Harmony_trial@cellColData@listData[["Clusters3"]] <- temp_cluster_names
# table(E105_DM_Harmony_trial$Clusters3,E105_DM_Harmony_trial$Sample)
# 
# p1 <- plotEmbedding(E105_DM_Harmony_trial, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP", pal = pal)
# p1
# plotPDF(p1,
#         name = "Plot-UMAP-Harmony_Low_Latest-Manual-Annotations_celltype.pdf",
#         ArchRProj = E105_DM_Harmony_trial,
#         addDOC = FALSE, width = , height = 5)
# 
# pal2 <- paletteDiscrete(values = E105_DM_Harmony_trial$Clusters3)
# pal2
# p2 <- plotEmbedding(E105_DM_Harmony_trial, colorBy = "cellColData", name = "Clusters3", embedding = "UMAP", pal = pal2)
# p2
# plotPDF(p2,
#         name = "Plot-UMAP-Harmony_Low_Latest-Manual-Annotations_celltype-condition.pdf",
#         ArchRProj = E105_DM_Harmony_trial,
#         addDOC = FALSE, width = , height = 5)

saveArchRProject(ArchRProj = E105_DM_Harmony_trial, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony_trial_Latest", load = TRUE)
E105_DM_Harmony_trial


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