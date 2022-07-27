### Script07: Filtering the candidate enhancers from detected DARs in AHF2(C15)
## The aim of this script is to identify the candidate enhancers with deferential accessibility between conditions and with an association with DEGs in the same cluster. 
## Strategy 1. Extract the peak2GeneLinks list which includes peak locations and associated genes for each of peaks 
##          2. Extract the DAR list and merge it with the peak2GeneLinks list using peak/DAR location data 
##          3. Take a DEG list 
##          4. Merge the DEG list with the DAR-peak2GeneLinks list.

### Load libraries and data ####
E105_DM_Harmony2 <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony2")
E105_DM_Harmony2 <- addPeak2GeneLinks(ArchRProj = E105_DM_Harmony2, reducedDims = "Harmony")

p2geneDF <- metadata(E105_DM_Harmony2@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]

markerTest_15 <- getMarkerFeatures(ArchRProj = E105_DM_Harmony2, useMatrix = "PeakMatrix", groupBy = "Clusters_Cd",
                                   testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"),
                                   useGroups = paste0("C15_STZ"), bgdGroups = paste0("C15_Veh"))


### Analyze DARs ####
## More accessible DARs 
markerList_C15_Open <- getMarkers(markerTest_15, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList_C15_Open_DF <-  markerList_C15_Open$C15_STZ
markerList_C15_Open_DF$peakName <- paste0(markerList_C15_Open_DF$seqnames, "_", markerList_C15_Open_DF$start, "_", markerList_C15_Open_DF$end)
markerList_C15_Open_DF_2 <- as.data.frame(markerList_C15_Open_DF)
p2geneDF_2 <- as.data.frame(p2geneDF)
DAR_Open_P2g_integrated <- inner_join(markerList_C15_Open_DF_2, p2geneDF_2, by="peakName")

## Less accessible DARs 
markerList_C15_Close <- getMarkers(markerTest_15, cutOff = "FDR <= 0.05 & Log2FC <= -1")
markerList_C15_Close_DF <-  markerList_C15_Close$C15_STZ
markerList_C15_Close_DF$peakName <- paste0(markerList_C15_Close_DF$seqnames, "_", markerList_C15_Close_DF$start, "_", markerList_C15_Close_DF$end)
markerList_C15_Close_DF_2 <- as.data.frame(markerList_C15_Close_DF)
DAR_Close_P2g_integrated <- inner_join(markerList_C15_Close_DF_2, p2geneDF_2, by="peakName")

## Load the DEG lists from AHF2 cluster
AHF2_DEG <- read.csv(file = "./DEGlists/*AHF2_DE_STZ_vs_Vehicle.csv")
AHF2_DEG_UP <- filter(AHF2_DEG, p_val_adj <= 0.05 & avg_log2FC > 0)
AHF2_DEG_DOWN <- filter(AHF2_DEG, p_val_adj <= 0.05 & avg_log2FC < 0)
names(AHF2_DEG_UP)[1] <- "geneName"
names(AHF2_DEG_DOWN)[1] <- "geneName"
colnames(AHF2_DEG_UP)

## For Open DARs and Upregulated DEG
DEGUP_DAR_Open_P2g_integrated <- inner_join(AHF2_DEG_UP, DAR_Open_P2g_integrated, by="geneName")
write.csv(DEGUP_DAR_Open_P2g_integrated, "DEGUP_DAR_Open_P2g_integrated.csv")

## For Closed DARs and Downregulated DEG
DEGDOWN_DAR_Close_P2g_integrated <- inner_join(AHF2_DEG_DOWN, DAR_Close_P2g_integrated, by="geneName")
write.csv(DEGDOWN_DAR_Close_P2g_integrated, "DEGDOWN_DAR_Close_P2g_integrated.csv")

## For Open DARs and Downregulated DEG
DEGDOWN_DAR_Open_P2g_integrated <- inner_join(AHF2_DEG_DOWN, DAR_Open_P2g_integrated, by="geneName")
write.csv(DEGDOWN_DAR_Open_P2g_integrated, "DEGDOWN_DAR_Open_P2g_integrated.csv")

## For Close DARs and Upregulated DEG
DEGUP_DAR_Close_P2g_integrated <- inner_join(AHF2_DEG_UP, DAR_Close_P2g_integrated, by="geneName")
write.csv(DEGUP_DAR_Close_P2g_integrated, "DEGUP_DAR_Close_P2g_integrated.csv")

### used the same framework for DAR-DEG analyses in PA2NC (C19).


### sessioninfo ####
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