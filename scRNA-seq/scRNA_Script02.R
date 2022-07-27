### Script02: QC checks, filtering, and normalization
### Load libraries and data ####
library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
set.seed(7)
E105 <- readRDS(file = "../data/rds/01_rawObject_noFilters_scoresAdded_06-28-2021.RDS")

### Set QC cutoffs
mito_cutoffs <- c(1,6)
ribo_cutoffs <- c(10, 30)
nFeat_cutoffs <- c(4000, 9500)
nCount_cutoffs <- c(10000, 100000)

### Subset out low QC cells and cluster ####
E105  <- subset(E105, subset = percent.mt > min(mito_cutoffs) & percent.mt < max(mito_cutoffs) &
                  percent.ribo > min(ribo_cutoffs) & percent.ribo <  max(ribo_cutoffs) &
                  nFeature_RNA > min(nFeat_cutoffs) & nFeature_RNA <  max(nFeat_cutoffs) &
                  nCount_RNA > min(nCount_cutoffs) & nCount_RNA < max(nCount_cutoffs))

E105 <- SCTransform(E105, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000, 
                    variable.features.rv.th = 1.3, vars.to.regress = c("percent.mt","percent.ribo"), 
                    return.only.var.genes = TRUE)
E105 <- RunPCA(E105, verbose = FALSE, npcs = 100)
E105 <- harmony::RunHarmony(E105, group.by.vars = "gem.group", assay.use="SCT")
# ElbowPlot(E105, ndims = 50)
dims <- 1:32 
E105 <- RunUMAP(E105, reduction = "harmony", dims = dims) 
E105 <- FindNeighbors(E105, reduction = "harmony", dims = dims)
E105 <- FindClusters(E105, resolution = 0.6)
### Plot UMAPs
DimPlot(E105, label = TRUE)
DimPlot(E105, split.by = "gem.group")
### Plot QC parameters per cluster
VlnPlot(E105, features = "nFeature_RNA")
VlnPlot(E105, features = "nCount_RNA")
VlnPlot(E105, features = "percent.mt")


### Find marker genes and annotate celltypes #### 
## Find cluster markers
allmarkers <- FindAllMarkers(E105, logfc.threshold = 0.4, only.pos = TRUE, min.pct = 0.45, return.thresh = 1e-10)
## Add annotations based on the maker genes above
E105 <- RenameIdents(E105,
                     `0` = "NeuralCrest",
                     `1` = "NeuralCrest",
                     `2` = "Endoderm",
                     `3` = "Meso/CardiopharyngealProgenitor",
                     `4` = "Meso/CardiopharyngealProgenitor",
                     `5` = "Endoderm",
                     `6` = "Meso/CardiopharyngealProgenitor",
                     `7` = "Meso/CardiopharyngealProgenitor",
                     `8` = "Cardiomyocyte",
                     `9` = "Blood",
                     `10` = "Meso/CardiopharyngealProgenitor",
                     `11` = "Meso/CardiopharyngealProgenitor",
                     `12` = "Epicardium",
                     `13` = "Cardiomyocyte",
                     `14` = "Meso/CardiopharyngealProgenitor",
                     `15` = "Endothelium",
                     `16` = "NeuralCrest",
                     `17` = "Endoderm",
                     `18` = "Endothelium",
                     `19` = "Endothelium",
                     `20` = "Cardiomyocyte",
                     `21` = "NeuralCrest",
                     `22` = "NeuralCrest",
                     `23` = "Meso/CardiopharyngealProgenitor",
                     `24` = "Endoderm",
                     `25` = "Endoderm",
                     `26` = "Blood"
)

E105@active.ident <- factor(E105@active.ident,
                            levels=c("Cardiomyocyte","Meso/CardiopharyngealProgenitor","Epicardium","NeuralCrest","Endothelium","Endoderm",
                                     "Blood"))

## save RDS
saveRDS(E105, "../data/rds/02_entire_broad-labeled_08-04-2021.RDS")
## recalculate marker genes if desired & plot heatmap of top N marker genes per cluster
allmarkers <- FindAllMarkers(E105, logfc.threshold = 0.5, only.pos = TRUE, min.pct = 0.5, return.thresh = 1e-20)
top5 <- allmarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
## capture lowest cluster cell count to subset heatmap groups
maxcells  <- min(table(Idents(E105)))
## plot
DoHeatmap(subset(E105, downsample = maxcells), features = top5$gene, size = 6, raster = F) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + 
  guides(color=FALSE) + 
  theme(axis.text.y = element_text(size = 16))


### sessioninfo ####
# R version 4.0.4 (2021-02-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.6.6
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# other attached packages:
#   [1] biomaRt_2.46.3              RColorBrewer_1.1-2          enrichR_3.0                 ArchR_1.0.1                
# [5] magrittr_2.0.1              rhdf5_2.34.0                Matrix_1.3-2                data.table_1.14.0          
# [9] SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.4        
# [13] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.0         MatrixGenerics_1.2.1       
# [17] matrixStats_0.58.0          dplyr_1.0.5                 stringr_1.4.0               gridExtra_2.3              
# [21] ggplot2_3.3.5               SeuratObject_4.0.0          Seurat_4.0.2               
# loaded via a namespace (and not attached):
#   [1] circlize_0.4.13.1001   BiocFileCache_1.14.0   plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2        
# [6] splines_4.0.4          listenv_0.8.0          scattermore_0.7        digest_0.6.27          htmltools_0.5.2       
# [11] fansi_0.4.2            memoise_2.0.0          tensor_1.5             cluster_2.1.1          ROCR_1.0-11           
# [16] ComplexHeatmap_2.6.2   globals_0.14.0         askpass_1.1            spatstat.sparse_1.2-1  prettyunits_1.1.1     
# [21] colorspace_2.0-0       rappdirs_0.3.3         blob_1.2.1             ggrepel_0.9.1          xfun_0.31             
# [26] crayon_1.4.1           RCurl_1.98-1.2         jsonlite_1.7.2         spatstat.data_2.0-0    survival_3.2-7        
# [31] zoo_1.8-9              glue_1.4.2             polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.36.0       
# [36] XVector_0.30.0         leiden_0.3.7           GetoptLong_1.0.5       DelayedArray_0.16.2    Rhdf5lib_1.12.1       
# [41] future.apply_1.7.0     shape_1.4.5            abind_1.4-5            scales_1.1.1           DBI_1.1.1             
# [46] miniUI_0.1.1.1         Rcpp_1.0.6             progress_1.2.2         viridisLite_0.3.0      xtable_1.8-4          
# [51] clue_0.3-58            reticulate_1.18        spatstat.core_1.65-5   bit_4.0.4              htmlwidgets_1.5.3     
# [56] httr_1.4.2             ellipsis_0.3.1         ica_1.0-2              XML_3.99-0.5           pkgconfig_2.0.3       
# [61] dbplyr_2.1.0           uwot_0.1.10            deldir_0.2-10          utf8_1.2.1             AnnotationDbi_1.52.0  
# [66] tidyselect_1.1.0       rlang_0.4.10           reshape2_1.4.4         later_1.1.0.1          cachem_1.0.4          
# [71] munsell_0.5.0          tools_4.0.4            RSQLite_2.2.4          generics_0.1.0         ggridges_0.5.3        
# [76] evaluate_0.14          fastmap_1.1.0          yaml_2.2.1             goftest_1.2-2          bit64_4.0.5           
# [81] knitr_1.31             fitdistrplus_1.1-3     purrr_0.3.4            RANN_2.6.1             pbapply_1.4-3         
# [86] future_1.21.0          nlme_3.1-152           mime_0.10              xml2_1.3.2             compiler_4.0.4        
# [91] curl_4.3               plotly_4.9.3           png_0.1-7              spatstat.utils_2.0-0   tibble_3.1.0          
# [96] stringi_1.5.3          lattice_0.20-41        vctrs_0.3.6            pillar_1.5.1           lifecycle_1.0.0       
# [101] rhdf5filters_1.2.0     spatstat.geom_1.65-5   lmtest_0.9-38          GlobalOptions_0.1.2    RcppAnnoy_0.0.18      
# [106] cowplot_1.1.1          bitops_1.0-6           irlba_2.3.3            httpuv_1.5.5           patchwork_1.1.1       
# [111] R6_2.5.0               promises_1.2.0.1       KernSmooth_2.23-18     parallelly_1.24.0      codetools_0.2-18      
# [116] MASS_7.3-53.1          assertthat_0.2.1       openssl_1.4.3          rjson_0.2.20           withr_2.4.1           
# [121] sctransform_0.3.2      harmony_1.0            GenomeInfoDbData_1.2.4 hms_1.0.0              mgcv_1.8-34           
# [126] grid_4.0.4             rpart_4.1-15           tidyr_1.1.3            rmarkdown_2.14         Cairo_1.5-12.2        
# [131] Rtsne_0.15             shiny_1.6.0       