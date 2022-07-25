### Script03_Analyses of Mesoderm-CM subset

### load libraries and set seed
library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
set.seed(7)

### load RDS from Script02
E105 <- readRDS(file ="../data/rds/02_entire_broad-labeled_6-28-2021.RDS")

### Subset Mesoderm and CM cells
Idents(E105) <- "SCT_snn_res.0.6"
levels(Idents(E105))
clusters2keep <- c(3,4,6,7,8,10,11,12,13,14,20,23)
Meso <- subset(E105, idents = clusters2keep)

### Run SCTransform to renormalize data
Meso <- SCTransform(Meso, assay = "RNA", new.assay.name = "SCT", variable.features.n = 3000, variable.features.rv.th = 1.3, 
                    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),return.only.var.genes = TRUE)
### Run PCA
Meso <- RunPCA(Meso, verbose = FALSE, npcs = 100)
### Run Harmony
Meso <- harmony::RunHarmony(Meso, group.by.vars = "gem.group", assay.use="SCT")
### Elbowplot, just to check
ElbowPlot(Meso, ndims = 50)
dims <- 1:30
Meso <- RunUMAP(Meso, reduction = "harmony", dims = dims)
Meso <- FindNeighbors(Meso, reduction = "harmony", dims = dims)
Meso <- FindClusters(Meso, resolution = c(0.6,1.4), random.seed = 7)
### Plot UMAPs and other QC parameters
Idents(Meso) <- "SCT_snn_res.0.6"
DimPlot(Meso, reduction = "umap", label = TRUE)
DimPlot(Meso, reduction = "umap", split.by = "condition", label = TRUE)
VlnPlot(Meso, features = c("nFeature_RNA"), pt.size = 0.1,ncol = 1)
VlnPlot(Meso, features = c("nCount_RNA"), pt.size = 0.1,ncol = 1)
VlnPlot(Meso, features = c("percent.mt"), pt.size = 0.1,ncol = 1)
Idents(Meso) <- "SCT_snn_res.1.4"
DimPlot(Meso, reduction = "umap", label = TRUE)
DimPlot(Meso, reduction = "umap", split.by = "condition", label = TRUE)
VlnPlot(Meso, features = c("nFeature_RNA"), pt.size = 0.1,ncol = 1)
VlnPlot(Meso, features = c("nCount_RNA"), pt.size = 0.1,ncol = 1)
VlnPlot(Meso, features = c("percent.mt"), pt.size = 0.1,ncol = 1)
### save RDS
saveRDS(Meso, file = "../data/rds/Mesoderm_CM_Epi_R-cellcycle_unlabeled_06-29-2021.RDS")

### Extract marker genes to do annotation
Idents(Meso) <- "SCT_snn_res.0.6"
allmarkers_1 <- FindAllMarkers(Meso, logfc.threshold = 0.5, only.pos = TRUE, min.pct = 0.5, return.thresh = 1e-20)
Idents(Meso) <- "SCT_snn_res.1.4"
allmarkers_2 <- FindAllMarkers(Meso, logfc.threshold = 0.5, only.pos = TRUE, min.pct = 0.5, return.thresh = 1e-20)

### Plot UMAP
DimPlot(Meso, label = TRUE) + NoLegend()
### Extract maker genes and plot a heatmap
allmarkers <- FindAllMarkers(Meso, logfc.threshold = 0.3, only.pos = TRUE, min.pct = 0.3, return.thresh = 1e-6)
write.csv(allmarkers, file = "../results/cluster_markers/Meso_07-02-2021.csv")
top5 <- allmarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
maxcells  <- min(table(Idents(Meso)))
DoHeatmap(subset(Meso, downsample = maxcells), features = top5$gene, size = 4)+ 
  theme(axis.text.y = element_text(size = 12))

### Differential gene expression analysis between conditions per each sub cell type
identlist <- levels(Idents(Meso))
print(identlist)
numvec <- as.vector(1:length(identlist))
for(i in numvec) {
  ident2sub <- identlist[i]
  print(ident2sub)
  sub_seuratobj <- subset(Meso, idents = ident2sub)
  sub_seuratobj <- SCTransform(sub_seuratobj, verbose = T, return.only.var.genes = FALSE) # moved vars.to.regress = c("percent.ribo", "percent.mt") to first SCT on whole obj
  saveRDS(sub_seuratobj, file = paste("../data/rds/2021-07-02_Mesosubset/",ident2sub, "_SCTnormSubset_top3000.RDS", sep = ""))
  Idents(sub_seuratobj) <- "condition"
  maxcells <- as.numeric(min(table(Idents(sub_seuratobj))))
  # Do DE analysis between condition
  foundmarkers <- FindMarkers(sub_seuratobj, ident.1 = "STZ", ident.2 = "Vehicle", logfc.threshold = 0.125, verbose = FALSE, max.cells.per.ident = maxcells)
  write.csv(foundmarkers, file = paste("../results/Mesoderm_CM_Epi_R-cellcycle/DEresults/2021-07-02/rawdata/",ident2sub, "_DE_STZ_vs_Vehicle.csv", sep = ""))
}




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