#### Script02_Identification of DARs with Motif enrichment analyses

# load libraries and set seed
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 24) 
set.seed(1)

# Re-load ArchR project!
E105_DM_Harmony <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony")

### Step1; Peak call per cluster+condition(Veh or STZ)
# Add labels of "clusters + condition" information and store them in "Clusters_Cd"
temp_cluster_names <- paste0(E105_DM_Harmony$Clusters,"_",substr(E105_DM_Harmony$Sample, 1,3))
length(unique(temp_cluster_names))
E105_DM_Harmony$Clusters_Cd <- temp_cluster_names
E105_DM_Harmony@cellColData@listData[["Clusters_Cd"]] <- temp_cluster_names
table(E105_DM_Harmony$Clusters,E105_DM_Harmony$Sample)
table(E105_DM_Harmony$Clusters_Cd,E105_DM_Harmony$Sample)

### peak calls using "Clusters_Cd" and get maker peaks
E105_DM_Harmony <- addGroupCoverages(ArchRProj = E105_DM_Harmony, groupBy = "Clusters_Cd", minCells = 40, maxCells = 500, force = T)
E105_DM_Harmony <- addReproduciblePeakSet(ArchRProj = E105_DM_Harmony, groupBy = "Clusters_Cd", pathToMacs2 = "/usr/local/bin/macs2",
                                          peaksPerCell = 500, maxPeaks = 150000, minCells = 25, force = T)
getPeakSet(E105_DM_Harmony)
merged_peak_setGR <- getPeakSet(E105_DM_Harmony)

E105_DM_Harmony2 <- addPeakMatrix(E105_DM_Harmony)
getAvailableMatrices(E105_DM_Harmony2)
saveArchRProject(ArchRProj = E105_DM_Harmony2, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony2", load = TRUE)
E105_DM_Harmony2 <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony2")

markersPeaks <- getMarkerFeatures(ArchRProj = E105_DM_Harmony2, useMatrix = "PeakMatrix", groupBy = "Clusters_Cd",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
# record maker peaks
numv <- length(markerList@listData)
for (i in c(1:2,5:17, 19:40, 42:numv)) {
  foundmarkers <- markerList@listData[i]
  gr <- foundmarkers[[1]]
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr),
                   scores=c(rep("", length(gr))),
                   strands=strand(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/MarkerPeaks/per_ClusterCon/MarkerPeak_",names(markerList@listData)[i],".bed"), quote=F, sep="\t", row.names=F, col.names=F)
  # save csv files to be used for running GREAT analyses
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/MarkerPeaks/per_ClusterCon/MarkerPeak_",names(markerList@listData)[i],"_GREAT.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
}
# save
saveArchRProject(ArchRProj = E105_DM_Harmony, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony", load = TRUE)

# Plot a HeatMap of marker peaks
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE, labelMarkers = NULL
) 
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

### Step2; Identification of differential accessible regions (DARs) between conditions and motif enrichment analysis in the detected DARs
E105_DM_Harmony2 <- addMotifAnnotations(ArchRProj = E105_DM_Harmony2, motifSet = "homer", name = "Motif", force = T)
celltype <- unique(E105_DM_Harmony2$Clusters)
for (i in celltype) {print(i)}
for (i in celltype) {
  markerTest <- getMarkerFeatures(
    ArchRProj = E105_DM_Harmony2, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters_Cd",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = paste0(i,"_STZ"),
    bgdGroups = paste0(i,"_Veh")
  )
  pma <- plotMarkers(seMarker = markerTest, name = paste0(i,"_STZ"), cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
  pv <- plotMarkers(seMarker = markerTest, name = paste0(i,"_STZ"), cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
  plotPDF(pma, pv, name = paste0(i,"_STZ-vs-Veh-DCA-MA-Volcano_Harmony"), width = 5, height = 5, ArchRProj = E105_DM_Harmony2, addDOC = FALSE)

  # enriched motifs in more accessible in STZ
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = E105_DM_Harmony2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )
  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  # enriched motifs in less accessible in STZ
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = E105_DM_Harmony2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC <= -1"
  )
  motifsDo
  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))
  ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  plotPDF(ggUp, ggDo, name = paste0(i,"_STZ-vs-Veh-DCA-Motifs-Enriched_Harmony"), width = 5, height = 5, ArchRProj = E105_DM_Harmony2, addDOC = FALSE)
}
saveArchRProject(ArchRProj = E105_DM_Harmony2, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony2", load = TRUE)

# store detected DARs
celltype <- c("C1", "C5", "C11", "C14", "C15", "C16", "C19", "C20","C22", "C23") #"C1", "C5", "C11"
for (i in celltype) {
  markerTest <- getMarkerFeatures(
    ArchRProj = E105_DM_Harmony2, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters_Cd",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = paste0(i,"_STZ"),
    bgdGroups = paste0(i,"_Veh")
  )
  pma <- plotMarkers(seMarker = markerTest, name = paste0(i,"_STZ"), cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
  pv <- plotMarkers(seMarker = markerTest, name = paste0(i,"_STZ"), cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
  ggAlignPlots(pma, pv, type = "h")
  
  DAR_UP <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
  DAR_UP <- DAR_UP[[1]]
  gr <- DAR_UP
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr),
                   scores=c(rep("", length(gr))),
                   strands=strand(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/DAR_clusters/DAR_",i,"_STZvsVeh_UP.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/DAR_clusters/DAR_",i,"_STZvsVeh_GREAT_UP.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
  DAR_DOWN <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & Log2FC <= -1", returnGR = TRUE)
  DAR_DOWN <- DAR_DOWN[[1]]
  gr <- DAR_DOWN
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr),
                   scores=c(rep("", length(gr))),
                   strands=strand(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/DAR_clusters/DAR_",i,"_STZvsVeh_DOWN.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/DAR_clusters/DAR_",i,"_STZvsVeh_GREAT_DOWN.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}


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