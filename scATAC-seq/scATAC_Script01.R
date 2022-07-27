### Script01: prepare ArchR project
### Load libraries and data ####
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 24) 
set.seed(1)


### Import data ####
inputFiles <- c("./E105fragments/Vehicle01_fragments.tsv.gz",
                "./E105fragments/Vehicle02_fragments.tsv.gz",
                "./E105fragments/Vehicle03_fragments.tsv.gz",
                "./E105fragments/STZ01_fragments.tsv.gz",
                "./E105fragments/STZ02_fragments.tsv.gz",
                "./E105fragments/STZ03_fragments.tsv.gz")
names(inputFiles) <- c("Veh01", "Veh02", "Veh03", "STZ01", "STZ02", "STZ03")
inputFiles
## add reference
addArchRGenome("mm10")
## create arrrow files (with modifications)
starttime <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 7,
  minFrags = 19952, 
  maxFrags = 1e+6,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles

## infer doubelts
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

## create ArchR project
E105_DM <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "E105_DM",
  copyArrows = TRUE 
)
E105_DM

## save ArchRProject
saveArchRProject(ArchRProj = E105_DM, outputDirectory = "/path/to/ArchrOutputDir/E105_DM", load = TRUE)


### restart R, delete fragments files & the initial arrow files, and reload your project ####
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 24) 
set.seed(1)
E105_DM <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM")
## Plot QC matrix for the ArchRProject
df <- getCellColData(E105_DM, select = c("log10(nFrags)", "TSSEnrichment"))
df
p<- ggPoint(x=df[,1], y=df[,2], colorDensity = TRUE, continuousSet = "sambaNight",
            xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment",
            xlim = c(log10(500), quantile(df[,1], probs=0.99)),
            ylim = c(0, quantile(df[,2], probs=0.99))) +
  geom_hline(yintercept = 7, lty="dashed") +
  geom_vline(xintercept = 4.3, lty="dashed")
p 
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = E105_DM, addDOC = FALSE)


## plot sample statistics ####
## TSS score ridge plot
p1 <- plotGroups(ArchRProj=E105_DM, groupBy="Sample", colorBy="cellColData",
                 name="TSSEnrichment", plotAs="ridges")
p1
## TSS score violin plot
p2 <- plotGroups(ArchRProj=E105_DM, groupBy="Sample", colorBy="cellColData",
                 name="TSSEnrichment", plotAs="violin", alpha=0.4, addBoxPlot=TRUE)
p2
## log10(unique nuclear fragments) ridge plot
p3 <- plotGroups(ArchRProj=E105_DM, groupBy="Sample", colorBy="cellColData",
                 name="log10(nFrags)", plotAs="ridges")
p3
## log10(unique nuclear fragments) violin plot
p4 <- plotGroups(ArchRProj=E105_DM, groupBy="Sample", colorBy="cellColData",
                 name="log10(nFrags)", plotAs="violin", alpha=0.4, addBoxPlot=TRUE)
p4

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = E105_DM, addDOC = FALSE, width = 4, height = 4)
rm(p, p1,p2,p3,p4)

## plot sample fragment size distribution and TSS enrichment profiles
## fragment size distributions
p1 <- plotFragmentSizes(ArchRProj = E105_DM)
p1
## TSS enrichment profiles
p2 <- plotTSSEnrichment(ArchRProj = E105_DM)
p2
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = E105_DM, addDOC = FALSE, width = 5, height = 5)
rm(p1,p2)

## save
saveArchRProject(ArchRProj = E105_DM, outputDirectory = "/path/to/ArchrOutputDir/E105_DM", load = TRUE)


### Clustering (with batch correction) ####
## restart R and reload your project
library(ArchR)
library(dplyr)
library(SeuratWrappers)
library(pheatmap)
addArchRThreads(threads = 24)
set.seed(1)
E105_DM <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM")
E105_DM_Harmony <- filterDoublets(ArchRProj = E105_DM)
E105_DM_Harmony
rm(E105_DM)
E105_DM_Harmony <- addIterativeLSI(ArchRProj = E105_DM_Harmony, useMatrix = "TileMatrix", name = "IterativeLSI",
                                   iterations = 2, clusterParams = list(resolution=c(2), sampleCells=10000, maxClusters = 6, n.start=10),
                                   varFeatures = 25000, dimsToUse = 1:30 , LSIMethod = 2, seed = 1, force = TRUE)
## add harmony batch correction
E105_DM_Harmony <- addHarmony(ArchRProj = E105_DM_Harmony, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
E105_DM_Harmony <- addClusters(input = E105_DM_Harmony, reducedDims = "Harmony", method = "Seurat",
                       name = "Clusters", resolution = 1, force = TRUE, maxClusters = 25)
table(E105_DM_Harmony$Clusters)
cM <- confusionMatrix(paste0(E105_DM_Harmony$Clusters), paste0(E105_DM_Harmony$Sample))
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p
plotPDF(p, name = "Plot-ConfusionMatrix-Heatmap_Harmony.pdf", ArchRProj = E105_DM_Harmony, addDOC = FALSE, width = 5, height = 5)

## visualize 2D UMAP embedding
E105_DM_Harmony <- addUMAP(ArchRProj = E105_DM_Harmony, reducedDims = "Harmony", name = "UMAP",
                   nNeighbors = 40, minDist = 0.4, metric = "cosine" , force = TRUE)
## UMAP colored by “Sample”:
p1 <- plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
## UMAP colored by “Clusters”:
p2 <- plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
## save plots 
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters_Harmony.pdf",
        ArchRProj = E105_DM_Harmony, addDOC = FALSE, width = 5, height = 5)
## save the project
saveArchRProject(ArchRProj = E105_DM_Harmony, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony", load = TRUE)



### scRNA-seq integration and cluster ID using non-batch corrected clustering ####
## restart R and reload your project
library(ArchR)
library(dplyr)
library(SeuratWrappers)
addArchRThreads(threads = 24) 
set.seed(1)
E105_DM_Harmony <- loadArchRProject(path = "/path/to/ArchrOutputDir/E105_DM_Harmony")

### cluster identification ####
## identifying Marker Genes
## to identify marker genes based on gene scores, we call the getMarkerFeatures() function with useMatrix = "GeneScoreMatrix"
markersGS <- getMarkerFeatures(ArchRProj = E105_DM_Harmony, useMatrix = "GeneScoreMatrix", 
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon", maxCells = 500)
## extract the marker genes we got from the above step. You can access to each cluster data using "$", like $C1.
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")
numv <- length(markerList@listData)
for (i in 1:numv) {
  foundmarkers <- markerList@listData[i]
  write.csv(foundmarkers, file = paste0("./E105_DM_Harmony/MarkerGS/C", i, "_MarkerGS.csv"))
}
## plot heatmap of makerGS with some labels of marker genes
markerGenes <- c("Hoxb1", "Osr1", "Tbx5", "Fgf8", "Isl1", "Mef2c", "Dlx5", "Twist1", "Barx1","Foxd1","Sox2","Tlx1",
                 "Tbx18", "Wt1", "Tcf21", "Irx4", "Myh7", "Rgs5", "Tnnt2", "Rspo3", "Emcn", "Pecam1","Alx3","Gata4",
                 "Epcam", "Bmp2", "Tbx20", "Ttn", "Foxa2", "Foxi2", "Sox10","Sox7", "Kdr", "Gsta1", "Hand2","Gata6",
                 "Tfap2b", "Tubb3", "Neurod4", "Lyz2", "Runx1", "Coro1a", "Foxe1", "Nkx2-1", "Hhex", "Isl1","Tbx1",
                 "Tcf15", "Pax3","Myf5", "Hba-x", "Hba-a1", "Klf1", "Afp", "Ttr", "Apoa2", "Shox2", "Meox1",
                 "Six1", "Six2", "Mef2c", "Fgb", "Rbp4","Dlx2", "Postn")
heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25", 
                               labelMarkers = markerGenes, transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = E105_DM_Harmony, addDOC = FALSE)

## visualize marker genes on an UMAP embedding
p <- plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "GeneScoreMatrix", name = markerGenes, 
                   embedding = "UMAP",quantCut = c(0.01, 0.95),imputeWeights = NULL)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = E105_DM_Harmony, addDOC = FALSE, width = 5, height = 5)


### marker gene imuputation with MAGIC ####
## We can use MAGIC to impute gene scores by smoothing signal across nearby cells. 
# this greatly improves the visual interpretation of gene scores.
E105_DM_Harmony <- addImputeWeights(E105_DM_Harmony)
p <- plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "GeneScoreMatrix", name = markerGenes, 
                   embedding = "UMAP",imputeWeights = getImputeWeights(E105_DM_Harmony))
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = E105_DM_Harmony, addDOC = FALSE, width = 5, height = 5)
rm(p, p2)
## save
saveArchRProject(ArchRProj = E105_DM_Harmony, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony", load = TRUE)


### maker peaks ####
## peak calls using "Cluster" 
E105_DM_Harmony <- addGroupCoverages(ArchRProj = E105_DM_Harmony, groupBy = "Clusters", minCells = 40, maxCells = 500, force = T)
E105_DM_Harmony <- addReproduciblePeakSet(ArchRProj = E105_DM_Harmony, groupBy = "Clusters", pathToMacs2 = "/usr/local/bin/macs2",
                                          peaksPerCell = 500, maxPeaks = 150000, minCells = 25, force = T)
getPeakSet(E105_DM_Harmony)
merged_peak_setGR <- getPeakSet(E105_DM_Harmony)
E105_DM_Harmony2 <- addPeakMatrix(E105_DM_Harmony)
getAvailableMatrices(E105_DM_Harmony2)
markersPeaks <- getMarkerFeatures(ArchRProj = E105_DM_Harmony2, useMatrix = "PeakMatrix", groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks3 <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.000000000001 & Log2FC >= 2",
  transpose = TRUE, labelMarkers = NULL
) 
ncol(heatmapPeaks3@matrix) # 28982 -> less than 32K -> fine
ComplexHeatmap::draw(heatmapPeaks3, heatmap_legend_side = "bot", annotation_legend_side = "bot")

numv <- length(markerList@listData)
for (i in 1:numv) {
  foundmarkers <- markerList@listData[i]
  write.csv(foundmarkers, file = paste0("./E105_DM_Harmony/MarkerPeaks/C", i, "_MarkerPeaks.csv"))
}

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
#numv <- length(markerList@listData)
numv <- c(1,3:25)
for (i in numv) {
  foundmarkers <- markerList@listData[i]
  gr <- foundmarkers[[1]]
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr),
                   scores=c(rep("", length(gr))),
                   strands=strand(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/MarkerPeaks/MarkerPeak_C",i,".bed"), quote=F, sep="\t", row.names=F, col.names=F)
  # the following csv is for running GREAT analysis
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr),
                   names=1:length(gr))
  write.table(df, file=paste0("/path/to/ArchrOutputDir/E105_DM_Harmony/MarkerPeaks/MarkerPeak_C",i,"_GREAT.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
}
saveArchRProject(ArchRProj = E105_DM_Harmony, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony", load = TRUE)

### motif enrichment in marker peaks (heatmaps) ####
## based on Cluster
plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
table(E105_DM_Harmony$Clusters)
markersPeaks <- getMarkerFeatures(ArchRProj = E105_DM_Harmony,
                                  useMatrix = "PeakMatrix",
                                  groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  testMethod = "wilcoxon"
                                  )
E105_DM_Harmony <- addMotifAnnotations(ArchRProj = E105_DM_Harmony, motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,
                                   ArchRProj = E105_DM_Harmony,
                                   peakAnnotation = "Motif",
                                   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
                                   )
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_cisbp", width = 16, height = 8, ArchRProj = E105_DM_Harmony, addDOC = FALSE)

## motif enrichment in marker peaks cont'd
## based on labels by "clusters + condition" stored in "Clusters_Cd"
temp_cluster_names <- paste0(E105_DM_Harmony$Clusters,"_",substr(E105_DM_Harmony$Sample, 1,3))
length(unique(temp_cluster_names))
E105_DM_Harmony$Clusters_Cd <- temp_cluster_names
E105_DM_Harmony@cellColData@listData[["Clusters_Cd"]] <- temp_cluster_names
table(E105_DM_Harmony$Clusters,E105_DM_Harmony$Sample)
table(E105_DM_Harmony$Clusters_Cd,E105_DM_Harmony$Sample)

plotEmbedding(ArchRProj = E105_DM_Harmony, colorBy = "cellColData", name = "Clusters_Cd", embedding = "UMAP")
table(E105_DM_Harmony$Clusters_Cd)
markersPeaks <- getMarkerFeatures(ArchRProj = E105_DM_Harmony,
                                  useMatrix = "PeakMatrix",
                                  groupBy = "Clusters_Cd",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  testMethod = "wilcoxon"
)
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks,
                                   ArchRProj = E105_DM_Harmony,
                                   peakAnnotation = "Motif",
                                   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_Cluster-condition_cisbp", width = 16, height = 10, ArchRProj = E105_DM_Harmony, addDOC = FALSE)
## save
saveArchRProject(ArchRProj = E105_DM_Harmony, outputDirectory = "/path/to/ArchrOutputDir/E105_DM_Harmony", load = TRUE)


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