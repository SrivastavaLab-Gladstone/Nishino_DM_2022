## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)

## load librarires 
library(Seurat)
library(WGCNA)
library(lmerTest)
library(data.table)
library(enrichR)

## set up WGCNA output directory
setwd("/path/to/desired/project/src_folder/")
WGCNAbasedir <- "path/to/results/WGCNA/"
dir.create(WGCNAbasedir)

## load Seurat object
WGCNA_SeuratObj <- readRDS("path/to/SeuratObject.RDS")
WGCNA_SeuratObj

## review loaded Seurat object
# DimPlot(WGCNA_SeuratObj, label = TRUE)
# DimPlot(WGCNA_SeuratObj, split.by = "gem.group")
# table(WGCNA_SeuratObj@meta.data$gem.group)
# table(Idents(WGCNA_SeuratObj))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## save variable features determinaed via SCT
varFeatures <- rownames(WGCNA_SeuratObj@assays$SCT@scale.data)
levels(Idents(WGCNA_SeuratObj))
WGCNA_SeuratObj_expData <- t(as.matrix(GetAssayData(WGCNA_SeuratObj)))
WGCNA_SeuratObj_expData_sub <- WGCNA_SeuratObj_expData[,varFeatures] ## only use variable genes in analysis


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Run WGCNA
allowWGCNAThreads(nThreads = 12) ## (`nThreads = NULL` == automatically determined)

net <- blockwiseModules(WGCNA_SeuratObj_expData_sub, power = 10,
                        corType = "bicor", # use robust correlation
                        networkType = "signed", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.15,
                        numericLabels = FALSE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

blockwiseModules()

saveRDS(net, file = paste(WGCNAbasedir, "RDSname.RDS", sep = ""))

## or simply load the precomputed `net` object:
# net <- readRDS(paste(WGCNAbasedir, "RDSname.RDS"))
# table(net$colors)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Convert labels to colors for plotting
mergedColors = net$colors
## Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## save module color:gene dictionary
# levels(as.factor(mergedColors))
colorinfo <- net$colors
write.csv(colorinfo, paste(WGCNAbasedir, "gene2module_dict.csv"), sep = "")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## add module scores for all WGCNA module of interest, e.g. turquoise and blue
ModuleGenes_turquoise <- colnames(WGCNA_SeuratObj_expData)[net$colors == "turquoise"]
ModuleGenes_blue <- colnames(WGCNA_SeuratObj_expData)[net$colors == "blue"]
## gene lists must be elements of a list of character vectors
ModuleScoreGenes <- list(ModuleGenes_turquoise, ModuleGenes_blue)
names(ModuleScoreGenes) <- c("WGCNA_turquoise", "WGCNA_blue")
## add module scores
WGCNA_SeuratObj <- AddModuleScore(WGCNA_SeuratObj, 
                                  features = ModuleScoreGenes, 
                                  name = names(ModuleScoreGenes))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## GO analysis
enrichrOut <- enrichr(ModuleGenes, databases = c("GO_Biological_Process_2021")) #, "GO_Molecular_Function_2021", "GO_Cellular_Component_2021"))
plotEnrich(enrichrOut$GO_Biological_Process_2021, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "WGCNA Signature: GO_Biological_Process_2021")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
### LMERTEST
## pull metdata w/ module scores
data4lmer <- WGCNA_SeuratObj@meta.data
head(data4lmer)
## manual type coercion for lmerTest
data4lmer$ConditionToTest <- as.factor(data4lmer$ConditionToTest) # e.g. WT, KO1, and KO2
data4lmer$random_effect_group <- as.factor(data4lmer$random_effect_group)
data4lmer$data4lmer$WGNCAmodule_color1 <- as.numeric(data4lmer$WGNCAmodule_color1)
## run LMER test
lmer_output <- lmer(formula = WGNCAmodule_color1~ConditionToTest + (1|random_effect_group), data = data4lmer)
summary(lmer_output)
TestSummary <- summary(lmer_output)
TestSummary <- data.frame(TestSummary$coefficients)
## pull p-value for condition of interest - edit as needed
raw_pval <- TestSummary$Pr...t..[rownames(TestSummary) %like% "KO2"]
## n_comparisons = no. of groups you're comparing
n_comparisons <- length(levels(data4lmer$ConditionToTest))
print("###################################")
print("Adj p-val, FDR")
print(stats::p.adjust(p = raw_pval, method = "fdr", n = n_comparisons))
print("###################################")




### sessioninfo
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS 12.2
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] enrichR_3.0           data.table_1.14.0     lmerTest_3.1-3        lme4_1.1-26           Matrix_1.3-2          WGCNA_1.70-3          fastcluster_1.1.25    dynamicTreeCut_1.63-1
# [9] Seurat_4.0.1.9006     SeuratObject_4.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.2.1       Hmisc_4.5-0           plyr_1.8.6            igraph_1.3.0          lazyeval_0.2.2        splines_4.0.5         listenv_0.8.0         scattermore_0.7      
# [9] ggplot2_3.3.6         digest_0.6.29         foreach_1.5.2         htmltools_0.5.2       GO.db_3.12.1          fansi_1.0.3           magrittr_2.0.3        checkmate_2.0.0      
# [17] memoise_2.0.0         tensor_1.5            cluster_2.1.1         doParallel_1.0.17     ROCR_1.0-11           globals_0.15.1        matrixStats_0.58.0    spatstat.sparse_2.0-0
# [25] jpeg_0.1-8.1          colorspace_2.0-3      blob_1.2.1            ggrepel_0.9.1         xfun_0.30             dplyr_1.0.8           crayon_1.5.1          jsonlite_1.8.0       
# [33] spatstat.data_2.1-0   impute_1.64.0         survival_3.2-10       zoo_1.8-9             iterators_1.0.14      glue_1.6.2            polyclip_1.10-0       gtable_0.3.0         
# [41] leiden_0.3.7          future.apply_1.9.0    BiocGenerics_0.36.1   abind_1.4-5           scales_1.1.1          DBI_1.1.1             miniUI_0.1.1.1        Rcpp_1.0.8.3         
# [49] viridisLite_0.4.0     xtable_1.8-4          htmlTable_2.1.0       reticulate_1.24       spatstat.core_2.0-0   foreign_0.8-81        bit_4.0.4             preprocessCore_1.52.1
# [57] Formula_1.2-4         stats4_4.0.5          htmlwidgets_1.5.3     httr_1.4.2            RColorBrewer_1.1-3    ellipsis_0.3.2        ica_1.0-2             pkgconfig_2.0.3      
# [65] nnet_7.3-15           uwot_0.1.10           deldir_1.0-6          utf8_1.2.2            tidyselect_1.1.2      rlang_1.0.2           reshape2_1.4.4        later_1.1.0.1        
# [73] AnnotationDbi_1.52.0  munsell_0.5.0         tools_4.0.5           cachem_1.0.4          cli_3.2.0             generics_0.1.3        RSQLite_2.2.12        ggridges_0.5.3       
# [81] stringr_1.4.0         fastmap_1.1.0         goftest_1.2-2         knitr_1.32            bit64_4.0.5           fitdistrplus_1.1-3    purrr_0.3.4           RANN_2.6.1           
# [89] pbapply_1.5-0         future_1.26.1         nlme_3.1-152          mime_0.10             compiler_4.0.5        rstudioapi_0.13       curl_4.3              plotly_4.9.3         
# [97] png_0.1-7             spatstat.utils_2.3-0  tibble_3.1.6          statmod_1.4.35        stringi_1.7.6         lattice_0.20-41       nloptr_1.2.2.2        vctrs_0.4.0          
# [105] pillar_1.7.0          lifecycle_1.0.1       spatstat.geom_2.4-0   lmtest_0.9-38         RcppAnnoy_0.0.18      cowplot_1.1.1         irlba_2.3.5           httpuv_1.5.5         
# [113] patchwork_1.1.1       R6_2.5.1              latticeExtra_0.6-29   promises_1.2.0.1      KernSmooth_2.23-18    gridExtra_2.3         IRanges_2.24.1        parallelly_1.32.0    
# [121] codetools_0.2-18      boot_1.3-27           MASS_7.3-53.1         assertthat_0.2.1      rjson_0.2.20          sctransform_0.3.2     S4Vectors_0.28.1      mgcv_1.8-34          
# [129] parallel_4.0.5        grid_4.0.5            rpart_4.1-15          minqa_1.2.4           tidyr_1.2.0           Rtsne_0.15            numDeriv_2016.8-1.1   Biobase_2.50.0       
# [137] shiny_1.6.0           base64enc_0.1-3      
