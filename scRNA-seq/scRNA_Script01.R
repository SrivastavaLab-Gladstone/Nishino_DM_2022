### Script01_Data import and initiate Seurat object with metadata

### Import library
library(dplyr)
library(Seurat)
library(biomaRt)
## define "marts"
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

### load dataset - Tomo STZ+Vehicle E10.5 data (prespliced-mapped data)
data.E105 <- Read10X(data.dir = "../data/filtered_feature_bc_matrix/") # this is your cellranger aggr output
# Initialize the Seurat object with the raw (non-normalized data).
E105 <- CreateSeuratObject(counts = data.E105, project = "E10.5_STZvVeh", min.cells = 3, min.features = 200)
rm(data.E105)

### Add metadata
## annotate the aggregation csv and use that to add metadata.
AggrSheet <- read.csv("../data/aggregation.csv") # this is your csv from cellranger-aggr
AggrSheet$gem.group <- c("Veh01", 
                         "Veh02", 
                         "Veh03", 
                         "STZ01", 
                         "STZ02", 
                         "STZ03")
AggrSheet$condition <- c("Vehicle", 
                         "Vehicle", 
                         "Vehicle", 
                         "STZ", 
                         "STZ", 
                         "STZ")
AggrSheet$condition_sex <- c("Veh_F", 
                             "Veh_M", 
                             "Veh_M", 
                             "STZ_M", 
                             "STZ_M", 
                             "STZ_M")
AggrSheet$cellID <- as.numeric(seq.int(nrow(AggrSheet))) ### adds an "index" column (1, 2, 3, etc.)
## drop useless columns
AggrSheet$molecule_h5 <- NULL
## inspect this carefully to make sure all of the info is accurate
AggrSheet
cellID <- as.numeric(gsub(".*-","", (colnames(x = E105)))) ### this variable name becomes the colname in the DF
names(cellID) <- colnames(x = E105)
## coerce named list to dataframe
metadata2add <- data.frame(cellID)
## fill out metadata DF based on annotated aggr.csv
metadata2add$rownames_for_metadata <- rownames(metadata2add) 
metadata2add <- merge(metadata2add, AggrSheet, by="cellID", all.x = TRUE, no.dups = FALSE, )
rownames(metadata2add) <- metadata2add$rownames_for_metadata
## drop columns we don't want in the metadata
metadata2add$cellID <- NULL
metadata2add$rownames_for_metadata <- NULL
head(metadata2add)
## AddMetaData
E105 <- AddMetaData(E105, metadata = metadata2add)
head(E105@meta.data)
# chenge order
E105@meta.data$condition <- factor(E105@meta.data$condition, levels = c("Vehicle", "STZ"))
levels(E105@meta.data$condition)

### Add a bit more metadata, % features (mito, ribo)
# calc % mt
E105[["percent.mt"]] <- PercentageFeatureSet(E105, pattern = "^mt-")
# calc % ribo
E105[["percent.ribo"]] <- PercentageFeatureSet(E105, pattern = "^Rp[sl]")
# cell cycle scoring
### Assign cell cycle scores
## cc gene list via Seurat tutorial: https://www.dropbox.com/s/dl/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip
cc.genes <- readLines(con = "../data/regev_lab_cell_cycle_genes.txt")
## convert Human to Mouse
cc.genes <- convertHumanGeneList(cc.genes)
convertedCCgenes = getLDS(attributes = c("hgnc_symbol"), 
                          filters = "hgnc_symbol", 
                          values = cc.genes ,
                          mart = human,
                          attributesL = c("mgi_symbol"),
                          martL = mouse,
                          uniqueRows=T)
convertedCCgenes <- unique(convertedCCgenes[, 2])

## segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:45]
g2m.genes <- cc.genes[46:100]
## assign cell cycle scores
E105 <- CellCycleScoring(E105, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(E105@meta.data)

### Save data
saveRDS(E105, file = "../data/rds/01_rawObject_noFilters_scoresAdded_06-28-2021.RDS")



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