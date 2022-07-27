### Script06: PA2 and AHF2 Alx3+ subset analyses
### Load libraries and data ####
library(Seurat)
set.seed(7)

### PA2 subset analysis ####
## Subset and cluster
NC <- readRDS(file = "/path/to/NC_Subset.RDS")
DimPlot(NC)
PA2subset <- subset(NC, idents = "PA2")
# rm(NC)
PA2subset <- SCTransform(PA2subset, assay = "RNA", new.assay.name = "SCT", variable.features.n = 4000, variable.features.rv.th = 1.3, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), return.only.var.genes = FALSE)
PA2subset <- RunPCA(PA2subset, verbose = FALSE, npcs = 100)
PA2subset <- harmony::RunHarmony(PA2subset, group.by.vars = "gem.group", assay.use="SCT")
dims <- 1:22
PA2subset <- RunUMAP(PA2subset, reduction = "harmony", dims = dims)
PA2subset <- FindNeighbors(PA2subset, reduction = "harmony", dims = dims)
PA2subset <- FindClusters(PA2subset, resolution = 0.16, random.seed = 7)
DimPlot(PA2subset, label = TRUE) + DimPlot(PA2subset, label = TRUE, group.by = "condition")

## DEG in cluster0 between conditions (STZ vs VEH)
PA2subset@meta.data$celltype_PA2subset <- PA2subset@active.ident
celltype_PA2subset_Con_vec <- paste0(PA2subset@active.ident,"_", substr(PA2subset@meta.data$condition, 1,3))
PA2subset@meta.data$celltype_PA2subset_Con <- celltype_PA2subset_Con_vec
PA2subset <- SetIdent(PA2subset, value = "celltype_PA2subset_Con")
table(PA2subset@active.ident)
Cluster0_SvsV <- FindMarkers(object = PA2subset, ident.1 = "0_STZ", ident.2 = "0_Veh", max.cells.per.ident = 616,logfc.threshold = 0.125)
cluster0_stzVveh <- cluster0_stzVveh[cluster0_stzVveh$p_val_adj < 0.05 & abs(cluster0_stzVveh$avg_log2FC) > 0.25,]

## end NC subset analysis ##




### Alx3 positive AHF2 subset analysis ####
## Subset and cluster
Meso <- readRDS(file = "/path/to/Meso_Subset.RDS")
Alx3PosAHF2 <- subset(x = Meso, subset = Alx3 > 0, idents = "*AHF2")
# rm(Meso)
Alx3PosAHF2 <- SCTransform(Alx3PosAHF2, assay = "RNA", new.assay.name = "SCT", variable.features.n = 4000, variable.features.rv.th = 1.3, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), return.only.var.genes = FALSE)
Alx3PosAHF2 <- RunPCA(Alx3PosAHF2, verbose = FALSE, npcs = 100)
Alx3PosAHF2 <- harmony::RunHarmony(Alx3PosAHF2, group.by.vars = "gem.group", assay.use="SCT")
dims <- 1:30
Alx3PosAHF2 <- RunUMAP(Alx3PosAHF2, reduction = "harmony", dims = dims)
Alx3PosAHF2 <- FindNeighbors(Alx3PosAHF2, reduction = "harmony", dims = dims)
Alx3PosAHF2 <- FindClusters(Alx3PosAHF2, resolution = 0.16, random.seed = 7)
DimPlot(Alx3PosAHF2, label = TRUE) + DimPlot(Alx3PosAHF2, label = TRUE, group.by = "condition")

## Check cell populations
AHF2sub_poptable <- table(Alx3PosAHF2@active.ident, Alx3PosAHF2@meta.data$gem.group)
DimPlot(Alx3PosAHF2, label = FALSE, cols = c('0' = "#66C2A5", '1' = "#5E4FA2",'2' = "#9E0142"))
DimPlot(object = Alx3PosAHF2, group.by = 'condition', cols = c('STZ' = '#D6604D', 'Vehicle' = '#4393C3'), label = FALSE, pt.size = 0.8,shuffle = TRUE)
FeaturePlot(Alx3PosAHF2, features = c("Crabp1", "Crabp2"), keep.scale = "all", order = T, min.cutoff = 0, blend = T, blend.threshold = 1, cols = c("#117733",  "#CC6677"), split.by = "condition")
## Subset on Hoxb1 expression if desired
# Alx3PosAHF2_Hoxb1Pos <- subset(x = Alx3PosAHF2, subset = Hoxb1 > 0)
# Alx3PosAHF2_Hoxb1Nega <- subset(x = Alx3PosAHF2, subset = Hoxb1 == 0)

## Run DEG between STZ vs VEH
Idents(Alx3PosAHF2) <- "condition"
Alx3PosAHF2_DEGs_SvsV <- FindMarkers(Alx3PosAHF2 , ident.1 = "STZ", ident.2 = "Vehicle", max.cells.per.ident = min(table(Alx3PosAHF2@meta.data$condition)))

## end AHF2 subset analysis ##