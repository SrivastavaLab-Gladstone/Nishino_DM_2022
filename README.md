## Introduction

This is the code used by Dr. Tomohiro Nishino and Mr. Angelo Pelonero in the Srivastava Lab for the submission "Single Cell Multimodal Analyses Reveal Epigenomic and Transcriptomic Basis for Birth Defects in Maternal Diabetes."

Manuscript is currently available on bioRxiv: [DOI 2022.07.25.501463](https://doi.org/10.1101/2022.07.25.501463)

## Analysis
All data was processed and analyzed using Cellranger, Seurat, ArchR and supporting packages as detailed in provided scripts. See 10x Genomics documenation for [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) and [Cellranger ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac) usage.

Analysis order:

1. Process scRNA/scATAC 10x Genomics Cellranger v5.0.0 & Cellranger-atac v2.0.0 pipelines:
   - `cellranger count` & `cellranger-atac count`
   - `cellranger aggr`
2. Analyze scRNA seq data with Seurat v4.0.2 using scripts 1-5 in `scRNA-seq/*/` folder:
    - `scRNA-seq/scRNA_Script01.R`: scRNA data read-in and processing
    - `scRNA-seq/scRNA_Script02.R`: scRNA data QC filtering and clustering
    - `scRNA-seq/scRNA_Script03.R`: Mesodermal cell subset analyses
    - `scRNA-seq/scRNA_Script04.R`: Neural-crest cell subset analyses
    - `scRNA-seq/scRNA_Script05.R`: WGCNA and related statistical analysis framework
    - `scRNA-seq/scRNA_Script06.R`: PA2 & AHF2 Alx3+ subset analyses
3. Analyze scATAC data with ArchR v1.0.1 using scripts 1-7 in `scATAC-seq/` folder
    - `scATAC-seq/scATAC_Script01.R`: scATAC data read-in, processing, and ArchR project creation
    - `scATAC-seq/scATAC_Script02.R`: Identification of differentially accessible regions between clusters + motif enrichment analsysis
    - `scATAC-seq/scATAC_Script03.R`: scRNA+scATAC data integration
    - `scATAC-seq/scATAC_Script04.R`: Neural Creast and Mesoderm subset analyses
    - `scATAC-seq/scATAC_Script05.R`: Identification of differentially accessible regions between treament conditions + motif enrichment analsysis
    - `scATAC-seq/scATAC_Script06.R`: ChromVAR analysis
    - `scATAC-seq/scATAC_Script07.R`: Identification of candidate enhancers

## Data Availability

All sequencing data *will be* available via GEO/SRA: [link-provided-when-data-released](https://www.ncbi.nlm.nih.gov/sra)

