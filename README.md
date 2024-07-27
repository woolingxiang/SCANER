### Description
We developed a marker-free approach, Seed-Cluster based Approach for NEoplastic cells Recognition (SCANER), to identify malignant cells on the basis of genomic instability induced violent gene expression variations. By examining different cancer types, we demonstrated that the SCANER showed remarkable robustness to extensive dropout events (0.1~0.25), and capable of precisely distinguishing tumor from non-tumor cells in diverse tumor tissues sequenced by Smart-seq2 or 10X Genomics.

### Installation
#### Option One:
devtools::install_github('woolingxiang/SCANER')

#### Option Two:
renv::init('.')
renv::install('woolingxiang/SCANER')
renv::snapshot(type='all')

### Usage
#### main function: identificationTumorCellsX
#### parameters:
* mat.f: An expression profile derived from scRNA-seq data. 
* cluster.tag: A vector containing tags for each cluster, e.g., c(1,1,1,2,3,3,3,2,2,2,4,5,4,4,4,4,5,5,5,5,...). The order of tags should correspond to the order of cells in expression profile. 
* seurat_obj: Instead, you can also input a Seurat object (< V5.0). 
* chrpq: A data.frame format data containing gene information on chromosomes. You can get this data by SCANER::CpRMAP_GRCh38 or SCANER::CpRMAP_GRCh37
* 
