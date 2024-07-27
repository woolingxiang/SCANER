### Description
We developed a marker-free approach, Seed-Cluster based Approach for NEoplastic cells Recognition (SCANER), to identify malignant cells on the basis of genomic instability induced violent gene expression variations. By examining different cancer types, we demonstrated that the SCANER showed remarkable robustness to extensive dropout events (0.1~0.25), and capable of precisely distinguishing tumor from non-tumor cells in diverse tumor tissues sequenced by Smart-seq2 or 10X Genomics.

### Installation
#### Option One:
```{r}
devtools::install_github('woolingxiang/SCANER')
```
#### Option Two:
```{r}
renv::init('.')
renv::install('woolingxiang/SCANER')
renv::snapshot(type='all')
```

### Usage
#### main function: identificationTumorCellsX
#### parameters:
* mat.f: An expression profile derived from scRNA-seq data. 
* cluster.tag: A vector containing tags for each cluster, e.g., c(1,1,1,2,3,3,3,2,2,2,4,5,4,4,4,4,5,5,5,5,...). The tags should correspond to the order of cells in the expression profile. 
* seurat_obj: Instead, you can also input a Seurat object (< V5.0). 
* chrpq: A data.frame format data containing gene information on chromosomes. You can get this data by SCANER::CpRMAP_GRCh38 or SCANER::CpRMAP_GRCh37. 
* GE.frac: Cutoff for the proportion of expressed genes. Default: 0.5. 
* smooth.ceiling: smoothing expression profile, avoiding the influence of abnormal expression on CNV inferring. Default: 1.0. 
* windows: The number of neighboring genes when tiling. Default: 100.
* cluster.NUM: Cutoff for excluding clusters with a limited number of cells. Default: 30. 
* putativeT.frac: Consider a cluster as a seed when the number of candidate cells within the cluster exceeds a given cutoff. Default: 0.6. 
* putativeT.cor: Consider a cell as a putative tumor cell when the correlation coefficient between the cell and seed exceeds a given cutoff. Recommend: 0.5. 
* inferCNVBy: Infer CNV based on tumor cells (putativeTumor) or tumor clusters (putativeTumor2). Default: putativeTumor2. 
* verbose: Show the detailed information of analysis. 

### Example
#### Quick start
```{r}
# suppose you have a Seurat object
rs = identificationTumorCellsX(seurat_obj=seuratObj,chrpq=SCANER::CpRMAP_GRCh38, putativeT.cor=0.5)

# suppose you have an expression profile and corresponding cluster information
rs = identificationTumorCellsX(mat.f=expressionMatrix,cluster.tag=clusterVector,chrpq=SCANER::CpRMAP_GRCh38, putativeT.cor=0.5)

# see the results
head(rs$information)
```
#### Create your own chromosome dictionary
```{r}
# you can donwload a mapDataFrame from Ensembl, it should include the following columns: 'ENSG','GENE','CHR','PQ','START','END'
CpRMAP = adjust.chrpq.map(mapDataFrame)
rs = identificationTumorCellsX(seurat_obj=seuratObj,chrpq=CpRMAP, putativeT.cor=0.5)
```
### Plot CNV
```{r}
CNVPlot(rs,putativeTumor=2,genomeSizeFile='hg38',outPDF='./scanerCNV.pdf')
```

