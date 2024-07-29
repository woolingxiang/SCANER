### Description
We developed a marker-free approach, Seed-Cluster based Approach for NEoplastic cells Recognition (SCANER), to identify malignant cells on the basis of genomic instability induced violent gene expression variations. By examining different cancer types, we demonstrated that the SCANER showed remarkable robustness to extensive dropout events (0.1~0.25), and capable of precisely distinguishing tumor from non-tumor cells in diverse tumor tissues sequenced by Smart-seq2 or 10X Genomics. Please visit [Tutorial-for-Using-SCANER](https://github.com/woolingxiang/SCANER/wiki/Tutorial-for-Using-SCANER) for a detailed tutorial. 
![WX20240729](https://github.com/user-attachments/assets/9b3a9623-5a58-43e7-bcae-497360ff6b4c)


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

---

### Usage
#### main function: identificationTumorCellsX
#### parameters:

* **mat.f**: An expression profile derived from scRNA-seq data.

* **cluster.tag**: A vector containing tags for each cluster, e.g., c(1,1,1,2,3,3,3,2,2,2,4,5,4,4,4,4,5,5,5,5,...). These tags should correspond to the order of cells in the expression profile. Alternatively, when the function is set to use *seurat_obj*, this parameter only accepts a single string which should be the name of a column in the meta.data of the Seurat object, such as *seurat_clusters*.

* **seurat_obj**: Alternatively, a Seurat object (< V5.0) can be used.

* **chrpq**: A data frame containing gene information on chromosomes, obtainable via SCANER::CpRMAP_GRCh38 or SCANER::CpRMAP_GRCh37.

* **GE.frac**: The cutoff for the proportion of expressed genes. Default is 0.05.

* **smooth.ceiling**: This smooths the expression profile to avoid the influence of abnormal expression on CNV inference. Default is 1.0.

* **windows**: The number of neighboring genes considered during tiling. Default is 100.

* **cluster.NUM**: The cutoff for excluding clusters with a limited number of cells. Default is 30.

* **putativeT.frac**: A cluster is considered a seed when the number of candidate cells within the cluster exceeds a given cutoff. Default is 0.6.

* **putativeT.cor**: A cell is considered a putative tumor cell when the correlation coefficient between the cell and seed exceeds a given cutoff. Recommended is 0.5.

* **inferCNVBy**: CNV inference is based on tumor cells (putativeTumor) or tumor clusters (putativeTumor2). Default is putativeTumor2.

* **verbose**: Displays detailed information about the analysis.

---

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
# you can donwload a mapDataFrame from Ensembl,
# it should include the following columns: 'ENSG','GENE','CHR','PQ','START','END'
CpRMAP = adjust.chrpq.map(mapDataFrame)
rs = identificationTumorCellsX(seurat_obj=seuratObj,chrpq=CpRMAP, putativeT.cor=0.5)
```
#### Plot CNV
```{r}
CNVPlot(rs,putativeTumor=2,genomeSizeFile='hg38',outPDF='./scanerCNV.pdf')
```

