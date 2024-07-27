CNVPlot <- function(scaner.rs = NULL, putativeTumor=2, genomeSizeFile='hg19', seed = 12345, outPDF='./scanerCNV.pdf') {
  ## scaner.rs: input the result of identificationTumorCells
  ## putativeTumor, 1 or 2
  ## outPDF: CNV plot(pdf format)

  library(ComplexHeatmap)
  library(circlize)
  #library(stringr)
  set.seed(seed)

  #### 1. CNV result
  CNV <- as.matrix(scaner.rs$cnv.adjust)
  CNV1 <- CNV[substr(rownames(CNV),1,2) %in% c(paste("0", 1:9, sep=""), 10:22, 23), ]
  if (putativeTumor == 1) {
    CNV2 <- cbind(CNV1[, scaner.rs$information$putativeTumor == "tumor"],
                  CNV1[, scaner.rs$information$putativeTumor == "non-tumor"])
  } else {
    CNV2 <- cbind(CNV1[, scaner.rs$information$putativeTumor2 == "tumor"],
                  CNV1[, scaner.rs$information$putativeTumor2 == "non-tumor"])
  }
  chrInfo <- as.numeric(substr(rownames(CNV2),1,2))

  # chrLen <- read.table(genomeSizeFile)  # Genome Size
  if(genomeSizeFile=='hg19') {
    chrLen=hg19.chrom.sizes
  } else if (genomeSizeFile=='hg38') {
    chrLen=hg38.chrom.sizes
  } else {
    cat(">>>Unsupported genomeSizeFile, still using hg19<<<")
    chrLen=hg19.chrom.sizes
  }

  colnames(chrLen) <- c("Chr", "Length")
  chrLen$Chr <- as.character(chrLen$Chr)
  chrLen1 <- chrLen[chrLen$Chr %in% c(paste("chr", 1:22, sep=""), "chrX"), ]
  chrLen1$Chr[which(chrLen1$Chr == "chrX")] <- "chr23"
  chrLen1$times <- chrLen1$Length / min(chrLen1$Length)

  #### 2. annotation
  # anno.col = data.frame(chr=substr(rownames(CNV2),1,2))
  # anno.col$chr <- as.character(anno.col$chr)
  # anno.col <- anno.col[anno.col$chr %in% c(paste("0", 1:9, sep=""), 10:22, 23), ,drop = FALSE]
  # anno.col$chr[anno.col$chr == "23"] <- "X"
  # anno.col1 <- HeatmapAnnotation(df = anno.col)

  anno.row <- scaner.rs$information
  anno.row1 <- anno.row[colnames(CNV2), ]
  anno.row1$Cluster <- as.factor(anno.row1$Cluster)
  anno.row2 <- rowAnnotation(df = anno.row1, 
                             col = list(putativeTumor2 = c('tumor'='#C1212D', 'non-tumor'='gray60')))
  # color
  #   col_fun = colorRamp2(c(min(CNV1), min(CNV1)/2, 0.02, max(CNV1)/2, max(CNV1)),
  #     c("blue", "skyblue", "white", "lightcoral", "red"))
  #   col_fun = colorRamp2(c(-0.2, -0.1, -0.02, 0, 0.02, 0.1, 0.2),
  #     c("blue", "skyblue", "white", "white", "white", "lightcoral", "red"))
  col_fun = colorRamp2(c(-0.15, -0.1, -0.02, 0, 0.02, 0.1, 0.15),
                       c("blue", "skyblue", "white", "white", "white", "lightcoral", "red"))

  ## make plot
  ht_list <- list()
  ht_opt(heatmap_border = TRUE)
  for (i in 1:23) {
    if(i == 1){
      anno_row <- anno.row2
    } else {
      anno_row <- NULL
    }
    ht <- Heatmap(as.matrix(t(CNV2[chrInfo %in% i, , drop=FALSE])), col = col_fun,
                  name = "CNV", cluster_rows = FALSE, cluster_columns = FALSE,
                  show_row_names = FALSE, show_column_names = FALSE,
                  width = unit(chrLen1$times[i] * 1.5, "mm"),
                  #top_annotation = anno.col1,
                  left_annotation = anno_row)
    ht_list[[paste("ht", i, sep="")]] <- ht
  }
  ht_list2 <- Reduce("+", ht_list)
  
  pdf(outPDF, width = 9, height = 5)
  draw(ht_list2, ht_gap = unit(0.8, "mm"))
  dev.off()
  
  return(ht_list2)
  
}
