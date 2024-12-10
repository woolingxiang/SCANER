#!/usr/bin/Rscript
#############################################

identificationTumorCellsX <- function(mat.f=NULL,cluster.tag=NULL,seurat_obj=NULL,
                                      chrpq=NULL,GE.frac=0.05,smooth.ceiling=1,
                                      windows=100,cluster.NUM=30,cores=4,
                                      putativeT.frac=0.6,putativeT.cor=NULL,
                                      inferCNVBy='putativeTumor2',verbose=T){
  
  # initiating ...
  # out.dir=''
  # GE.frac=0.05
  # smooth.ceiling=1
  # g.version='grch37'
  # windows=100
  # cluster.NUM=30
  # putativeT.frac=0.6
  # putativeT.cor=0.5
  
  
  # handle with Seurat Object, require Seurat version < 5
  if(!is.null(seurat_obj)){
    cat('[INFO]',date(),'Detect as a seurat object\n')
    mat.f = as.matrix(seurat_obj[['RNA']]@data)
    if(length(cluster.tag)==1 & length(which(colnames(seurat_obj@meta.data)==cluster.tag))==1){
      cluster.tag = seurat_obj@meta.data[[cluster.tag]]
      names(cluster.tag) = rownames(seurat_obj@meta.data)
    }
  }
  
  # handle with matrix file
  if(ncol(mat.f) != length(cluster.tag)){
    stop(paste0("Sorry, difference in sample size between profile (",
                ncol(mat.f),") and cluster information (",length(cluster.tag),")\n"))
  }
  
  # create object for subsequent analysis
  cat('[INFO]',date(),'create a standard object for analysis\n')
  info.f = data.frame(Cluster=cluster.tag,row.names=names(cluster.tag))
  info.f = as.data.frame(as.matrix(info.f)[colnames(mat.f),])
  colnames(info.f) = 'Cluster'
  
  # two criteria should be kept in mind
  cat('[INFO]',date(),'criteria evaluating ...\n')
  stat2 <- apply(mat.f,1,function(x){length(which(x>0))/length(x)})		# prop. of cells expressing the indicated gene
  tst.expr.f <- mat.f[which(stat2>GE.frac),]					# remove the genes rarely expressed in cells
  c1.detectedFrac = apply(tst.expr.f,2,function(x){length(which(x>0))/length(x)})	# prop. of genes detected in the indicated cell
  dat.f <- apply(tst.expr.f,2,function(x,m){x[which(x>m)]=m;x},m=smooth.ceiling)	# smooth profile
  cp  = cnv.preproc(dat.f,chrpq,cp.default.model,windows=windows,mtd=mean,cores=cores)	# calculate raw copy number
  c2.cnvMad = apply(cp$cnv,2,mad)							# instability of copy number
  info.f$putativeSeed = factor(ifelse(c1.detectedFrac>mean(c1.detectedFrac) & c2.cnvMad>mean(c2.cnvMad),'seed','non-seed'),levels=c('non-seed','seed'))
  
  # select the candidated cluster as tumor seeds
  cat('[INFO]',date(),'choose tumor seeds\n')
  stat = table(info.f$Cluster,info.f$putativeSeed)
  stat.frac = stat[,2]/rowSums(stat)			# prop. of cells marked as tumor-seed in each cluster
  tmp = stat.frac[which(rowSums(stat)>cluster.NUM)]	# remove the clusters consisting of few cells
  seed = unique(c(names(tmp)[which(tmp>putativeT.frac)],names(tmp)[which.max(tmp)]))	# at least, one cluster should be marked as tumor-seed [CAUTION!]
  seed.cnv = apply(cp$cnv[,which(info.f$Cluster%in%seed)],1,mean)				# create a cnv seed based on cluster-seed
  pool = apply(cp$cnv,2,function(x,r){
    cor.test(x,r)$estimate
  },r=seed.cnv)
  
  # generalized seed
  cat('[INFO]',date(),'generalize tumor seeds ...\n')
  if(is.null(putativeT.cor)){
    dis.seed = pool[which(info.f$Cluster %in% seed)]				# get all correlation coeff. in cluster-seed
    exclude.frac = ifelse(putativeT.frac>0.5,1-putativeT.frac,0.5)			# in common, putativeT.frac should be greater than 0.5
    putativeT.cor = qnorm(exclude.frac,mean(dis.seed),sd(dis.seed))			# get the cutoff of correlation coeff. for excluding non-tumor cells
    putativeT.cor = ifelse(putativeT.cor>1 | putativeT.cor< -1, 0.5, putativeT.cor)	# in case of exceeding valid range, set a default value
  }
  info.f$putativeTumor = factor(ifelse(pool>putativeT.cor,'tumor','non-tumor'),levels=c('non-tumor','tumor'))
  
  # generalized cluster
  mark = apply(table(info.f$putativeTumor,info.f$Cluster),2,which.max)-1
  info.f$putativeTumor2 = 'non-tumor'
  info.f[which(info.f$Cluster %in% names(mark)[which(mark>0)]),]$putativeTumor2 = 'tumor'
  info.f$putativeTumor2 = factor(info.f$putativeTumor2,levels=c('non-tumor','tumor'))
  
  pool2 = pool
  pool2[which(pool[info.f$putativeTumor2=='tumor']<putativeT.cor)] = putativeT.cor
  pool2[which(pool[info.f$putativeTumor2=='non-tumor']>putativeT.cor)] = putativeT.cor
  
  baseline = apply(cp$cnv[,which(info.f[[inferCNVBy]]=='non-tumor')],1,mean)
  cnv = apply(cp$cnv,2,function(x){x-mean(x)})
  cnv.adjust = apply(cnv,2,function(x){x-baseline+mean(baseline)})
  
  
  if(verbose){
    cat('[LOG PARAMETERS] ==================================================\n')
    cat('[PARAM SMOOTH CEILING]:',smooth.ceiling,'\n')
    cat('[PARAM WINDOWS]:',windows,'\n')
    cat('[PARAM FRACTION OF EXPRESSED GENES CUTOFF]:',GE.frac,'\n')
    cat('[PARAM PUTATIVE TUMOR CORRELATION CUTOFF]:',putativeT.cor,'\n')
    cat('[PARAM PUTATIVE TUMOR FRACTION CUTOFF]:',putativeT.frac,'\n')
    cat('[PARAM CELL NUMBER OF CLUSTER CUTOFF]:',cluster.NUM,'\n')
    cat('[LOG STATISTICS] ==================================================\n')
    cat('[STAT NUMBER OF GENES AND CELLS]:',dim(mat.f),'\n')
    cat('[STAT CLUSTERS]:',levels(as.factor(info.f$Cluster)),'\n')
    cat('[STAT NUMBER OF CELLS PER CLUSTER]:',table(as.factor(info.f$Cluster)),'\n')
    cat('[STAT FRACTION OF EXPRESSED GENES PER CLUSTER]:',
        tapply(c1.detectedFrac,info.f$Cluster,mean),'\n')
    cat('[STAT MAD OF CN PER CLUSTER]:',
        tapply(c2.cnvMad,info.f$Cluster,mean),'\n')
    cat('[STAT TUMOR CLUSTER MASK]:',mark,'\n')
    cat('[STAT TUMOR CLUSTER SEED FRAC]:',stat.frac,'\n')
    cat('[STAT TUMOR CLUSTER SEED]:',seed,'\n')
    cat('[STAT TUMOR PURITY (putativeTumor/putativeTumor2)]:',
        table(as.factor(info.f$putativeTumor))[2]/sum(table(as.factor(info.f$putativeTumor))),
        table(as.factor(info.f$putativeTumor2))[2]/sum(summary(as.factor(info.f$putativeTumor2))),'\n')
  }
  
  
  return(list(
    information = info.f,
    seed.frac = stat.frac,
    cnv.raw = cp,
    cnv.basline = baseline,
    cnv.adjust = cnv.adjust,
    pool = pool,
    pool2 = pool2
  ))
  
  
}






