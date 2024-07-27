#!/usr/bin/Rscript
#############################################

identificationTumorCells <- function(mat.f,cluster.tag,chrpq=NULL,GE.frac=0.05,smooth.ceiling=1,
 					 windows=100,cluster.NUM=30,
					 putativeT.frac=0.6,putativeT.cor=NULL,
					 inferCNVBy='putativeTumor2',verbose=T){

# initiating ...
#out.dir=''
#GE.frac=0.05
#smooth.ceiling=1
#g.version='grch37'
#windows=100
#cluster.NUM=30
#putativeT.frac=0.6
#putativeT.cor=0.5

#Args <- commandArgs(T)
#if(length(Args)==0){
#  return(NULL)
#}else{
#  for(cmd in Args) eval(parse(text=cmd))
#}

if(ncol(mat.f) != length(cluster.tag)){
	cat("[ERROR]: Sorry, difference in sample size between mat.f (",
		ncol(mat.f),") and cluster.tag (",length(cluster.tag),")\n")
	return(NULL)
}

info.f = data.frame(Cluster=cluster.tag,row.names=names(cluster.tag))
info.f = as.data.frame(as.matrix(info.f)[colnames(mat.f),])
colnames(info.f) = 'Cluster'

# source('CpR.R')
# load('CpRMAP.RData')

#
# chrpq=c()
# if(g.version=='grch37'){
# 	chrpq=GRCh37.chrpq
# }else if(g.version=='grch38'){
# 	chrpq=GRCh38.chrpq
# }else{
# 	cat("[ERROR]: Sorry, cannot recognize version (",g.version,")\n")
# 	return(NULL)
# }

stat2 <- apply(mat.f,1,function(x){length(which(x>0))/length(x)})
tst.expr.f <- mat.f[which(stat2>GE.frac),]
dat.f <- apply(tst.expr.f,2,function(x,m){x[which(x>m)]=m;x},m=smooth.ceiling)
tmp = dat.f
cp  = cnv.preproc(tmp,chrpq,cp.default.model,windows=windows,mtd=mean)
c1.detectedFrac = apply(tst.expr.f,2,function(x){length(which(x>0))/length(x)})
c2.cnvMad = apply(cp$cnv,2,mad)
info.f$putativeSeed = ifelse(c1.detectedFrac>mean(c1.detectedFrac) & c2.cnvMad>mean(c2.cnvMad),'seed','non-seed')

stat = table(info.f$Cluster,info.f$putativeSeed)
stat.frac = stat[,2]/rowSums(stat)
tmp = stat.frac[which(rowSums(stat)>cluster.NUM)]
seed = unique(c(names(tmp)[which(tmp>putativeT.frac)],names(tmp)[which.max(tmp)]))
seed.cnv = apply(cp$cnv[,which(info.f$Cluster%in%seed)],1,mean)
pool = apply(cp$cnv,2,function(x,r){
  cor.test(x,r)$estimate
},r=seed.cnv)
info.f$putativeTumor = ifelse(pool>putativeT.cor,'tumor','non-tumor')

table(info.f$putativeTumor,info.f$Cluster)
mark = apply(table(info.f$putativeTumor,info.f$Cluster),2,which.max)-1
info.f$putativeTumor2 = 'non-tumor'
info.f[which(info.f$Cluster %in% names(mark)[which(mark>0)]),]$putativeTumor2 = 'tumor'

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






