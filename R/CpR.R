#!/usr/bin/Rscript
###########################################################
#


# map = read.table("ensembl88.38_chr.map.txt",header=T,sep='\t')[,-2]
# map.uq = unique(map)
# tmp = table(map.uq$Associated.Gene.Name)
# map.fn = map.uq[which(map.uq[,2] %in% names(tmp[tmp==1])),]
# map.fn = adjust.chr.map(map.fn)
# map.chr = map.fn[order(map.fn$CHR,map.fn$START,map.fn$END),]
#
# map = read.table("ensembl88.38_chrpq.map.txt",header=T,sep='\t')[,-2]
# map.uq = unique(map)
# tmp = table(map.uq$Gene.name)
# map.fn = map.uq[which(map.uq[,6] %in% names(tmp[tmp==1])),]
# map.fn = adjust.chrpq.map(map.fn[,c(1,6,2,5,3,4)])
# map.chrpq = map.fn[order(map.fn$CHR,map.fn$START,map.fn$END),]
#
# map.chrpqb = adjust.chrpqb.map(map.chrpq)
# map.chrpqb = map.chrpqb[order(map.chrpqb$CHR,map.chrpqb$START,map.chrpqb$END),]
#
# save(file='MAP.RData',list=c('map.chr','map.chrpq'))
#
# tcga.n.expr = tcga.expr[,grep('11$',colnames(tcga.expr))]
# base.norm = apply(tcga.n.expr,1,mean)
#
# save(file='BAS.RData',list=c('base.norm'))
#
# map = read.table('mart_export-6.txt',header=T,sep='\t')
# map.uq = unique(map)
# tmp = table(map.uq$Gene.name)
# map.fn = map.uq[which(map.uq[,2] %in% names(tmp[tmp==1])),]
# map.chr = adjust.chr.map(map.fn[,-4])
# GRCh37.chr = map.chr[order(map.chr$CHR,map.chr$START,map.chr$END),]
# map.chr = adjust.chrpq.map(map.fn)
# GRCh37.chrpq = map.chr[order(map.chr$CHR,map.chr$START,map.chr$END),]
# GRCh37.chrpqb = adjust.chrpqb.map(GRCh37.chrpq)


# make sure the input file with following format:
#
# ENSG ID, Gene Symbol, Chromosome, Position Start & End
#
adjust.chr.map <- function(map){
	colnames(map) = c('ENSG','GENE','CHR','START','END')
	rownames(map) = map$GENE
	map$CHR = as.vector(map$CHR)
	map[which(map$CHR==1),]$CHR='01'
	map[which(map$CHR==2),]$CHR='02'
	map[which(map$CHR==3),]$CHR='03'
	map[which(map$CHR==4),]$CHR='04'
	map[which(map$CHR==5),]$CHR='05'
	map[which(map$CHR==6),]$CHR='06'
	map[which(map$CHR==7),]$CHR='07'
	map[which(map$CHR==8),]$CHR='08'
	map[which(map$CHR==9),]$CHR='09'
	map[which(map$CHR==10),]$CHR='10'
	map[which(map$CHR==11),]$CHR='11'
	map[which(map$CHR==12),]$CHR='12'
	map[which(map$CHR==13),]$CHR='13'
	map[which(map$CHR==14),]$CHR='14'
	map[which(map$CHR==15),]$CHR='15'
	map[which(map$CHR==16),]$CHR='16'
	map[which(map$CHR==17),]$CHR='17'
	map[which(map$CHR==18),]$CHR='18'	
	map[which(map$CHR==19),]$CHR='19'
	map[which(map$CHR==20),]$CHR='20'
	map[which(map$CHR==21),]$CHR='21'
	map[which(map$CHR==22),]$CHR='22'
	map[which(map$CHR=='X'),]$CHR='23'
	map[which(map$CHR=='Y'),]$CHR='24'
	map$CHR = as.factor(map$CHR)
	return(map)
}

adjust.chrpq.map <- function(map){
	colnames(map) = c('ENSG','GENE','CHR','PQ','START','END')
	rownames(map) = map$GENE
	ch = map$CHR
	pq = substr(map$PQ,1,1)
	map$CHR = as.vector(map$CHR)
	for(c in levels(ch)){
		if(length(which(map$CHR==c & pq=='p'))>0){
			h = ifelse(c=='X' | c=='x',23,ifelse(c=='Y' | c=='y',24,c))
			h = ifelse(as.numeric(h)<10,paste(0,h,'p',sep=''),paste(h,'p',sep=''))
			map[which(map$CHR==c & pq=='p'),]$CHR=h
		}
		if(length(which(map$CHR==c & pq=='q'))>0){
			h = ifelse(c=='X' | c=='x',23,ifelse(c=='Y' | c=='y',24,c))
			h = ifelse(as.numeric(h)<10,paste(0,h,'q',sep=''),paste(h,'q',sep=''))
			map[which(map$CHR==c & pq=='q'),]$CHR=h
		}
	}
	map$CHR = as.factor(map$CHR)
	return(map)
}

# input file derived from output of adjust.chrpq.map
adjust.chrpqb.map <- function(map){
  colnames(map) = c('ENSG','GENE','CHR','PQ','START','END')
  rownames(map) = map$GENE
  map$CHR <- apply(map[,c('CHR','PQ')],1,function(x){paste(gsub('(p|q)','',x[1]),x[2],sep='')})
  map$CHR = as.factor(map$CHR)
  return(map)
}

# main function, cnv preprocessing
cnv.preproc <- function(expr,map,model,...){
	cm.genes = intersect(rownames(expr),rownames(map))
	cm.map = map[cm.genes,]
	cm.map$CHR <- as.vector(cm.map$CHR)
	cm.map = cm.map[order(cm.map$CHR,cm.map$START,cm.map$END),]
	cm.expr = expr[rownames(cm.map),]
	barcode = paste(cm.map$CHR,'_',sep='')
	rs.md = apply(cm.expr,2,function(s,b){
		unlist(tapply(s,b,function(c){
			model(c,...)$v
		}))
	},b=barcode)
	chrs = as.factor(gsub('_.*$','',rownames(rs.md)))		
	return(list(cnv=rs.md,chr=chrs,map=cm.map))
}

# main function, cnv smoothing	
cnv.smooth <- function(cp,model,min=1,...){
  barcode = paste(cp$chr,'_',sep='')
	rs.sm = apply(cp$cnv,2,function(x,b){
		unlist(tapply(x,b,function(y,z){
			if(length(y)<=min){
				rep(z,times=length(y))
			}else{
				model(y,z,...)$v
			}
		},z=mean(x)))
	},b=barcode)
	return(list(cnv=rs.sm,chr=cp$chr))
}
	
# main function, cnv calling
cnv.calling <- function(cp,model,...){
  barcode <- paste(cp$chr,'_',sep='')
	rs = apply(cp$cnv,2,function(x,b,me){
	  unlist(tapply(x,b,function(y,z){
	    model(y,list(sampleMean=z,all.sampleMean=me),...)$v
	  },z=mean(x)))
	},b=barcode,me=mean(cp$cnv))
	return(list(cnv=rs,chr=cp$chr))
}
	
# core function, processing data
# must have x and m
cp.default.model <-  function(x,m=NULL,windows=100,win.frac=NULL,mtd=mean){
  if(!is.null(win.frac))
    windows=round(win.frac*length(x))
	rs <- c()
	if(length(x)<=windows){
		rs=mtd(x)
	}else{
		rs=sapply(1:(length(x)-windows+1),function(i){
		  mtd(x[i:(i+windows-1)])
		})	
	}
	return(list(v=rs))
}

# must have x and m
cp.wintile.model <-  function(x,m=NULL,windows=10,win.frac=NULL){
  if(!is.null(win.frac))
    windows=round(win.frac*length(x))
	rs=sapply(1:ceiling(length(x)/windows),function(i){
		mean(head(x[((i-1)*windows+1):length(x)],windows))
	})	
	return(list(v=rs))
}

# must have x and m
cs.whittaker.model <- function(x,m=NULL,...){
	library(pracma)
	# rs = head(tail(whittaker(c(rep(m,times=length(x)),x,rep(m,times=length(x))),...),2*length(x)),length(x))
	rs = whittaker(x,...)
	return(list(v=rs))
}

# must have x and m
cc.default.model <- function(x,m=NULL,placehold=0.01,...){
  rs = log2(x/m$sampleMean+placehold)
  return(list(v=rs))
}

