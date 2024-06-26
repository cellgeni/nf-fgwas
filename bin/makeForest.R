#!/usr/bin/env -S Rscript --vanilla

library(qvalue)

#######################################################################
#                          HANDLE PARAMETERS                          #
#######################################################################

# fixed parameters
cell_type_file <- "tss_cell_type_exp.txt"

# parse command line arguments
file_chunks <- stringr::str_sort(commandArgs(trailingOnly = TRUE), numeric = TRUE)
file_chunks <- file_chunks[!startsWith(file_chunks, "-")]  # remove Rscript options like `--vanilla`

# derive parameters
cellt_df <- read.table(cell_type_file, sep="\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
celltype <- colnames(cellt_df)[3:(ncol(cellt_df) - 1)]


# tss=read.table(tss_file, as.is=T)
# GWAS=snakemake$wildcards$runID
# RESULT_DIR=snakemake$wildcards$out_suffix
# OUTPUT_SUFFIX=snakemake$wildcards$out_suffix
# celltype=scan(snakemake$input$cell_type_list,"",sep="\n")


#######################################################################
#                             RUN SCRIPT                              #
#######################################################################

res <- NULL
for(fc in file_chunks){
	res <- rbind(
		res,
		readBin(
			file(
				file.path(fc, "feature_level.bin"),
				"rb"
			),
			"double",
			1000
		)
	)
}

i <- 2
res <- apply(res, 1, function(xx){
	mat <- matrix(xx[4:12], 3)
	H <- try(solve(mat))

	if (is.character(H)) {
		c(0, 1)
	} else {
		c(xx[i], sqrt(H[i, i]))
	}
})
res <- t(res)

res <- cbind(
	res[,1], 
	res[,1] - 1.96*res[,2], 
	res[,1] + 1.96*res[,2], 
	pchisq(
		(res[,1] / res[,2])^2, 
		1, 
		lower=F
	)
)

#res=cbind(res, qvalue(res[,4], lambda=0)$qval)
res <- cbind(
	res, 
	p.adjust(res[,4], meth="BH")
)

rownames(res) <- celltype
colnames(res) <- c("log OR","CI (lower 95%)","CI (upper 95%)","P-value","FDR")


###############################################################################
#                             PLOTTING FUNCTION	                              #
###############################################################################

Forest <- function(x, q.cutoff=1){
	#par(mar=c(3,12,1,1))
	#xt=apply(x,2,sum);
	#nam=rownames(x);

	tfall=as.matrix(x[,2:5]) #t(apply(x,1,function(y){tf=fisher.test(cbind(xt,xt-y)); c(tf$est,tf$conf,tf$p)}));
	rownames(tfall)=x[,1]
	print(tfall)
	col = c((qvalue(tfall[,4], lambda=0)$qval<0.1)+1)[qvalue(tfall[,4], lambda=0)$qval<q.cutoff&tfall[,1]>0]
	#col = c((p.adjust(tfall[,4],meth="BH")<0.5)+1)[qvalue(tfall[,4], lambda=0)$qval<q.cutoff&tfall[,1]>0]
	#tfall=tfall[qvalue(tfall[,4], lambda=0)$qval<q.cutoff&tfall[,1]>0,]
	tfall=tfall[p.adjust(tfall[,4],meth="BH")<q.cutoff&tfall[,1]>0,]
	col = col[rev(order(tfall[,1]))]
	tfall = tfall[rev(order(tfall[,1])),]
	if(is.null(rownames(tfall))){rownames(tfall)=seq(nrow(tfall))}
	tfall=tfall[rownames(tfall)!="empty",]
	xlim=c(tfall[,2:3])
	if(xlim[1]>0){xlim[1]=0}
	xlim=range(xlim[xlim>(-Inf)&xlim<Inf&!is.na(xlim)])
	plot(tfall[,1],rev(seq(nrow(tfall))),pch=15,xlim=xlim,axes=F,xlab="Enrichment",ylab="")
	segments(tfall[,2],rev(seq(nrow(tfall))),tfall[,3],rev(seq(nrow(tfall))))
	par(xpd=NA)
	text(rep(xlim[1],nrow(tfall)),rev(seq(nrow(tfall))),rownames(tfall),pos=2, col=col)
	par(xpd=F)
	axis(1)
	abline(v=0,lty=2)
	tfall
}


#######################################################################
#                             SAVE FILES                              #
#######################################################################

# result table

print(res)

write.table(
	res, 
	col=T, 
	row=T, 
	sep="\t", 
	file="enrichment.tsv"
)

# forest plot

pdf("forest_plot.pdf", width=9, height=14)
par(mar=c(3,14,1,1))

Forest(data.frame(celltype,res),1)

dev.off()



#print(res)
#res[is.na(res[,4]),4]=1

#res=cbind(res, qvalue(res[,4], lambda=0)$qval)
#rownames(res)=celltype
#colnames(res)=c("log OR","CI (lower 95%)","CI (upper 95%)","P-value","FDR")
#write.table(res, col=T, row=T, sep="\t", file=paste(GWAS,"/forest_",GWAS,".tsv",sep=""))

