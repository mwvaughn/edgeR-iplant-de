#!/usr/bin/env Rscript

# edgeR using GLM
#
# Provides the ability to do replicated (or even unreplicated)
# pairwise comparisons of tag libraries
# using General Linear Models

# Example invocation
#./run_edgeR_glm.R --input join_tool_output_fake.txt --namecol 1 --factors "1,1,2,2" \
# --minsumcount 5 --mtc "BH" --fdr 0.05 --pair "1,2"

# Example infile
#Gene	AA	AA	BB	BB	CC	CC
#AC148152.3_FG005	12	10	219	222	219	222
#AC148152.3_FG008	40	37	15	18	15	36

# Auto-set up dependencies
# This is for API deployment
#
#.libPaths(getwd())
#install.packages("getopt", repos="http://lib.stat.cmu.edu/R/CRAN", lib=getwd())
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR", destdir=getwd(), lib=getwd())
#.libPaths(getwd())

# Use pre-installed libraries. This is a UA local cluster deployment
.libPaths("/usr/local3/bin/edgeR-1.00/R")
#install.packages("getopt", repos="http://lib.stat.cmu.edu/R/CRAN")
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")

library("getopt")
library("edgeR")

example.help <-function(){
	cat("Help is on its way!\n")
	quit(status=0)
	
}

args<-commandArgs(TRUE)

options<-matrix(c(	'input','i',1,"character",
					'namecol','c',1,"integer",
					'factors','x',1,"character",
					'minsumcount','k',1,"integer",
					'dispersion','d',1,"double",
					'pair','a',1,"character",					
					'mtc','m',1,"character",
					'fdr','p',1,"double"), ncol=4,byrow=TRUE)
					
ret.opts<-getopt(options,args)
if(is.null(ret.opts) || !is.null(ret.opts$help)){
	example.help()
}

tableFile <- ret.opts$input
tableRowNames <- ret.opts$namecol

# Provide later if needed
tableDelim <- "\t"
tableQuote <- "\""
tableDec <- "."
tableComment <- "#"
tableHeader <- FALSE

# Comma separated list of factors for each column
groupFactor <- ret.opts$factors
rowsumFilter <- ret.opts$minsumcount

# Dispersion
# -1 means auto-calculate it
# Otherwise, this is a float
if (is.null(ret.opts$dispersion)) {
	declaredDispersion <- -1
} else {
	declaredDispersion <- ret.opts$dispersion
}

# Pairing behavior
# Default behavior is to compare first two factors
if (is.null(ret.opts$pair)) {
	pairFactor <- NULL
} else {
	pairFactor <- unlist(strsplit(ret.opts$pair, ","))
}

# multiple testing method
mtcMethod <- ret.opts$mtc

# threshold for FDR filtering
minFDR <- ret.opts$fdr

x <- read.delim(tableFile, row.names=tableRowNames, sep=tableDelim, quote=tableQuote, dec=tableDec, comment=tableComment)

# User passes in factor string "1,2" for example
# How do I turn that into a matrix or vector automatically
gf <- unlist(strsplit(groupFactor, ","))
group <- factor(gf)

# Set up factored design matrix
design <- model.matrix(~group)

# Set of the DGEList 
D <- DGEList(counts=x,group=group)

# Ref column has to be numeric. Same as for step 1
# I don't know how to make use of this in downstream analysis
# Authors of edgeR seem to think its usually not needed so we
# can come back to this

# Normalize
D <- calcNormFactors(D)

# estimateCommonDisp(object, tol=1e-06, rowsum.filter=5)
# rowsum.filter is user-specified. It's the 
# value for the filtering out of low abundance tags in 
# the estimation of the common dispersion

# method (low-level function) used to estimated the trended dispersions. Possible values are "bin.spline" (default), "bin.loess" (which both result in a call to dispBinTrend), "power" (call to dispCoxReidPowerTrend), or "spline" (call to dispCoxReidSplineTrend).
D <- estimateGLMTrendedDisp(D,design, method="bin.spline")
D <- estimateGLMTagwiseDisp(D,design)

# GLM Fit
# Allow user to declare their own dispersion value or 
# use the one calculated by edgeR

if (declaredDispersion >= 0) {
	cat("GLM using user-declared dispersion\n")
	#et <- exactTest(D, dispersion=declaredDispersion, pair=pairFactor)
	fit <- glmFit(D,design,dispersion=declaredDispersion)
	lrt <- glmLRT(D,fit,coef=2)
} else {
	cat("GLM using edgeR-calculated dispersion\n")
	#et <- exactTest(D, pair=pairFactor)
	fit <- glmFit(D,design,dispersion=declaredDispersion)
	lrt <- glmLRT(D,fit,coef=2)
}

# pass in adjust string: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none"
tt <- topTags(lrt, n=dim(D$counts)[1], adjust.method=mtcMethod )
write.table(tt, file=paste("edgeR-glm-all", "txt", sep=".") , sep="\t", quote=FALSE)

# Filter toptable by FDR
tt2 <- tt$table[tt$table[4] <= minFDR,1:5]
write.table(tt2, file=paste("edgeR-glm-significant", "txt", sep=".") , sep="\t", , quote=FALSE)
