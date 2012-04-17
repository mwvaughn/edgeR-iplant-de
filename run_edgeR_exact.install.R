#!/usr/bin/Rscript

# Run as root

.libPaths("/usr/local3/bin/edgeR-1.00/R")
install.packages("getopt", repos="http://lib.stat.cmu.edu/R/CRAN")
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

