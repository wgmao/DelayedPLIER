# DelayedPLIER

## Example

<b>#We first create a toy example using an available large-scale scRNA-seq</b>
library(TENxPBMCData)
sce <- TENxPBMCData("pbmc4k")
mat <- counts(sce)
rownames(mat) <- rowData(sce)$ENSEMBL_ID 
colnames(mat) <- colData(sce)$Sequence
writeHDF5Array(mat, filepath = "counts.hdf5", name = "count")
saveRDS(list(row.names = rowData(sce)$ENSEMBL_ID , col.names = colData(sce)$Sequence), file = "dimnames.RDS")
<b>#In the latest version of HDF5Array, dimension names can be written in to the hdf5 file directly </b>
<b>#writeHDF5Array(mat, filepath = "counts.hdf5", name = "count", with.dimnames=T) </b>

<b>Load all functions in funcs.R</b>
source("funcs.R") #glmnet and HDF5Array will be loaded 
setRealizationBackend("HDF5Array") #supportedRealizationBackends(), getRealizationBackend()

<b>read in the file we just created. The name can be verified via h5ls("counts.hdf5") </b>
sce <- DelayedArray(seed = HDF5ArraySeed(filepath = "counts.hdf5", name = "count"))


#Read the dimension names
dimnamaes <- readRDS("dimnames.RDS")
rownames(sce) <- dimnamaes$row.names
colnames(sce) <- dimnamaes$col.names


#remove all NA values and genes that are not differentially expressed
sce[is.na(sce)] <- 0
data <- sce[which(DelayedMatrixStats::rowSds(sce) >0),]

#read in prior information matrix
priorMat <- readRDS("canonicalPathways_ENSG.RDS")


#run the DelayedPLIER
#PLIER.res$residual, PLIER.res$Z and PLIER.res$B are stored as three on-disk files instead of variables in the memory. PLIER() function will realize and write
#these three to hard disk based on the output_path. Other values are still in the memory can be directly accessed via PLIER.res$'name of the values'
ptm <- proc.time()
PLIER.res <- PLIER(data, priorMat, output_path = "output/")
print(proc.time()-ptm)



## Tutorial

### Basic Functions & Tutorial
1. https://www.bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/UsingHDF5Array.html
2. https://bioconductor.github.io/BiocWorkshops/effectively-using-the-delayedarray-framework-to-support-the-analysis-of-large-datasets.html#
3. http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/DelayedArrayWorkshop__Effectively_using_the_DelayedArray_framework_for_users/


### Available Delayed Operators
1. DelayedArray-utils
https://www.rdocumentation.org/packages/HDF5Array/versions/1.0.2/topics/DelayedArray-utils
2. DelayedMatrixStats
https://www.bioconductor.org/packages/release/bioc/vignettes/DelayedMatrixStats/inst/doc/DelayedMatrixStatsOverview.html
3. BiocSingular
https://master.bioconductor.org/packages/devel/bioc/manuals/BiocSingular/man/BiocSingular.pdf
