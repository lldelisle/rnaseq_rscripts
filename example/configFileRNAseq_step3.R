### Required for all steps ###
RNAseqFunctionPath <- "~/rnaseq_rscripts/RNAseqFunctions.R"
# This file should be a tabulated file with at least one column called
# 'sample'.
samplesPlan <- "~/rnaseq_rscripts/example/samplesPlan.txt"


#### STEP 3 - PLOTS #### Required You can put here either the FPKM norm values (subsetted or not)
#### or the count norm values obtained after DESeq2
tableWithNormalizedExpression <- "~/rnaseq_rscripts/outputs/DESeq2/DESeq2AnalysisDFLPerGeno.txt"
# In case you are using a file with both raw counts and FPKM you need to choose
# which values you want to plot.  If set to T, only columns called FPKM_sample
# will be used.
useFPKM <- T
# Optional
outputFolder <- "~/rnaseq_rscripts/outputs/plots/step3/"
# By default pdf is used as output. If set to T, png will be used.
usePng <- T
# You can provide color for each value of each factor in you samples plan to
# have constistant graphs for PCA, correlation and genes.
fixedColors <- list(Line = c(Wt = "green", `del(attP-Rel5)d9lac` = "purple"), Replicate = c(`1` = "mediumturquoise",
  `2` = "orchid", `3` = "peru"))

### Common to PCA and clustering ### In DESeq2 they restrict to the 500 most
### variant genes. If you want to keep all genes, comment the line or put
### 1000000.
restrictToNMoreVariantGenes <- 500

### PCA ### Put here the number of PC you want to see (0=do not perform PCA,
### 1=Only look at first component, 2=look at the 2 first etc...)
nbOfPC <- 3
# If the nbOfPC is greater than 1, you will have a barchart of each PC and you
# may want to use different parameters to identify your samples using the
# column names of the samples plan.
PCA1D <- list(fill = "Line", alpha = "Tissue", linetype = "Replicate")
# Possible personalizations are : fill is for the color of the bar, alpha is
# for the transparency, color is for the color of the border of the bar
# linetype is for the type of border of the bar If you do not want to use one
# of the parameter, just remove it from the list.

# If the nbOfPC is greater than 2, you will have projection in 2 dimension and
# to identify your sample you may want to use the column names of the samples
# plan.
PCA2D <- list(fill = "Line", shape = "Tissue", color = "Replicate")
# Possible personalizations are : color is for the color of the symbol, alpha
# is for the transparency, shape is for the shape of the symbol. You can also
# choose 2 colors fill and color. If so, the fill will be inside and the color
# will be the border. If you do not want to use one of the parameter, just
# remove it from the list.

# Do you want to have the contribution of each gene to each PC (T=yes, F=no).
getGeneContributionToPCA <- T

### Clustering ###

# Do you want to perform a correlation matrix and clustering (T=yes, F=no)
plotMatrixAndClustering <- T

### Genes ###
# One gene per line. The first line of the gene file should correspond to a
# column in the expression file.
fileWithGenes <- "~/rnaseq_rscripts/example/genesHoxDandAround.txt"
# By default, the title of the plot is the id provided in the fileWithGenes but
# you can add a meaning full name like gene_short_name if it is provided in the
# tableWithNormalizedExpression.
# geneIDToAdd<-'gene_short_name'
# By default, the values of expression plotted are log2(1+expression) (when T,
# if F the raw expression will be plotted.)
useLogExpression <- T
# By default, each gene is plotted on an adjusted scale. If
# useSameYmaxForAllGenes is T, all genes will be plotted with the same y axis.
useSameYmaxForAllGenes <- T
# A factor which will be used as x axis.
xaxisForGenes <- "Tissue"
plotGenesPara <- list(color = "Line", shape = "Replicate")
# Possible personalizations are :
# color is for the color of the symbol, 
# alpha is for the transparency,
# shape is for the shape of the symbol.
# You can also choose 2 colors fill and color. If so, the fill will be inside and the color will be the border.
# If you do not want to use one of the parameter, just remove it from the list.

# If you only want a heatmap and not one gene per one gene. Put it to T.
doNotPlotGeneByGene <- F
# If you want to have a heatmap with the genes provided in the list. All values
# of plotGenesPara will be used to annotate the samples.
addGlobalHeatmap <- T
# By default, genes are clustered by euclidean distance and complete
# clustering. If you want to keep the original order. Put keepGeneOrder to T.
keepGeneOrder <- F
# By default, samples are not clustered.
clusterSamples <- F
