### Required for all steps ###
RNAseqFunctionPath <- "~/rnaseq_rscripts/RNAseqFunctions.R"

#### STEP 4 - PLOTS RESULTS OF DESEQ2 ####
# Required You can put here the result of step2 or any table with at least 2
# columns (padj and log2FoldChange) for Volcano and 3 columns (padj,
# log2FoldChange and baseMean) for MAP.
tableWithResultsOfDifferentialExpression <- "~/rnaseq_rscripts/outputs/DESeq2/DESeq2AnalysisDFLPerGeno.txt"

outputFolder <- "~/rnaseq_rscripts/outputs/plots/step4"
# By default pdf is used as output. If set to T, png will be used.
usePng <- T

### Common for Volcano and MAP ###
# significant thresholds:
# Put the maximum value of ajusted p-value to be significant (default 0.05)
maxPAdj <- 0.05
# Put the minimum log2 fold change to be significant (default 1.5)
minLFC <- 1.5
# colors
colOfNonSignificant <- "grey"
colOfSignificant <- "blue"
# Do you want to click on the plot to be able to identify some genes. (T=yes,
# F=no)
click <- T
# If you want to click, provide here the name of the column which can be used
# to label the gene.
geneID <- "gene_name"
# You can also provide a list of genes you want to identify on the plots. One
# gene per line. The first line of the gene file should correspond to a column
# in the tableWithResultsOfDifferentialExpression file.
fileWithGenes <- "~/rnaseq_rscripts/example/genesHoxDandAround.txt"
colOfCircle <- "red"

### Volcano ###
# If you want to zoom on bigger p-values because you have some genes with very
# low p-values, you can put here a value to restrict the plot to p-values
# higher than 10^-maxYVolcano. Put NA or comment line if you do not need to zoom
# the initial plot.
maxYVolcano <- NA

### MAP ###
# If you want to zoom on smaller log2 fold changes because you have some genes
# with very high log2 fold changes, you can put here minimum and maxiumum log2
# fold change values to restrict the plot to these values. Put NA or comment
# line if you do not need to zoom the initial plot, put c(minValue,maxValue)
# if you want to restrict.
ylimMAP <- NA
