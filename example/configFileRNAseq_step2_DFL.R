### Required for all steps ###
RNAseqFunctionPath <- "~/rnaseq_rscripts/RNAseqFunctions.R"
# This file should be a tabulated file with at least one column called
# 'sample'.
samplesPlan <- "~/rnaseq_rscripts/example/samplesPlanDFL.txt"

#### STEP 2 - DESEQ 2 ANALYSIS ####
# Required
tableWithCounts <- "~/rnaseq_rscripts/outputs/mergedTables/AllHTSeqCounts_subset.txt"
# Specify here the name of the column which contains the gene IDs (they need to
# be unique).
geneIDColCounts <- "Ens_ID"
# For the DESeq2 analysis you need to specify a factor on which you want to do
# the analysis: This needs to be a name of a column of the samplesPlan file.
factor <- "Line"

# Optional
# This can be table from cufflinks or cuffdiff or Biomart to annotate genes.
# You will need to choose a file with at least one column with the Ensembl Gene IDs.
tableWithAnnotations <- "~/rnaseq_rscripts/outputs/mergedTables/AllCufflinks_Simplified_norm.txt"
# Specify here the name of the column which contains the gene IDs (it must
# match with the content of the geneID from the table with counts).
geneIDColInAnnotations <- "gene_id"
# You can also provide a gtf:
gtfFile <- "~/rnaseq_rscripts/data/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.94_ExonsOnly_UCSC_withIsland3.gtf.gz"
# Default test is Wald but you can change to likelihood ratio test (LRT) with
# reduced formula ~1. Put F to keep Wald and put T to use LRT.
changeTest <- F
outputDESeqTable <- "~/rnaseq_rscripts/outputs/DESeq2/DESeq2AnalysisDFLPerGeno.txt"
# If you want to have another table with only significant genes abs(l2FC) > 1.5
# and corrected p-value < 0.05
outputSignificantTable <- T
