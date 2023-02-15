### Required for all steps ###
RNAseqFunctionPath <- "~/rnaseq_rscripts/RNAseqFunctions.R"
# This file should be a tabulated file with at least one column called
# 'sample'. Optionnaly, the paths to the counts tables and FPKM tables can be
# provided under the column called: htseq_count_file and cufflinks_file.
samplesPlan <- "~/rnaseq_rscripts/example/samplesPlan.txt"

#### STEP 1 - MERGE TABLES ####
# If the merged tables are not already generated:
outputFolderForStep1 <- "~/rnaseq_rscripts/outputs/mergedTables/"
# Needed for DESeq2: Do you want to merge counts? T=yes F or commented=no
mergeCounts <- T
# Optional: subset the count table Do you want to remove some genes from the
# count table
subsetCounts <- T
# If the table with counts have already been generated and you just want to
# remove some genes.
# initialTableWithCount<-'~/rnaseq_rscripts/example/mergedTables/AllHTSeqCounts.txt'
# If you provide the initialTableWithCount you need to provide the name of the
# column with the ensembl id.
# geneIDColInInitialTable<-'Ens_ID'
# List of genes id to remove (one per line with no header).
genesToRmFromCounts <- "~/rnaseq_rscripts/example/MTgenes.txt"
# Optional:
mergeFPKM <- T
# By default cufflinks split the transcripts which do not overlap in different
# locus and so different lines, put T if you want to sum the FPKM for non
# overlapping transcripts (put F if not).
oneLinePerEnsemblID <- T
# Optional: subset the FPKM table Do you want to remove some genes from the
# FPKM table
subsetFPKM <- T
chrToRemove <- c("chrM")
# Anouk method: Genes that have the less variable rank should have the same
# expression.
normFPKMWithAnoukMethod <- T
# If the table with FPKM have already been generated and you just want to
# normalize it.
# initialTableWithFPKM<-'~/rnaseq_rscripts/example/mergedTables/AllCufflinks_Simplified.txt'
# Usually, it is recommanded to remove mitochondrial genes before doing the
# normalization. In some cases, it can also be useful to remove the sex
# chromosomes (put c('chrX','chrY','chrM')).  If you do not want to remove any
# gene put NA or comment the line.
chrToRemoveBeforeNormWithAnoukMethod <- c("chrM")
# Default is 1000, you can change here.
nbOfGenesWithAnoukMethod <- 1000
# If you want to keep the genes used in the normalization from Anouk, they will
# be written in a file.
keepGenesUsedForNorm <- F
