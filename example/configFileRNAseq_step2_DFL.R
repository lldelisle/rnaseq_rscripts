### Required for all steps ###
RNAseqFunctionPath<-"~/rnaseq_rscripts/RNAseqFunctions.R"
samplesPlan<-"~/rnaseq_rscripts/example/samplesPlanDFL.txt" #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 2 - DESEQ 2 ANALYSIS ###
#Required
tableWithCounts<-"~/rnaseq_rscripts/outputs/mergedTables/AllHTSeqCounts_subset.txt"
geneIDColCounts<-"Ens_ID" #Specify here the name of the column which contains the gene IDs (they need to be unique).
#For the DESeq2 analysis you need to specify a factor on which you want to do the analysis:
factor<-"Line" #This needs to be a name of a column of the samplesPlan file.

#Optional
tableWithAnnotations<-"~/rnaseq_rscripts/outputs/mergedTables/AllCufflinks_Simplified_norm.txt" #This can be table from cufflinks or cuffdiff or Biomart to annotate genes. You will need to choose a file with at least one column with the Ensembl Gene IDs.
geneIDColInAnnotations<-"gene_id"  #Specify here the name of the column which contains the gene IDs (it must match with the content of the geneID from the table with counts).
#You can also provide a gtf:
gtfFile<-"~/rnaseq_rscripts/data/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.94_ExonsOnly_UCSC_withIsland3.gtf.gz"
changeTest<-F #Default test is Wald but you can change to likelihood ratio test (LRT) with reduced formula ~1. Put F to keep Wald and put T to use LRT.
outputDESeqTable<-"~/rnaseq_rscripts/outputs/DESeq2/DESeq2AnalysisDFLPerGeno.txt"
outputSignificantTable<-T #If you want to have another table with only significant genes abs(l2FC)>1.5 and corrected p-value<0.05
