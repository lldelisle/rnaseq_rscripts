### Required for all steps ###
RNAseqFunctionPath<-"~/Dropbox/scripts/toShare/RNAseqNewVersion/RNAseqFunctions_v05.R"
samplesPlan<-"/run/media/ldelisle/LD_2/RNAseqMarian/samplesPlanWithPath.txt" #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 2 - DESEQ 2 ANALYSIS ###
#Required
tableWithCounts<-"/run/media/ldelisle/LD_2/RNAseqMarian/mergedTables/AllHTSeqCounts_subset.txt"
geneIDColCounts<-"Ens_ID" #Specify here the name of the column which contains the gene IDs (they need to be unique).
#For the DESeq2 analysis you need to specify a factor on which you want to do the analysis:
factor<-"tissue" #This needs to be a name of a column of the samplesPlan file.

#Optional
tableWithAnnotations<-"/run/media/ldelisle/LD_2/RNAseqMarian/mergedTables/AllCufflinks_Simplified_norm.txt" #This can be table from cufflinks or cuffdiff or Biomart to annotate genes. You will need to choose a file with at least one column with the Ensembl Gene IDs.
geneIDColInAnnotations<-"gene_id"  #Specify here the name of the column which contains the gene IDs (it must match with the content of the geneID from the table with counts).
#You can also provide a gtf:
gtfFile<-"/run/media/ldelisle/LD_2/RNAseqMarian/FilteredTranscriptsOfMus_musculus.GRCm38.90_ExonsOnly_UCSC.gtf"
changeTest<-F #Default test is Wald but you can change to likelihood ratio test (LRT) with reduced formula ~1. Put F to keep Wald and put T to use LRT.
outputDESeqTable<-"/run/media/ldelisle/LD_2/RNAseqMarian/DESeq2/DESeq2AnalysisForFactorTissue.txt"
outputSignificantTable<-T #If you want to have another table with only significant genes abs(l2FC)>1.5 and corrected p-value<0.05
