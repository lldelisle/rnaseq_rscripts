### Required for all steps ###
RNAseqFunctionPath<-"~/Dropbox/scripts/toShare/RNAseqNewVersion/RNAseqFunctions_v05.R"
samplesPlan<-"/run/media/ldelisle/LD_2/RNAseqMarian/samplesPlanWithPath.txt" #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 1 - MERGE TABLES ###
#If the merged tables are not already generated:
outputFolderForStep1<-"/run/media/ldelisle/LD_2/RNAseqMarian/mergedTables/"
#Needed for DESeq2:
mergeCounts<-T #Do you want to merge counts? T=yes F or commented=no
#Optional: subset the count table
subsetCounts<-T #Do you want to remove some genes from the count table
#initialTableWithCount<-"/run/media/ldelisle/LD_2/RNAseqMarian/mergedTables/AllHTSeqCounts.txt" #If the table with counts have already been generated and you just want to remove some genes.
#geneIDColInInitialTable<-"Ens_ID" #If you provide the initialTableWithCount you need to provide the name of the column with the ensembl id.
genesToRmFromCounts<-"~/Dropbox/scripts/toShare/RNAseqNewVersion/MTgenes.txt" #List of genes id to remove (one per line with no header).
#Optional:
mergeFPKM<-T
oneLinePerEnsemblID<-T #By default cufflinks split the transcripts which do not overlap in different locus and so different lines, put T if you want to sum the FPKM for non overlapping transcripts (put F if not).
normFPKMWithAnoukMethod<-T #Anouk method=Genes that have the less variable rank should have the same expression.
#initialTableWithFPKM<-"/run/media/ldelisle/LD_2/RNAseqMarian/mergedTables/AllCufflinks_Simplified.txt" #If the table with FPKM have already been generated and you just want to normalize it.
chrToRemoveBeforeNormWithAnoukMethod<-c("chrM") #Usually, it is recommanded to remove mitochondrial genes before doing the normalization. In some cases, it can also be useful to remove the sex chromosomes (put c("chrX","chrY","chrM")). If you do not want to remove any gene put NA or comment the line.
nbOfGenesWithAnoukMethod<-1000 #Default is 1000, you can change here.
keepGenesUsedForNorm<-F #If you want to keep the genes used in the normalization from Anouk, they will be written in a file.

