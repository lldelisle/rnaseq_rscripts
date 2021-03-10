# Documentation

## Inputs
### From fastq to counts/FPKM
You just got your RNAseq samples sequenced, you have your fastq.
I recommand you to map them with STAR with our custom gtf from Ensembl (see repository toolboxForMutantAndWTGenomes) and getting the gene count from it. You can also use cufflinks to generate FPKM.
You should have one count file per sample and one FPKM file per sample.

For Duboule lab members, use last workflow in galaxy with STAR.

### Other possible inputs
- One table where each row is a gene and each column is a sample. If they are counts you can use (step2), if they are FPKM or like it, you can use step3.
- One table from a differential expression where each row is a gene and with at least 2 columns (padj and log2FoldChange), you can use step4 for Volcano plot, if you also have a baseMean column you can also do MA Plot.

### Samplesplan
The samplesplan is a tabulated file which have one line per sample and have multiple information that can be used for plots and analysis.

The simplest samplesplan has only one column: one with the samples names called 'sample'.
If you want to do a differencial analysis, you need another one with the name of the group on which you want to do the differential analysis.

You can add in this samplesplan multiple information that you can then use to annotate your samples on plots.

In addition, if you want to use your individual counts/FPKM files, you need to add 2 columns named htseq_count_file and cufflinks_file which should contains the full path of each file.

An example of samplesplan is available in the repository.

An easy way to do your samplesplan is:
- Open the samplesplan provided as an example (in the example folder) with excel.
- Modify the columns to fit with your experiment, add or remove extra rows, extra columns.
- Save it as tabulated delimited file.

**Warning**: Each sample name should begin with a letter, Each column name also...

### Config files
Config files are R files which will be sourced by the script to get the value of the different parameters. So, it is a R syntax: `parameter <- value` or `parameter = value`.

## step1
### Goal
The first script has 4 parts:
- Merge the tables from HTSeq-count.
- Remove from this table genes you do not want to include in the analysis (often MT genes).
- Merge the tables from Cufflinks.
- Normalize this table with a median scaling procedure method (see [Brawand, et al 2011](https://doi.org/10.1038/nature10532)). As this method was introduced to us by Anouk, we also call it AnoukMethod.

### Prepare the config file
I highly recommand to use the template configFileRNAseq_step1.R in the example folder.
- RNAseqFunctionPath: put the path for the script `RNAseqFunctions.R`
- samplesPlan: put the path for the samples plan file which correspond to your experiment.
- (optionnal) outputFolderForStep1: the folder where the merged tables will be (if not provided, the folder of the samples plan is used).
- mergeCounts: put T if you want to merge the counts.
- subsetCounts: put T if you want to remove from the merged table of counts some genes. Then you will need:
    - genesToRmFromCounts: the path for the list of ensembl ID of genes you want to exclude from the table. If you want to remove genes from a whole chromosome, you can use the script getGeneListFromChrAndGTF.R available in the [toolBoxForMutantAndWTGenomes repository](https://github.com/lldelisle/toolBoxForMutantAndWTGenomes).
    - If you had put mergeCounts<-F, you need to provide the path of the table with all counts in initialTableWithCount as well as the geneIDColInInitialTable which is the column name where you have the genes ids.
- mergeFPKM: put T if you want to merge the FPKM values from cufflinks. Then you will need:
    - oneLinePerEnsemblID (By default cufflinks split the transcripts which do not overlap in different locus and so different lines, put T if you want to sum the FPKM for non overlapping transcripts (put F if not).)
- normFPKMWithAnoukMethod: put T if you want to renormalize the FPKM (FPKM are already normalized to million mapped reads) with Anouk method: Genes that have the less variable rank should have the same expression. Then you will need:
    - chrToRemoveBeforeNormWithAnoukMethod: the list in R format (c("chrA","chrB")) of chromosome you want to remove from the normalization process. Warning, these genes will still be part of the output file.
    - nbOfGenesWithAnoukMethod: the number of genes you want to use to evaluate the correction factor which should be applied.
    - keepGenesUsedForNorm: T if you want to write in a file the genes which have been considered as house-keeping regarding the rank.
    - If you had mergeFPKM<-F, you need to provide the path of the table with all FPKM in initialTableWithFPKM.

### Launch it
Either you use command line:
```
Rscript step1-generateTables.R example/configFileRNAseq_step1.R
```

Or you open step1-generateTables.R in RStudio, click on the source button (top right).
A window will be opened and you will need to choose the config file.
You will see messages or errors.
When you see back the `>` symbol. This means it is done.

## step2
### Goal
The second script do the differential analysis with DESeq2.

### Prepare the config file
I highly recommand to use the template configFileRNAseq_step2_DFL.R in the example folder.
- RNAseqFunctionPath: put the path for the script `RNAseqFunctions.R`
- samplesPlan: put the path for the samples plan file which correspond to your experiment (which can be smaller than the original samplesPlan).
- tableWithCounts: put the path for the table with all counts (one column per sample called by the name of the sample and a column for gene identification). You can use AllHTSeqCounts.txt or AllHTSeqCounts_subset.txt if you want to exclude some genes/chr from analysis.
- geneIDColCounts: If the table was not provided by step1. Put here the name of the column for gene identification.
- factor: put the name of the column of the samples plan you want to use to do the analysis.
- changeTest: Default test is Wald but you can change to likelihood ratio test (LRT) with reduced formula ~1. Put F to keep Wald and put T to use LRT.
- outputDESeqTable: put the path where you want to have the output table of the script.
- outputSignificantTable: If you want to have another table with only significant genes abs(l2FC)>1.5 and corrected p-value<0.05 put T else put F.
- By default, the output table contains only the ensembl id to identify genes. It can be useful to add more information. To do so, you have different options:
    - gtfFile: provide the path for a gtf file, you will get coordinates, gene_biotype and gene_name.
    - tableWithAnnotations: provide the path for a table with at least one column with the ensembl ids. It can be the FPKM table or a table from biomart.
    - geneIDColInAnnotations: If you provided a tableWithAnnotations, put here the name of the column which contains the ensembl ids.

### Launch it
Either you use command line:
```
Rscript step2-DESeq2.R example/configFileRNAseq_step2_DFL.R
```

Or you open step2-DESeq2.R in RStudio, click on the source button (top right).
A window will be opened and you will need to choose the config file.
You will see messages or errors.
When you see back the `>` symbol. This means it is done.

### Description of output
In the output table you have:
- One column for ensembl ID
- (Optional) Different columns coming from your annotation file(s).
- One column per sample with the count of reads per gene normalized per DESeq2.
- The baseMean which is the mean expression over all samples
- The log2 fold change
- The log2 fold change Standard Error
- The stat value indiquate how far the difference of expression between the two groups is from the variability of this expression (in the Wald test). This value is used to calculate the p-value.
- The p-value
- The p-value adjusted for multiple tests with the Benjamini \& Hochberg method.
 
It is ordered by increasing adjusted p-value.

## step3
### Goal
The third script do three plots: PCA, Clustering and gene expession.

### Prepare the config file
I highly recommand to use the template configFileRNAseq_step3.R in the example folder.

- RNAseqFunctionPath: put the path for the last version of the script RNAseqFunctions
- samplesPlan: put the path for the samples plan file which correspond to your experiment
- tableWithNormalizedExpression: put the path for the table with either FPKM norm values or the count norm values obtained after DESeq2. If you want to exclude from the analyses the chromosomes you removed to run DESeq2, use the table from DESeq2 (you can still use FPKM if you chose as annotation file the merged table of FPKM). We expect that the name of the columns containing the expression values match the names of the samples if it is normalized counts and match FPKM_ followed by the names of the samples if it is FPKM.
- useFPKM: If your table contains columns with FPKM_ followed by the name of the samples instead of the name of your samples, put T (for example, if you are using the table from DESeq2 but you want to use FPKM values instead of counts). Otherwise, put F.
- outputFolder: put the path of the folder in which should be put the plots and the tables.
- usePng: the possible output formats are pdf (put F) or png (put T).
- fixedColors: If you want to keep a color code for all plots. Put here a list, each item of the list has the same of a column of your samplesplan and is a named vector (see configFile example).
- restrictToNMoreVariantGenes: if you want to do PCA and/or sample clustering only with the most variant genes. You can put here a number (for example, in DESeq2, they suggest to use 500).
- For PCA:
    - nbOfPC: put here the last component of the PCA to look at. Put 0 if you do not want to do any,  1=Only look at first component, 2=look at the 2 first etc...
    - PCA1D and PCA2D: put here the list of the attributes you want to see changing depending on the categories in the PCA1D (one component at the time) or PCA2D (2 components at the time).
    Please look at the example provided and test, to see what you can change.
- For clustering:
    - plotMatrixAndClustering: Do you want to perform a correlation matrix and clustering (T=yes, F=no). The factors used in PCA1D and PCA2D will be used to annotate the clustering.
- For Genes:
    - fileWithGenes: The path of the gene list. This need to be a text file where the first line is the name of the column in the file with expression then the gene names/IDs need to be one per line.
    - geneIDToAdd: By default, the title of the plot is the id provided in the fileWithGenes but you can add a meaning full name like gene_short_name if it is provided in the tableWithNormalizedExpression.
    - useLogExpression: By default, the values of expression plotted are log2(1+expression) (when T, if F the raw expression will be plotted.)
    - useSameYmaxForAllGenes: By default, each gene is plotted on an adjusted scale. If useSameYmaxForAllGenes is T, all genes will be plotted with the same y axis.
    - xaxisForGenes: By default, each sample is plotted on a seperate x but you can group them using one of the column of your samples plan. Put here the name of the column to use.
    - plotGenesPara: Put here a list like for PCA2D. More explainations in the configFile.
    - doNotPlotGeneByGene: If you do not want to have one plot per gene but just the heatmap put it to T.
    - addGlobalHeatmap: If you want to have a global heatmap with the genes provided in the list (the annotation used are XaxisForGenes and the one used in plotGenesPara).
    - keepGeneOrder: By default, genes are clustered by euclidean distance and complete clustering. If you want to keep the original order. Put keepGeneOrder to T.
    - clusterSamples: By default, samples are not clustered but you can change it.

### Launch it
Either you use command line:
```
Rscript step3-graphClusteringPCAGenes.R example/configFileRNAseq_step3.R
```

Or you open step3-graphClusteringPCAGenes.R in RStudio, click on the source button (top right).
A window will be opened and you will need to choose the config file.
You will see messages or errors.
When you see back the `>` symbol. This means it is done.

### Description of output
Of course you will have files with the plots;/tables you asked for but in addition, a file called the name of the configFile followed by _command lineLaunched.R is created.
It contains all commands launched to produce the plots and tables.
You can open this new file in RStudio and modify it as you wish especially for ggplot2 commands and then source it or run it line per line.

## step4
### Goal
The fourth script do two plots: Volcano and MAP.

### Prepare the config file
I highly recommand to use the template configFileRNAseq_step4DFL_hoxdAndaround.R in the example folder.
- RNAseqFunctionPath: put the path for the last version of the script RNAseqFunctions
- tableWithResultsOfDifferentialExpression: put the path for the table with results of DESeq2 or any table with at least 2 columns (padj and log2FoldChange) for Volcano and 3 columns (padj, log2FoldChange and baseMean) for MAP.
- outputFolder: put the path of the folder in which should be put the plots and the tables.
- usePng: the possible output formats are pdf (put F) or png (put T).
- maxPAdj: Put the maximum value of ajusted p-value to be significant (default 0.05)
- minLFC: Put the minimum log2 fold change to be significant (default 1.5)
- colOfNonSignificant: Put here the color you want for points which are not significant. The list of colors with names is available [here](http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf).
- colOfSignificant: Put here the color you want for points which are significant.
- click: Do you want to be able to click on the plot to be able to identify some genes. (T=yes, F=no).
- if click<-T:
    - geneID: provide here the name of the column which can be used to label the genes.
- colOfCircle: you can choose to make a circle around each selected gene, you can choose the color of it.
- fileWithGenes: You can provide the path to a list of gene you want to see on the volcano and MAP plots. This need to be a text file where the first line is the name of the column in the file with expression then the gene names/IDs need to be one per line.
- For Volcano:
    - maxYVolcano: If you want to zoom on bigger p-values because you have some genes with very low p-values, you can put here a value to restrict the plot to p-values higher than $10^{-maxYVolcano}$. Put NA or comment line if you do not need to zoom the initial plot.
- For MAP:
    - ylimMAP: If you want to zoom on smaller log2 fold changes because you have some genes with very high log2 fold changes, you can put here minimum and maxiumum log2 fold change values to restrict the plot to these values. Put NA or comment line if you do not need to zoom the initial plot, put c(minValue,maxValue) if you want to restrict.

### Launch it
Either you use command line (but you cannot click on the graphs):
```
Rscript step4-graphVolcanoAndMAP.R example/configFileRNAseq_step4DFL_hoxdAndaround.R
```

Or you open step4-graphVolcanoAndMAP.R in RStudio, click on the source button (top right).
A window will be opened and you will need to choose the config file.
You will see messages or errors.
If you put click <- T, then a graph will be display in the graph part of RStudio, you will be able to click on dots and when you finished, press ESC.
When you see back the `>` symbol. This means it is done.

## Additional info
### Gene contribution to PCA

- The PCA is a simplification of your data based on the variance (the mean level of expression is totally omitted): instead of having one dimension per gene you have one dimension per principal component. By default in R you have 8 components but the most important ones are usually the first, the second and sometimes the third.
- The most important things in a PCA are the coordinates of each sample in this new space (PC1, PC2 and sometimes PC3).
- You can also take advantage of this transformation to highlight genes that contributed the most to each of the component. For example, imagine you have 2 types of tissue and 2 genotypes. If the PC1 spread among tissue and PC2 among genotype, you may be interested to know which genes are hidden in the PC1 and which are hidden in the PC2.
- This information is in the geneContributionToPCA.txt table. For each gene you have the proportion of its contribution to each component. If you sort them for PC1 for example the genes with highest values are the one that have strong difference between the tissues and for PC2 difference between genotypes.

### Correlation matrix

- Distance:
    - By default:
        The matrix is based on the Euclidean distance or spearman correlation (1-cor).
    - To change it:
        Open the file created after the third script have been launched with all command lines. Look for the line containing: `sampleDists <- dist(t(rldata),method="euclidean")` and replace by `sampleDists <- dist(t(rldata),method="blablabla")`  (see [here](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/dist.html) for the possible methods).
- Clustering:
    - By default:
        The clustering method of the first cluster is "complete" this mean that the distance between two clusters is the maximum distance between individuals. Then it is ward.
    - To change it:
        Open the file created after the third script have been launched with all command lines.
        Look for the line containing: clustering_method="complete" and replace by clustering_method="blablabla" (see [here](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html) for the possible methods).