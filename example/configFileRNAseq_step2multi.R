### Required for all steps ###
RNAseqFunctionPath <- "~/rnaseq_rscripts/RNAseqFunctions.R"
# This file should be a tabulated file with at least one column called
# 'sample'.
# At this step it is super important to put the samples reference for your comparison before the others.
samplesPlan <- "~/rnaseq_rscripts/example/samplesPlan.txt"

#### STEP 2 - DESEQ 2 ANALYSIS ####
# Required
tableWithCounts <- "~/rnaseq_rscripts/outputs/mergedTables/AllHTSeqCounts_subset.txt"
# Specify here the name of the column which contains the gene IDs (they need to
# be unique).
geneIDColCounts <- "Ens_ID"
# For the DESeq2 analysis you need to specify a factor on which
# you want to do the analysis:
# This needs to be a name of a column of the samplesPlan file.
# The reference will be defined as the first value of this column
# among the subsetted samples.
# Here you want to do multiple analyses at once.
# You describe all what you want into the variable 'all.anayses':
# list(factorToStudy = list(loopingVariable = list(subsetting)))

# First super simple example:
# You want to compare all to a reference:
# all.analyses <- list("Condition" = list("Condition" = list()))
# The reference will be the first value in the "Condition" column.
# All pairwise comparisons between the different values
# of the "Condition" column and
# the first value in the "Condition" column will be computed.

# Second example:
# You want to compare mutant to WT among DFL samples
# And you want to compare PFL to DFL amount wt samples
all.analyses <- list(
  "Line" = list("Line" = list("Tissue" = "DFL")),
  "Tissue" = list("Tissue" = list("Line" = "Wt"))
)
# How the software will imprete this:
# First,
#   using Line as factor,
#   for each Line which is not the reference = first in samplesplan
#   (here there is only one but there could be more),
#   using only samples where Tissue is DFL.
#
# Then,
#   using Tissue as factor,
#   for each Tissue which is not the reference = first in samplesplan
#   (here only one),
#   using only samples where Line is Wt.

# Third example:
# Compare mutant to wt but in each tissue separately:
# all.analyses <- list(
#   "Line" = list(
#     "Line" = list("Tissue" = "Limbs"),
#     "Line" = list("Tissue" = "Trunk")
#   )
# )
# In the given example, it will use the factor 'Line' and do:
# 1.
#    For each Line which is not the reference (first one in samples plan),
#    comparison to the reference
#    using only samples where Tissue is Limbs.
# 2.
#    For each Line which is not the reference (first one in samples plan),
#    comparison to the reference
#    using only samples where Tissue is Trunk.


# Optional
# If you want to add covariate like library type, experimental batch etc...
# The names must match columns of the samples plan
# Uncomment:
# covariates <- list("covariate1", "covariate2")
# For example:
# covariates <- list("Replicate")

# If you want to filter genes based on a minimal FPKM mean value in at least one group
# Uncomment and update the path:
# tableWithFPKM <- "~/rnaseq_rscripts/outputs/mergedTables/AllCufflinks_Simplified.txt"
# min.mean.FPKM <- 1

# This can be table from cufflinks or cuffdiff or Biomart to annotate genes.
# This will be only added to the 2 final output files (summary.txt and all_summary.txt)
# You will need to choose a file with at least one column with the Ensembl Gene IDs.
tableWithAnnotations <- "~/rnaseq_rscripts/outputs/mergedTables/AllCufflinks_Simplified_norm.txt"
# Specify here the name of the column which contains the gene IDs (it must
# match with the content of the geneID from the table with counts).
geneIDColInAnnotations <- "gene_id"
# You can also provide a gtf:
gtfFile <- "~/rnaseq_rscripts/data/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.94_ExonsOnly_UCSC_withIsland3.gtf.gz"

# A directory with outputs
pathForDESeq2 <- "~/rnaseq_rscripts/outputs/DESeq2_multi/"
# Set the threshold for significance here
log2FC.threshold <- log2(1.5)
