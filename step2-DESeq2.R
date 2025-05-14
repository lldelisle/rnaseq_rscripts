options(stringsAsFactors = F)
rm(list = ls())

if (any(!c("colorspace", "DESeq2", "tools") %in% installed.packages())) {
  if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
  }
  devtools::install_github("lldelisle/usefulLDfunctions")
  library(usefulLDfunctions)
  safelyLoadAPackageInCRANorBioconductor("colorspace")
  safelyLoadAPackageInCRANorBioconductor("DESeq2")
  safelyLoadAPackageInCRANorBioconductor("tools")
} else {
  library(colorspace)
  library(DESeq2)
  library(tools)
}
if (length(commandArgs(TRUE)) > 0) {
  f <- commandArgs(TRUE)[1]
} else {
  cat("Select the config file.\n")
  # Ask for the config file
  f <- file.choose()
}
### Check the necessary options###
if (!file.exists(f)) {
  stop("This file does not exist.")
}
source(f)
if (!exists("samplesPlan")) {
  stop("The config file do not have samplesPlan definition.")
}
if (!file.exists(samplesPlan)) {
  stop("The file specified as samplesPlan does not exist:", samplesPlan)
}
samplesPlanDF <- read.delim(samplesPlan, check.names = FALSE)
if (!("sample" %in% colnames(samplesPlanDF))) {
  stop("The samplesPlan table do not contain a column called \"sample\".")
}
setwd(dirname(samplesPlan))

if (!exists("tableWithCounts")) {
  stop("The config file do not have tableWithCounts definition.")
}
if (!file.exists(tableWithCounts)) {
  stop("The file specified as tableWithCounts:", tableWithCounts, "does not exists.")
}
htseqCounts <- read.delim(tableWithCounts, check.names = FALSE)

if (!"Ens_ID" %in% colnames(htseqCounts)) {
  if (!exists("geneIDColCounts")) {
    stop("Ens_ID is not part of the column names and geneIDColCounts is not specified.")
  } else {
    if (!geneIDColCounts %in% colnames(htseqCounts)) {
      stop("Ens_ID is not part of the column names and geneIDColCounts:", geneIDColCounts,
        "neither.")
    } else {
      colnames(htseqCounts)[colnames(htseqCounts) == geneIDColCounts] <- "Ens_ID"
    }
  }
}

if (!exists("tableWithFPKM")) {
  fpkm.table <- NULL
} else if (!file.exists(tableWithFPKM)) {
  stop("The file specified as tableWithFPKM:", tableWithFPKM, "does not exists.")
} else {
  fpkm.table <- read.delim(tableWithFPKM, check.names = FALSE)
  if (!"gene_id" %in% colnames(fpkm.table)) {
    stop("gene_id is not part of the column names of the FPKM table")
  }
  rownames(fpkm.table) <- fpkm.table$gene_id
  if (!exists("min.mean.FPKM")) {
    min.mean.FPKM <- 1
  } else if (!is.numeric(min.mean.FPKM)) {
    stop("min.mean.FPKM must be numeric.")
  }
}

# The samplesplan may contain more values than the htseqCounts
sampleNamesWithValues <- intersect(colnames(htseqCounts), samplesPlanDF$sample)
if (length(sampleNamesWithValues) < 2) {
  cat("Samples in samplesplan are \n")
  cat(samplesPlanDF$sample)
  cat("\n")
  cat("The column names of the file are \n")
  cat(colnames(htseqCounts))
  cat("\n")
  stop("The samplesplan and the HTSeqcount table are not compatible")
}

if (!exists("RNAseqFunctionPath")) {
  stop("The RNAseqFunctionPath is not provided.")
} else {
  if (!file.exists(RNAseqFunctionPath)) {
    stop("The file provided in RNAseqFunctionPath:", RNAseqFunctionPath, " does not exists.")
  }
}
source(RNAseqFunctionPath)



# The htseqCounts dataframe have geneID as rownames
rownames(htseqCounts) <- htseqCounts$Ens_ID
htseqCounts[, "Ens_ID"] <- NULL

simplifiedSP <- simplifyDF(samplesPlanDF, sampleNamesWithValues)

possibleFactors <- setdiff(colnames(simplifiedSP), c("sample"))

if (length(possibleFactors) == 0) {
  cat("Will use sample as factor.\n")
  factor <- "sample"
} else if (length(possibleFactors) == 1) {
  cat("Will use ", possibleFactors, " as factor.\n")
  factor <- possibleFactors
} else if (!exists("factor")) {
  stop("To do a DESeq2 analysis, a factor (a group definition for samples) is required. Here there are multiple possibilities.")
} else {
  if (!factor %in% colnames(samplesPlanDF)) {
    stop("The factor specified in the config file:", factor, " is not part of the column names of the samplesPlan table.")
  } else if (!factor %in% colnames(simplifiedSP)) {
    stop("There is only one value for the factor ", factor, " in the samples present in the count table.")
  }
}


### DESEQ2 ANALYSIS ###

if (sum(is.na(htseqCounts)) > 0) {
  cat("There are NAs in the count table. They will be replaced by 0.\n")
  htseqCounts[is.na(htseqCounts)] <- 0
}

cmd <- paste("dds <- DESeqDataSetFromMatrix(countData = htseqCounts[,rownames(simplifiedSP)], colData = simplifiedSP, design = ~ ",
  factor, ")", sep = "")
eval(parse(text = cmd))

# Filtering:
if (is.null(fpkm.table)) {
  cat("Genes that are never expressed are removed.\n")
  dds <- dds[rowSums(counts(dds)) > 1, ]
} else {
  fpkm.means <- NULL
  for (factorValue in unique(simplifiedSP[, factor])) {
    current.samples <- rownames(
      subset(simplifiedSP,
             simplifiedSP[, factor] == factorValue)
    )
    fpkm.means <- cbind(fpkm.means, rowMeans(fpkm.table[, paste0("FPKM_", current.samples)]))
  }
  colnames(fpkm.means) <- paste0("FPKM_mean_", unique(simplifiedSP[, factor]))
  rownames(fpkm.means) <- rownames(fpkm.table)
  fpkm.means.max <- apply(fpkm.means, 1, max)
  genes.to.keep <- intersect(
    names(fpkm.means.max[fpkm.means.max > min.mean.FPKM]),
    rownames(dds)[rowSums(counts(dds)) > 1]
  )
  cat(paste("Genes with average FPKM per group below", min.mean.FPKM, " in all groups are removed.\n"))
  dds <- dds[genes.to.keep, ]
}
if (exists("changeTest")) {
  if (!is.logical(changeTest)) {
    cat("The value provided in changeTest is not logical. The default test (Wald will be used.\n)")
    changeTest <- F
  }
} else {
  changeTest <- F
}

if (changeTest) {
  dds <- DESeq(dds, test = "LRT", reduced = ~1)
} else {
  dds <- DESeq(dds)
}

res <- results(dds)
resOrdered <- res[order(res$padj), ]

### Annotation addition ###

if (exists("tableWithAnnotations")) {
  if (file.exists(tableWithAnnotations)) {
    ann <- read.delim(tableWithAnnotations, check.names = FALSE)
    if (exists("geneIDColInAnnotations")) {
      if (!geneIDColInAnnotations %in% colnames(ann)) {
        cat("The geneIDColInAnnotations provided: ", geneIDColInAnnotations,
          " is not part of the column names of the annotation file. Unable to match the annotation file with the results of DESeq2. The annotation file will not be added.\n")
        geneIDColInAnnotations <- NA
      }
    } else if ("gene_id" %in% colnames(ann)) {
      geneIDColInAnnotations <- "gene_id"
    } else if ("Ens_ID" %in% colnames(ann)) {
      geneIDColInAnnotations <- "Ens_ID"
    } else {
      cat("The geneIDColInAnnotations is not defined in the config file. Unable to match the annotation file with the results of DESeq2. The annotation file will not be added.\n")
      geneIDColInAnnotations <- NA
    }
    if (!is.na(geneIDColInAnnotations)) {
      annot.df <- ann[match(rownames(resOrdered), ann[, geneIDColInAnnotations]),
        ]
    }
  } else {
    cat("The annotation file specified:", tableWithAnnotations, " does not exists. It will not be added.\n")
  }
}

if (exists("gtfFile")) {
  if (file.exists(gtfFile)) {
    if (!"rtracklayer" %in% installed.packages()) {
      if (!"devtools" %in% installed.packages()) {
        install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
      }
      devtools::install_github("lldelisle/usefulLDfunctions")
      library(usefulLDfunctions)
      safelyLoadAPackageInCRANorBioconductor("rtracklayer")
    } else {
      library(rtracklayer)
    }
    cat("Reading gtf file...")
    gtf <- readGFF(gtfFile)
    cat("Done\n")
    start <- aggregate(gtf$start - 1, by = list(gene_id = gtf$gene_id), FUN = min)
    end <- aggregate(gtf$end, by = list(gene_id = gtf$gene_id), FUN = max)
    gtf.columns <- intersect(
      c("seqid", "strand", "gene_id", "gene_name", "gene_biotype", "gene_type"),
      colnames(gtf)
    )
    gtf <- unique(gtf[, gtf.columns])

    annot.gtf <- gtf[match(rownames(resOrdered), gtf$gene_id), ]
    annot.gtf$start <- start$x[match(rownames(resOrdered), start$gene_id)]
    annot.gtf$end <- end$x[match(rownames(resOrdered), end$gene_id)]
    annot.gtf <- annot.gtf[, c(gtf.columns, "start", "end")]
    if (exists("annot.df")) {
      annot.df <- cbind(annot.gtf, annot.df[, setdiff(colnames(annot.df), geneIDColInAnnotations)])
    } else {
      annot.df <- annot.gtf
    }
  } else {
    cat("The gtf file specified:", gtfFile, " does not exists. It will not be added.\n")
  }
}


if (exists("annot.df")) {
  resToExport <- cbind(annot.df, counts(dds, normalized = TRUE)[rownames(resOrdered),
    ], resOrdered)
} else {
  resToExport <- data.frame(Ens_ID = rownames(resOrdered), counts(dds, normalized = TRUE)[rownames(resOrdered),
    ], resOrdered, check.names = FALSE)
}

if (!is.null(fpkm.table)) {
  resToExport <- cbind(
    resToExport,
    fpkm.means[rownames(resOrdered), ]
  )
}

if (exists("outputDESeqTable")) {
  dir.create(dirname(outputDESeqTable), recursive = T, showWarnings = F)
} else {
  outputDESeqTable <- paste0(dirname(samplesPlan), "/DESeq2AnalysisForFactor",
    factor, ".txt")
}
write.table(resToExport, outputDESeqTable, sep = "\t", row.names = F, quote = F)
cat("The result table is written in ", outputDESeqTable, "\n")

if (exists("outputSignificantTable")) {
  if (is.logical(outputSignificantTable)) {
    if (outputSignificantTable) {
      outputSignificantTableFN <- paste0(dirname(outputDESeqTable), "/Significant_",
        basename(outputDESeqTable))
      write.table(subset(resToExport, padj < 0.05 & abs(log2FoldChange) > 1.5),
        outputSignificantTableFN, sep = "\t", row.names = F, quote = F)
      cat("The result table with significants is written in ", outputSignificantTableFN,
        "\n")
    }
  }
}
