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
      stop(
        "Ens_ID is not part of the column names and geneIDColCounts:", geneIDColCounts,
        "neither."
      )
    } else {
      colnames(htseqCounts)[colnames(htseqCounts) == geneIDColCounts] <- "Ens_ID"
    }
  }
}

rownames(htseqCounts) <- htseqCounts$Ens_ID

# The samplesplan may contain more values than the htseqCounts
sampleNamesWithValues <- intersect(samplesPlanDF$sample, colnames(htseqCounts))
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

if (sum(is.na(htseqCounts)) > 0) {
  cat("There are NAs in the count table. They will be replaced by 0.\n")
  htseqCounts[is.na(htseqCounts)] <- 0
}

samples.plan.df <- simplifyDF(samplesPlanDF, sampleNamesWithValues)

big.table.fn.long <- "summary_long.txt"

### DESEQ2 ANALYSIS ###

# Prepare a big table with the results of all DESeq2
big.annot <- data.frame(Ens_ID = htseqCounts$Ens_ID)
big.annot2 <- NULL

if (!dir.exists(pathForDESeq2)) {
  dir.create(pathForDESeq2, recursive = TRUE)
}
# I will loop over the list all.analyses:
# First is the factorToStudy
for (factorToStudy in names(all.analyses)) {
  print("FACTOR TO STUDY")
  print(factorToStudy)
  # Second indicates the looping variable
  for (i in 1:length(all.analyses[[factorToStudy]])) {
    loopingVariable <- names(all.analyses[[factorToStudy]][i])
    print("LOOPING VARIABLE")
    print(loopingVariable)
    # Third indicates the subsetting
    pre.samples.plan <- samples.plan.df
    subsetting.list <- all.analyses[[factorToStudy]][i][[loopingVariable]]
    subsetting.name <- ""
    for (rn in names(subsetting.list)) {
      print("SUBSET")
      print(rn)
      pre.samples.plan <- pre.samples.plan[pre.samples.plan[, rn] %in% subsetting.list[[rn]], ]
      subsetting.name <- paste0(subsetting.name, paste(subsetting.list[[rn]], collapse = "or"), "_")
    }
    looping.values <- unique(pre.samples.plan[, loopingVariable])
    ref.value <- NULL
    if (factorToStudy == loopingVariable) {
      ref.value <- intersect(levels(pre.samples.plan[, loopingVariable]), pre.samples.plan[, loopingVariable])[1]
      looping.values <- setdiff(looping.values, ref.value)
    }
    for (my.value in looping.values) {
      new.samples.plan <- pre.samples.plan[pre.samples.plan[, loopingVariable] %in% c(ref.value, my.value), ]
      # Drop levels for factorToStudy
      new.samples.plan[, factorToStudy] <- factor(new.samples.plan[, factorToStudy],
        levels = intersect(
          levels(new.samples.plan[, factorToStudy]),
          unique(new.samples.plan[, factorToStudy])
        )
      )
      base.filename <- paste0(factorToStudy, "_", subsetting.name, my.value)
      if (factorToStudy == loopingVariable) {
        base.filename <- paste0(base.filename, "vs", ref.value, "_")
      } else {
        base.filename <- paste0(base.filename, "_")
      }
      print(base.filename)
      # Run or read DESeq2 results with Wald test threshold of FC at 1.5
      if (!file.exists(file.path(pathForDESeq2, paste0(base.filename, "DESeq2significant.txt")))) {
        print(new.samples.plan)
        deseqAnaWithCovariates(htseqCounts, factorToStudy, NULL,
          file.path(pathForDESeq2, base.filename),
          new.samples.plan,
          LRT = FALSE,
          lfcT = log2FC.threshold,
          writeRLOG = FALSE,
          gene_id = "Ens_ID"
        )
        # theta = c(0.15, 0.99))
      } else {
        print("Exists")
      }
      # Add results to the dataframe
      all.res <- read.delim(file.path(pathForDESeq2, paste0(base.filename, "DESeq2Results.txt")))
      rownames(all.res) <- all.res[, "Ens_ID"]
      # Add results to the dataframe
      big.annot[, paste0(base.filename, "l2fc")] <-
        all.res$log2FoldChange[match(big.annot[, "Ens_ID"], all.res[, "Ens_ID"])]
      big.annot[, paste0(base.filename, "padj")] <- all.res$padj[match(big.annot[, "Ens_ID"], all.res[, "Ens_ID"])]
      big.annot[, paste0(base.filename, "signif")] <-
        with(
          all.res[match(big.annot[, "Ens_ID"], all.res[, "Ens_ID"]), ],
          !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold
        )
      tail(big.annot)
      all.res.fmt <- subset(
        all.res,
        select = intersect(
          c("Ens_ID", "gene_short_name", "baseMean", "log2FoldChange", "padj"),
          colnames(all.res)
        )
      )
      all.res.fmt$factor <- factorToStudy
      all.res.fmt$subsetting <- subsetting.name
      all.res.fmt$value <- my.value
      all.res.fmt$ref.value <- ref.value
      if (!is.null(ref.value) & !is.null(big.annot2) & !"ref.value" %in% colnames(big.annot2)) {
        big.annot2$ref.value <- NA
      }
      if (is.null(ref.value) & "ref.value" %in% colnames(big.annot2)) {
        all.res.fmt$ref.value <- NA
      }
      big.annot2 <- rbind(big.annot2, all.res.fmt)
      tail(big.annot2)
    }
  }
}

### Annotation addition ###

if (exists("tableWithAnnotations")) {
  if (file.exists(tableWithAnnotations)) {
    ann <- read.delim(tableWithAnnotations, check.names = FALSE)
    if (exists("geneIDColInAnnotations")) {
      if (!geneIDColInAnnotations %in% colnames(ann)) {
        cat(
          "The geneIDColInAnnotations provided: ", geneIDColInAnnotations,
          " is not part of the column names of the annotation file. Unable to match the annotation file with the results of DESeq2. The annotation file will not be added.\n"
        )
        geneIDColInAnnotations <- NA
      }
    } else if ("Ens_ID" %in% colnames(ann)) {
      geneIDColInAnnotations <- "Ens_ID"
    } else {
      cat("The geneIDColInAnnotations is not defined in the config file. Unable to match the annotation file with the results of DESeq2. The annotation file will not be added.\n")
      geneIDColInAnnotations <- NA
    }
    if (!is.na(geneIDColInAnnotations)) {
      annot.df <- ann
      colnames(annot.df)[colnames(annot.df) == geneIDColInAnnotations] <- "Ens_ID"
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
    start <- aggregate(list(start = gtf$start - 1), by = list(gene_id = gtf$gene_id), FUN = min)
    end <- aggregate(list(end = gtf$end), by = list(gene_id = gtf$gene_id), FUN = max)
    gtf.columns <- intersect(
      c("gene_id", "seqid", "strand", "gene_name", "gene_biotype", "gene_type"),
      colnames(gtf)
    )
    gtf <- unique(gtf[, gtf.columns])

    annot.gtf <- merge(gtf, start)
    annot.gtf <- merge(annot.gtf, end)
    annot.gtf <- annot.gtf[, c(gtf.columns, "start", "end")]
    colnames(annot.gtf)[1] <- "Ens_ID"
    if (exists("annot.df")) {
      annot.df <- merge(annot.gtf, annot.df, all = TRUE)
    } else {
      annot.df <- annot.gtf
    }
  } else {
    cat("The gtf file specified:", gtfFile, " does not exists. It will not be added.\n")
  }
}


if (exists("annot.df")) {
  big.annot <- merge(big.annot, annot.df, all.x = TRUE)
  big.annot2 <- merge(big.annot2, annot.df, all.x = TRUE)
}

write.table(
  big.annot,
  file.path(pathForDESeq2, "summary.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  big.annot2,
  file.path(pathForDESeq2, big.table.fn.long),
  sep = "\t", row.names = FALSE, quote = FALSE
)
