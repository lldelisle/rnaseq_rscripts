options(stringsAsFactors=F)
rm(list=ls())

library(tools)

if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggrepel")

if(length(commandArgs(TRUE))>0){
  f<-commandArgs(TRUE)[1]
} else{
  cat("Select the config file.\n")
  #Ask for the config file
  f <- file.choose()
}
###Check the common necessary options###
if(!file.exists(f)){
  stop("This file does not exist.")
}
source(f)

if(length(commandArgs(TRUE))>0){
  click<-F
}

if(!exists("tableWithResultsOfDifferentialExpression")){
  stop("The config file do not have tableWithResultsOfDifferentialExpression definition.")
}
if(!file.exists(tableWithResultsOfDifferentialExpression)){
  stop("The file specified as tableWithResultsOfDifferentialExpression:",tableWithResultsOfDifferentialExpression,"does not exists.")
}

df<-read.delim(tableWithResultsOfDifferentialExpression)

if (!all(c("padj","log2FoldChange") %in% colnames(df))) {
      stop("The file do not contains padj and log2FoldChange column names.")
}

if(!exists("RNAseqFunctionPath")){
    stop("The RNAseqFunctionPath is not provided.")
} else {
  if(!file.exists(RNAseqFunctionPath)){
    stop("The file provided in RNAseqFunctionPath:",RNAseqFunctionPath," does not exists.")
  }
}

source(RNAseqFunctionPath)

if(exists("maxPAdj")){
  if(!is.numeric(maxPAdj)){
    cat("maxPAdj is not numeric. 0.05 will be used.\n")
    maxPAdj<-0.05
  }
}else {
  cat("maxPAdj is not defined. 0.05 will be used.\n")
  maxPAdj<-0.05
}

if(exists("minLFC")){
  if(!is.numeric(minLFC)){
    cat("minLFC is not numeric. 1.5 will be used.\n")
    minLFC<-1.5
  }
}else {
  cat("minLFC is not defined. 1.5 will be used.\n")
    minLFC<-1.5
}

if(exists("colOfNonSignificant")){
  if(!isValidColor(colOfNonSignificant)){
    cat("colOfNonSignificant cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\n grey will be used.\n")
    colOfNonSignificant<-"grey"
  }
}else {
  cat("colOfNonSignificant is not defined. grey will be used.\n")
    colOfNonSignificant<-"grey"
}

if(exists("colOfSignificant")){
  if(!isValidColor(colOfSignificant)){
    cat("colOfSignificant cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\n blue will be used.\n")
    colOfSignificant<-"blue"
  }
}else {
  cat("colOfSignificant is not defined. blue will be used.\n")
    colOfSignificant<-"blue"
}

if(!exists("outputFolder")){
  outputFolder<-paste0(dirname(samplesPlan),"/plots")
}
dir.create(outputFolder,recursive=T,showWarnings=F)

if(exists("usePng")){
  if(!is.logical(usePng)){
    cat("usePng is not logical. pdf will be used as output.\n")
    usePng<-F
  }
} else {
  usePng<-F
}

if(exists("click")){
  if(is.logical(click)){
    if(click){
      if(exists("geneID")){
        if(!geneID%in%colnames(df)){
          cat("click is set to T but the geneID is not part of the names of the columns of the tableWithResultsOfDifferentialExpression.\nNo clicking will be done.\n")
          click<-F
        }
      } else {
        cat("click is set to T but the geneID is not part of the config file.\nNo clicking will be done.\n")
        click<-F
      }
    }
  } else {
    cat("click is not logical. It will be put to F.\n")
    click<-F
  }
} else {
  click<-F
}

if(exists("fileWithGenes")){
  if(file.exists(fileWithGenes)){
    dfGene<-read.delim(fileWithGenes)
    #colOfGeneID indicates in which column the gene list should be look for
    colOfGeneID <- colnames(dfGene)[1]
    if (!(colOfGeneID %in% colnames(df))) {
      cat("The first line of the gene file does not correspond to a column in the tableWithResultsOfDifferentialExpression file.\nIt will not be used.\n")
      rm(dfGene)
    } else {
      if(exists("colOfCircle")){
        if(!isValidColor(colOfCircle)){
          cat("colOfCircle cannot be interpreted by R. Go to http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf or use rgb to create a compatible color.\n No circle will be plotted around the genes in the list.\n")
          rm(colOfCircle)
        }
      }
    }
  } else {
    cat("The file provided for fileWithGenes:",fileWithGenes,"does not exists. It will not be used.\n")
  }
}

colorDots <- rep(colOfNonSignificant,length(df$padj))
colorDots[df$padj < maxPAdj & abs(df$log2FoldChange) > minLFC] <- colOfSignificant


### Volcano ###

if(usePng){
  png(paste0(outputFolder,"/Volcano.png"))
} else {
  pdf(paste0(outputFolder,"/Volcano.pdf"),title="Volcano")
}
# Remove the genes with no padj:
sub.df <- subset(df, ! is.na(padj))
sub.colorDots <- colorDots[!is.na(df$padj)]
#the first plot is drawn with all points (volcano plot is -log10(p-val)=f(log2FC))
plot(
  sub.df$log2FoldChange,-log10(sub.df$padj),pch = 16,col = sub.colorDots,cex = 0.5,
  xlab = "log2 Fold Change",ylab = "-log10 of adjusted p-value"
)
yMaxUsed<-par('usr')[4]
dev.off()

if(exists("maxYVolcano")){
  if(!is.na(maxYVolcano)){
    if(is.numeric(maxYVolcano)){
      #the plot is redrawn with the maximum set per the user
      if(usePng){
        png(paste0(outputFolder,"/Volcano_subset.png"))
      } else {
        pdf(paste0(outputFolder,"/Volcano_subset.pdf"),title="Volcano_subset")
      }
      plot(
        sub.df$log2FoldChange,-log10(sub.df$padj),pch = 16,col = sub.colorDots,cex = 0.5,
        ylim = c(0,maxYVolcano),
        xlab = "log2 Fold Change",ylab = "-log10 of adjusted p-value"
      )
      dev.off()
      yMaxUsed<-maxYVolcano
    }
  }
}

if(click){
  plot(
    sub.df$log2FoldChange,-log10(sub.df$padj),pch = 16,col = sub.colorDots,cex = 0.5,
    ylim = c(0,yMaxUsed),
    xlab = "log2 Fold Change",ylab = "-log10 of adjusted p-value"
  )
  cat("Click on points to print the name of the gene.\nClick on ESC to stop selecting points.\n")
  idx <-
    identify(sub.df$log2FoldChange,-log10(sub.df$padj), labels = sub.df[,geneID])
  write.table(sub.df[idx,],paste0(outputFolder,"/ValuesForGenesClickedInVolcano.txt"),sep="\t",quote=F,row.names=F)
  cat("The values for the points identified are in",paste0(outputFolder,"/ValuesForGenesClickedInVolcano.txt.\n"))
  if(usePng){
    dev.copy(png, paste0(outputFolder,"/Volcano_clicked.png"))
  } else {
    dev.copy(pdf, paste0(outputFolder,"/Volcano_clicked.pdf"),title="Volcano_clicked")
  }
  dev.off()
  library(ggplot2)
  library(ggrepel)

  cmd<-paste("ggplot(data = sub.df, aes(x = log2FoldChange, y = -log10(padj))) + theme_classic() + 
  geom_point(colour = sub.colorDots, size = 1) +
  xlab(\"log2 Fold Change\") +
  ylab(\"-log10 of adjusted p-value\") +
  ylim(y=0,yMaxUsed) +
  geom_text_repel(data = sub.df[idx,] ,aes(label =",geneID,"), 
       box.padding = unit(0.45, \"lines\"))")
  if(exists("colOfCircle")){
    cmd<-paste0(cmd,"+
    geom_point(data=sub.df[idx,],
             aes(log2FoldChange, y = -log10(padj)),shape=1,col=colOfCircle)")
  }
  if(usePng){
    png(paste0(outputFolder,"/Volcano_clicked_pretty.png"))
  } else {
    pdf(paste0(outputFolder,"/Volcano_clicked_pretty.pdf"),title="Volcano_clicked_pretty")
  }
  print(eval(parse(text = cmd)))
  dev.off()
}

if(exists("dfGene")){
  library(ggplot2)
  library(ggrepel)
  cmd<-paste("ggplot(data = sub.df, aes(x = log2FoldChange, y = -log10(padj))) + theme_classic() + 
  geom_point(colour = sub.colorDots, size = 1) +
  xlab(\"log2 Fold Change\") +
  ylab(\"-log10 of adjusted p-value\") +
  ylim(y=0,yMaxUsed) +
  geom_text_repel(data = sub.df[sub.df[,colOfGeneID]%in%dfGene[,1],] ,aes(label =",colOfGeneID,"), 
       box.padding = unit(0.45, \"lines\"))")
  if(exists("colOfCircle")){
    cmd<-paste0(cmd,"+
    geom_point(data=sub.df[sub.df[,colOfGeneID]%in%dfGene[,1],],
             aes(log2FoldChange, y = -log10(padj)),shape=1,col=colOfCircle)")
  }
  if(usePng){
    png(paste0(outputFolder,"/Volcano_listOfGenes_pretty.png"))
  } else {
    pdf(paste0(outputFolder,"/Volcano_listOfGenes_pretty.pdf"),title="Volcano_listOfGenes_pretty")
  }
  print(eval(parse(text = cmd)))
  dev.off()
  write.table(sub.df[sub.df[,colOfGeneID]%in%dfGene[,1],],paste0(outputFolder,"/ValuesForGenesIndfGenes.txt"),sep="\t",quote=F,row.names=F)
}

### MAP ###

if(!"baseMean"%in%colnames(df)){
  stop("No MAP can be done because there is no column called baseMean in the tableWithResultsOfDifferentialExpression.")
}

if(usePng){
  png(paste0(outputFolder,"/MAP.png"))
} else {
  pdf(paste0(outputFolder,"/MAP.pdf"),title="MAP")
}
#the first plot is drawn with all points (MAP plot is log2FC=f(meanExpression))
plot(
  log2(1+df$baseMean),df$log2FoldChange,pch = 16,col = colorDots,cex = 0.5,
  xlab = "log2(1+Mean counts)",ylab = "log2 Fold Change"
)
yLimUsed<-par('usr')[3:4]
dev.off()

if(exists("ylimMAP")){
  if(!all(is.na(ylimMAP))){
    if(all(is.numeric(ylimMAP))&length(ylimMAP)==2){
      #the plot is redrawn with the limit set per the user
      if(usePng){
        png(paste0(outputFolder,"/MAP_subset.png"))
      } else {
        pdf(paste0(outputFolder,"/MAP_subset.pdf"),title="MAP_subset")
      }
      plot(
        log2(1+df$baseMean),df$log2FoldChange,pch = 16,col = colorDots,cex = 0.5,
        ylim = ylimMAP,
        xlab = "log2(1+Mean counts)",ylab = "log2 Fold Change"
      )
      dev.off()
      yLimUsed<-ylimMAP
    }
  }
}

if(click){
  plot(
    log2(1+df$baseMean),df$log2FoldChange,pch = 16,col = colorDots,cex = 0.5,
    ylim = yLimUsed,
    xlab = "log2(1+Mean counts)",ylab = "log2 Fold Change"
  )
  cat("Click on points to print the name of the gene.\nClick on ESC to stop selecting points.\n")
  idx <-
    identify(log2(1+df$baseMean),df$log2FoldChange, labels = df[,geneID])
  write.table(df[idx,],paste0(outputFolder,"/ValuesForGenesClickedInMAP.txt"),sep="\t",quote=F,row.names=F)
  cat("The values for the points identified are in",paste0(outputFolder,"/ValuesForGenesClickedInMAP.txt.\n"))
  if(usePng){
    dev.copy(png, paste0(outputFolder,"/MAP_clicked.png"))
  } else {
    dev.copy(pdf, paste0(outputFolder,"/MAP_clicked.pdf"),title="MAP_clicked")
  }
  dev.off()
  library(ggplot2)
  library(ggrepel)

  cmd<-paste("ggplot(data = df, aes(x = log2(1+baseMean), y = log2FoldChange)) + theme_classic() + 
  geom_point(colour = colorDots, size = 1) +
  ylab(\"log2 Fold Change\") +
  xlab(\"log2(1+Mean counts)\") +
  ylim(y=yLimUsed[1],yLimUsed[2]) +
  geom_text_repel(data = df[idx,] ,aes(label =",geneID,"), 
       box.padding = unit(0.45, \"lines\"))")
  if(exists("colOfCircle")){
    cmd<-paste0(cmd,"+
    geom_point(data=df[idx,],
             aes(log2(1+baseMean),log2FoldChange),shape=1,col=colOfCircle)")
  }
  if(usePng){
    png(paste0(outputFolder,"/MAP_clicked_pretty.png"))
  } else {
    pdf(paste0(outputFolder,"/MAP_clicked_pretty.pdf"),title="MAP_clicked_pretty")
  }
  print(eval(parse(text = cmd)))
  dev.off()
}

if(exists("dfGene")){
  library(ggplot2)
  library(ggrepel)
  cmd<-paste("ggplot(data = df, aes(x = log2(1+baseMean), y = log2FoldChange)) + theme_classic() + 
  geom_point(colour = colorDots, size = 1) +
  ylab(\"log2 Fold Change\") +
  xlab(\"log2(1+Mean counts)\") +
  ylim(y=yLimUsed[1],yLimUsed[2]) +
  geom_text_repel(data = df[df[,colOfGeneID]%in%dfGene[,1],] ,aes(label =",colOfGeneID,"), 
       box.padding = unit(0.45, \"lines\"))")
  if(exists("colOfCircle")){
    cmd<-paste0(cmd,"+
    geom_point(data=df[df[,colOfGeneID]%in%dfGene[,1],],
             aes(log2(1+baseMean), log2FoldChange),shape=1,col=colOfCircle)")
  }
  if(usePng){
    png(paste0(outputFolder,"/MAP_listOfGenes_pretty.png"))
  } else {
    pdf(paste0(outputFolder,"/MAP_listOfGenes_pretty.pdf"),title="MAP_listOfGenes_pretty")
  }
  print(eval(parse(text = cmd)))
  dev.off()
}
