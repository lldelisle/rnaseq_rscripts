options(stringsAsFactors=F)
rm(list=ls())

library(tools)

if(!("ggplot2" %in% rownames(installed.packages()))) {
  install.packages("ggplot2")
}

if (!("pheatmap" %in% rownames(installed.packages()))) {
  source("http://bioconductor.org/biocLite.R") #adds bioconductor site as a package source
  biocLite("pheatmap") #downloads and install pheatmap package from bioconductor
}
if (!("RColorBrewer" %in% rownames(installed.packages()))) {
  source("http://bioconductor.org/biocLite.R") #adds bioconductor site as a package source
  biocLite("RColorBrewer") #downloads and install pheatmap package from bioconductor
}


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
if(!exists("samplesPlan")){
  stop("The config file do not have samplesPlan definition.")
}
if(!file.exists(samplesPlan)){
  stop("The file specified as samplesPlan does not exist:",samplesPlan)
}
samplesPlanDF<-read.delim(samplesPlan)
if(!("sample"%in%colnames(samplesPlanDF))){
  stop("The samplesPlan table do not contain a column called \"sample\".")
}
setwd(dirname(samplesPlan))

if(!exists("tableWithNormalizedExpression")){
  stop("The config file do not have tableWithNormalizedExpression definition.")
}
if(!file.exists(tableWithNormalizedExpression)){
  stop("The file specified as tableWithNormalizedExpression:",tableWithNormalizedExpression,"does not exists.")
}

expressionDF<-read.delim(tableWithNormalizedExpression)

#Because the samples plan may contain information about samples that are not in the data we restrict the samples plan
samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)

if(exists("useFPKM")){
  if(is.logical(useFPKM)){
    if(useFPKM){
      expressionDF[,samplesToPlot]<-list(NULL)
      samplesToPlot<-unlist(lapply(strsplit(colnames(expressionDF),"^FPKM_"),function(x){if(length(x)>1){return(x[2])}}))
      colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
    }
  } else {
    useFPKM<-F
  }
} else {
  useFPKM<-F
}


if (length(samplesToPlot) < 1) {
  if(!useFPKM){
    samplesToPlot<-unlist(lapply(strsplit(colnames(expressionDF),"^FPKM_"),function(x){if(length(x)>1){return(x[2])}}))
    if(length(samplesToPlot)<1){
      stop("The samplesPlan table is incompatible with the table with expression values.")
    } else {
      useFPKM<-T
      cat("FPKM values will be used.\n")
      colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
    }
  } else {
    stop("The samplesPlan table is incompatible with the table with expression values.")
  }
}

if(!exists("RNAseqFunctionPath")){
    stop("The RNAseqFunctionPath is not provided.")
} else {
  if(!file.exists(RNAseqFunctionPath)){
    stop("The file provided in RNAseqFunctionPath:",RNAseqFunctionPath," does not exists.")
  }
}
source(RNAseqFunctionPath)
factorizedSP<-simplifyDF(samplesPlanDF,samplesToPlot,keepFactorsWithOneValue=T)

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

#The data are restricted to the samples to plot and to non 0 values
data <- expressionDF[,samplesToPlot]
sumperline <- apply(data,1,sum)
cat("Only genes with at least non-null expression in one sample are considered.\n")
nonZdata <- data[sumperline != 0,]
    
#The data are transformed into log data
cat("Data are transformed to log2(1+expression).\n")
ldata <- log2(nonZdata + 1)

rldata <- ldata
if(exists("restrictToNMoreVariantGenes")){
  if(is.numeric(restrictToNMoreVariantGenes)){
    rldata <- ldata[order(apply(ldata,1,var),decreasing = T)[1:min(nrow(ldata),restrictToNMoreVariantGenes)],]
    cat("Only the most variant expression values are used.\n")
  } else {
    cat("The value put in restrictToNMoreVariantGenes is not numeric. The PCA will not be restricted to a subset of genes.\n")
  }
}

if(exists("nbOfPC")){
  if(is.numeric(nbOfPC)){
    if(nbOfPC>0){
      cat("Performing PCA...\n")      
      sample.pca <- prcomp(t(rldata),
                         center = TRUE,
                         scale. = FALSE)
      cat("The PCA is performed on center but not scaled data.\n")
      #the data with the PCA coordinates and the factors are grouped
      new.df <- data.frame(factorizedSP,sample.pca$x[samplesToPlot,])
      #var contains the variance for each PC
      var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)
      library("ggplot2")
      #####
      #Plot one PC per one PC:
      if(exists("PCA1D")){
        param<-paste0(getStringFromListAndSP(PCA1D,colnames(factorizedSP),possibleValues=c("fill","alpha","color","linetype")),")")
        if(!grepl("color",param) && grepl("linetype",param)){
          param<-paste0(param,", color=\"black\"")
        }
      } else {
        param<-")"
      }
      for(i in 1:nbOfPC){
        if(usePng){
          png(paste0(outputFolder,"/PC",i,".png"))
        } else {
          pdf(paste0(outputFolder,"/PC",i,".pdf"),title=paste0("PC",i))
        }
        cmd<-paste0("ggplot(new.df, aes(sample)) +
                      geom_bar(aes(weight=PC",i,param,",size=1.5) +
                      theme_grey(base_size = 20) +
                      theme(axis.text.x  = element_text(angle=90,vjust=0.5,hjust=1)) +
                      ylab(paste0(\"PC",i,": \",var[",i,"],\"% variance\"))")
        print(eval(parse(text = cmd)))
        dev.off()
      }
      if(nbOfPC>2){
        #####
        #Plot PCs 2 per 2:
        if(exists("PCA2D")){
          param<-getStringFromListAndSP(PCA2D,colnames(factorizedSP),possibleValues=c("fill","alpha","color","shape"))
          if(nchar(param)>0){
            param<-paste0("aes(",param,")")
            if(!grepl("shape",param) && grepl("fill",param)){
              param<-paste0(param,", shape=21, stroke=2")
            } else if(grepl("fill",param)){
              param<-paste0(param,", stroke=2")
              additionalParam<-paste0("+ scale_shape_manual(values=21:",20+length(levels(factorizedSP[,PCA2D$shape])),") +
              guides(fill=guide_legend(override.aes = list(shape = 21)),alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))")
            }
          }
        } else {
          param<-""
        }
        for(i in 1:(nbOfPC-1)){
          for(j in (i+1):nbOfPC){
            if(usePng){
              png(paste0(outputFolder,"/PC",i,"-PC",j,".png"))
            } else {
              pdf(paste0(outputFolder,"/PC",i,"-PC",j,".pdf"),title=paste0("PC",i,"-PC",j))
            }
            cmd<-paste0("ggplot(new.df, aes(PC",i,",PC",j,")) +
                          geom_point(",param,",size=3) +
                          theme_grey(base_size = 20) +
                          xlab(paste0(\"PC",i,": \",var[",i,"],\"% variance\"))+
                          ylab(paste0(\"PC",j,": \",var[",j,"],\"% variance\"))")
            if(exists("additionalParam")){
              cmd<-paste0(cmd,additionalParam)
            }
            print(eval(parse(text = cmd)))
            dev.off()
          }
        }
      }
      if(exists("getGeneContributionToPCA")){
        if(is.logical(getGeneContributionToPCA)){
          if(getGeneContributionToPCA){
            pcaDF<-cbind(expressionDF[rownames(sample.pca$rotation),setdiff(colnames(expressionDF),samplesToPlot)],sample.pca$rotation*sample.pca$rotation)
            write.table(pcaDF,paste0(outputFolder,"/geneContributionToPCA.txt"),sep = "\t",row.names = F,quote=F)
            cat("The table with gene contribution is written in ",outputFolder,"/geneContributionToPCA.txt.\n")
          }
        } else {
          cat("The value of getGeneContributionToPCA is not logical. No table will be written.\n")
        }
      }
    }
  } else {
    cat("nbOfPC is not numeric. No PCA will be performed.\n")
  }
}

if(exists("plotMatrixAndClustering")){
  if(is.logical(plotMatrixAndClustering)){
    if(plotMatrixAndClustering){
      cat("Performing correlation and clustering...\n")
      cat("The distances are euclidean distances.\n")
      cat("The clustering method is complete.\n")
      sampleDists <- dist(t(rldata))
      library("pheatmap")
      library("RColorBrewer")
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- samplesToPlot
      colnames(sampleDistMatrix) <- samplesToPlot
      colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      #an annot data frame is made only with the factors used in the pca
      cols<-NULL
      if(exists("PCA1D")){
        cols <- unique(c(cols,which(colnames(factorizedSP) %in% PCA1D)))
      }
      if(exists("PCA2D")){
        cols<-unique(c(cols,which(colnames(factorizedSP) %in% PCA2D)))
      }
      annot <- NULL
      if (length(cols) == 1) {
        fi <- cols
        annot <- data.frame(factorizedSP[,fi])
        colnames(annot) <- colnames(factorizedSP)[fi]
        rownames(annot) <- samplesToPlot
      } else if(length(cols)>1) {
        annot<-factorizedSP[,cols]
      }
      if(usePng){
        png(paste0(outputFolder,"/CorrelationMatrix.png"))
      } else {
        pdf(paste0(outputFolder,"/CorrelationMatrix.pdf"),title="CorrelationMatrix",onefile = FALSE)
      }  
      pheatmap(
        sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        cellwidth = 10,
        cellheight = 10,
        annotation = annot,
        col = colors
      )
      dev.off()
    }
  } else {
    cat("plotMatrixAndClustering is not logical. No clustering was performed.\n")
  }
}

if(exists("fileWithGenes")){
  if(file.exists(fileWithGenes)){
    dfGene<-read.delim(fileWithGenes)
    #colOfGeneID indicates in which column the gene list should be look for
    colOfGeneID <- colnames(dfGene)[1]
    if (!(colOfGeneID %in% colnames(expressionDF))) {
      stop("The first line of the gene file does not correspond to a column in the expression file.")
    }
    if(exists("xaxisForGenes")){
      if(!xaxisForGenes%in%colnames(factorizedSP)){
        cat("The factor provided as xaxisForGenes is not part of the samplesPlan table. The sample name will be used.\n")
        xaxisForGenes<-"sample"
      }
    } else {
      cat("The xaxisForGenes is not specified in the config file. The sample name will be used.\n")
      xaxisForGenes<-"sample"
    }
    ylab<-"Normalized counts"
    if(useFPKM){
      ylab<-"FPKM"
    }
    if(exists("useLogExpression")){
      if(is.logical(useLogExpression)){
        if(useLogExpression){
          ylab<-paste0("log2(1+",ylab,")")
          dataToPlot<-log2(1+data)
        } else {
          dataToPlot<-data
        }
      } else {
        cat("useLogExpression is not logicial. Will use log values.\n")
          ylab<-paste0("log2(1+",ylab,")")
          dataToPlot<-log2(1+data)
      }
    } else {
      cat("Will use log values.\n")
      ylab<-paste0("log2(1+",ylab,")")
      dataToPlot<-log2(1+data)
    }
    #Set the parameters:
    additionalParam<-""
    if(exists("plotGenesPara")){
      param<-getStringFromListAndSP(plotGenesPara,colnames(factorizedSP),possibleValues=c("fill","alpha","color","shape"))
      if(nchar(param)>0){
        param<-paste0("aes(",param,")")
        if(!grepl("shape",param) && grepl("fill",param)){
          param<-paste0(param,", shape=21, stroke=2")
        } else if(grepl("fill",param)){
          param<-paste0(param,", stroke=2")
          additionalParam<-paste0("+ scale_shape_manual(values=21:",20+length(levels(factorizedSP[,PCA2D$shape])),") +
          guides(fill=guide_legend(override.aes = list(shape = 21)),alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))")
        }
      }
    } else {
      param<-""
    }    
    #Set all additional parameters:
    #y lim
    if(exists("useSameYmaxForAllGenes")){
      if(is.logical(useSameYmaxForAllGenes)){
        if(useSameYmaxForAllGenes){
          additionalParam<-paste0(additionalParam,"+ expand_limits(y=c(0,",max(dataToPlot[expressionDF[,colOfGeneID]%in%dfGene[,1],samplesToPlot]),"))")
        } else {
          additionalParam<-paste0(additionalParam,"+ expand_limits(y=0)")
        }
      } else {
        cat("useSameYmaxForAllGenes is not logical. Will be set to F.\n")
        additionalParam<-paste0(additionalParam,"+ expand_limits(y=0)")
      }
    } else {
      additionalParam<-paste0(additionalParam,"+ expand_limits(y=0)")
    }
    #x label if too numerous:
    if(length(levels(factorizedSP[,xaxisForGenes]))>4){
      additionalParam<-paste0(additionalParam,"+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))")
    }
    library("ggplot2")
    for (i in 1:nrow(dfGene)) {
      # for each gene in the fg file the ggplot command line is applied and the plot is stored in a file in the directory of the gene list.
      geneID <- dfGene[i,1]
      if (geneID %in% expressionDF[,colOfGeneID]) {
        cat("plotting",geneID,"\n")
        subdf <- dataToPlot[expressionDF[,colOfGeneID] == geneID,]
        new.df <- cbind(factorizedSP,t(subdf[1,]))
        colnames(new.df)[ncol(new.df)]<-"y"
        if(usePng){
          png(paste0(outputFolder,"/",xaxisForGenes,"-",geneID,".png"))
        } else {
          pdf(paste0(outputFolder,"/",xaxisForGenes,"-",geneID,".pdf"),title=paste0(xaxisForGenes,"-",geneID))
        }
        cmd<-paste0("ggplot(new.df, aes(",xaxisForGenes,",y)) +
                      geom_point(",param,",size=3)  +
                      labs(title = geneID) +
                      theme_grey(base_size = 20) +
                      ylab(\"",ylab,"\")",additionalParam)
        print(eval(parse(text = cmd)))
        dev.off()
      } else{
        cat(paste(geneID,"is not in the file with expression values.\n"))
      }
    }
  } else {
    cat("The file provided as fileWithGenes does not exists.\n")
  }
}
