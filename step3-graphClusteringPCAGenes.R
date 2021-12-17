options(stringsAsFactors=F)
rm(list=ls())

library(tools)

if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")

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
fileWithAllCommands<-paste0(dirname(f),"/",basename(file_path_sans_ext(f)),"_commandLinesLaunched.R")
cat("options(stringsAsFactors=F)\n",file=fileWithAllCommands)
#cat(paste0("f <- \"",f,"\"\n"),file=fileWithAllCommands,append=T)
source(f)
#cat("source(f)\n",file=fileWithAllCommands,append=T)
if(!exists("samplesPlan")){
  stop("The config file do not have samplesPlan definition.")
}
if(!file.exists(samplesPlan)){
  stop("The file specified as samplesPlan does not exist:",samplesPlan)
}
cat(paste0("samplesPlan <- \"",samplesPlan,"\"\n"),file=fileWithAllCommands,append=T)
samplesPlanDF<-read.delim(samplesPlan)
cat("samplesPlanDF<-read.delim(samplesPlan)\n",file=fileWithAllCommands,append=T)
if(!("sample"%in%colnames(samplesPlanDF))){
  stop("The samplesPlan table do not contain a column called \"sample\".")
}

if(!exists("tableWithNormalizedExpression")){
  stop("The config file do not have tableWithNormalizedExpression definition.")
}
if(!file.exists(tableWithNormalizedExpression)){
  stop("The file specified as tableWithNormalizedExpression:",tableWithNormalizedExpression,"does not exists.")
}
cat(paste0("tableWithNormalizedExpression <- \"",tableWithNormalizedExpression,"\"\n"),file=fileWithAllCommands,append=T)

expressionDF<-read.delim(tableWithNormalizedExpression)
cat("expressionDF<-read.delim(tableWithNormalizedExpression)\n",file=fileWithAllCommands,append=T)
metaCols<-which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn])!="numeric"}))
cat("metaCols<-which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn])!=\"numeric\"}))\n",file=fileWithAllCommands,append=T)
#Because the samples plan may contain information about samples that are not in the data we restrict the samples plan
samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)
cat("samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)\n",file=fileWithAllCommands,append=T)

if(exists("useFPKM")){
  if(is.logical(useFPKM)){
    if(useFPKM){
      expressionDF[,samplesToPlot]<-list(NULL)
      colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
      samplesToPlot<-intersect(colnames(expressionDF),samplesPlanDF$sample)
      cat("expressionDF[,samplesToPlot]<-list(NULL)
colnames(expressionDF)<-gsub(\"^FPKM_\",\"\",colnames(expressionDF))
samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)\n",file=fileWithAllCommands,append=T)
    }
  } else {
    useFPKM<-F
  }
} else {
  useFPKM<-F
}


if (length(samplesToPlot) < 1) {
  if(!useFPKM){
    colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
    samplesToPlot<-intersect(colnames(expressionDF),samplesPlanDF$sample)
    if(length(samplesToPlot)<1){
      stop("The samplesPlan table is incompatible with the table with expression values.")
    } else {
      useFPKM<-T
      cat("FPKM values will be used.\n")
      cat("colnames(expressionDF)<-gsub(\"^FPKM_\",\"\",colnames(expressionDF))
samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)\n",file=fileWithAllCommands,append=T)
    }
  } else {
    stop("The samplesPlan table is incompatible with the table with expression values.")
  }
}

nSamples<-length(samplesToPlot)

if(!exists("RNAseqFunctionPath")){
  stop("The RNAseqFunctionPath is not provided.")
} else {
  if(!file.exists(RNAseqFunctionPath)){
    stop("The file provided in RNAseqFunctionPath:",RNAseqFunctionPath," does not exists.")
  }
}
cat(paste0("RNAseqFunctionPath <- \"",RNAseqFunctionPath,"\"\n"),file=fileWithAllCommands,append=T)
source(RNAseqFunctionPath)
cat("source(RNAseqFunctionPath)\n",file=fileWithAllCommands,append=T)
factorizedSP<-simplifyDF(samplesPlanDF,samplesToPlot,keepFactorsWithOneValue=T)
cat("factorizedSP<-simplifyDF(samplesPlanDF,samplesToPlot,keepFactorsWithOneValue=T)\n",file=fileWithAllCommands,append=T)



if(!exists("outputFolder")){
  outputFolder<-paste0(dirname(samplesPlan),"/plots")
  cat("outputFolder<-paste0(dirname(samplesPlan),\"/plots\")\n", file=fileWithAllCommands,append=T)
} else {
  cat(paste0("outputFolder <- \"",outputFolder,"\"\n"),file=fileWithAllCommands,append=T)
}
dir.create(outputFolder,recursive=T,showWarnings=F)
cat("dir.create(outputFolder,recursive=T,showWarnings=F)\n", file=fileWithAllCommands,append=T)
if(exists("usePng")){
  if(!is.logical(usePng)){
    cat("usePng is not logical. pdf will be used as output.\n")
    usePng<-F
  }
} else {
  usePng<-F
}

if(exists("fixedColors")){
  fixedColors<-checkedFixedColors(fixedColors,factorizedSP)
  cat("fixedColors<-",stringFromListOfNamedVec(fixedColors),"\n", file=fileWithAllCommands,append=T)
} else {
  fixedColors<-NULL
  cat("fixedColors<-NULL\n", file=fileWithAllCommands,append=T)
}

#The data are restricted to the samples to plot 
data <- expressionDF[,samplesToPlot]
#and to non 0 values for PCA and clustering
sumperline <- apply(data,1,sum)
cat("Only genes with at least non-null expression in one sample are considered.\n")
nonZdata <- data[sumperline != 0,]
#The data are transformed into log data
cat("Data are transformed to log2(1+expression).\n")
ldata <- log2(nonZdata + 1)

cat("data <- expressionDF[,samplesToPlot]
sumperline <- apply(data,1,sum)
nonZdata <- data[sumperline != 0,]
ldata <- log2(nonZdata + 1)\n", file=fileWithAllCommands,append=T)
rldata <- ldata
cat("rldata <- ldata\n", file=fileWithAllCommands,append=T)
if(exists("restrictToNMoreVariantGenes")){
  if(is.numeric(restrictToNMoreVariantGenes)){
    rldata <- ldata[order(apply(ldata,1,var),decreasing = T)[1:min(nrow(ldata),restrictToNMoreVariantGenes)],]
    cat("rldata <- ldata[order(apply(ldata,1,var),decreasing = T)[1:min(nrow(ldata),",restrictToNMoreVariantGenes,")],]\n", file=fileWithAllCommands,append=T)
    cat("Only the most variant expression values are used.\n")
  } else {
    cat("The value put in restrictToNMoreVariantGenes is not numeric. The PCA will not be restricted to a subset of genes.\n")
  }
}

if(exists("nbOfPC")){
  if(is.numeric(nbOfPC)){
    cat(paste0("nbOfPC <- ",nbOfPC,"\n"),file=fileWithAllCommands,append=T)
    if(nbOfPC>0){
      cat("Performing PCA...\n")      
      sample.pca <- prcomp(t(rldata),
                           center = TRUE,
                           scale. = FALSE)
      cat("sample.pca <- prcomp(t(rldata),
                         center = TRUE,
                         scale. = FALSE)\n", file=fileWithAllCommands,append=T)
      cat("The PCA is performed on center but not scaled data.\n")
      #the data with the PCA coordinates and the factors are grouped
      new.df <- data.frame(factorizedSP,sample.pca$x[samplesToPlot,])
      cat("new.df <- data.frame(factorizedSP,sample.pca$x[samplesToPlot,])\n", file=fileWithAllCommands,append=T)
      #var contains the variance for each PC
      var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)
      cat("var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)\n", file=fileWithAllCommands,append=T)
      cat("library(ggplot2)\n", file=fileWithAllCommands,append=T)
      #####
      #Plot one PC per one PC:
      if(exists("PCA1D")){
        param<-paste0(getStringFromListAndSP(PCA1D,colnames(factorizedSP),possibleValues=c("fill","alpha","color","linetype")),")")
        addParam<-gsub("proposedColors","fixedColors",getAddPara(PCA1D,fixedColors))
        if(!grepl("color",param) && grepl("linetype",param)){
          param<-paste0(param,", color=\"black\"")
        }
      } else {
        param<-")"
        addParam<-""
      }
      pdfSize<-max(7,6+0.12*nSamples)
      pngSize<-max(500,358+9*nSamples)
      for(i in 1:nbOfPC){
        if(usePng){
          png(paste0(outputFolder,"/PC",i,".png"),width=pngSize,height=pngSize)
        } else {
          pdf(paste0(outputFolder,"/PC",i,".pdf"),title=paste0("PC",i),width=pdfSize,height=pdfSize)
        }
        cmd<-paste0("ggplot(new.df, aes(sample)) +
                      geom_bar(aes(weight=PC",i,param,",size=1.5) +
                      theme_grey(base_size = 20) +
                      theme(axis.text.x  = element_text(angle=90,vjust=0.5,hjust=1)) +
                      ylab(paste0(\"PC",i,": \",var[",i,"],\"% variance\"))",addParam)
        print(eval(parse(text = cmd)))
        dev.off()
      }
      cat("for(i in 1:nbOfPC){\n", file=fileWithAllCommands,append=T)
      if(usePng){
        cat("png(paste0(outputFolder,\"/PC\",i,\".png\"),width=",pngSize,",height=",pngSize,")\n", file=fileWithAllCommands,append=T)
      } else {
        cat("pdf(paste0(outputFolder,\"/PC\",i,\".pdf\"),title=paste0(\"PC\",i),width=",pdfSize,",height=",pdfSize,")\n", file=fileWithAllCommands,append=T)
      }
      cat("cmd<-paste0(\"ggplot(new.df, aes(sample)) +
                      geom_bar(aes(weight=PC\",i,\"",gsub("\\\"","\\\\\"",param),",size=1.5) +
                      theme_grey(base_size = 20) +
                      theme(axis.text.x  = element_text(angle=90,vjust=0.5,hjust=1)) +
                      ylab(paste0(\\\"PC\",i,\": \\\",var[\",i,\"],\\\"% variance\\\"))",gsub("\\\"","\\\\\"",addParam),"\")\n", file=fileWithAllCommands,append=T)
      cat("print(eval(parse(text = cmd)))\n", file=fileWithAllCommands,append=T)
      cat("dev.off()\n", file=fileWithAllCommands,append=T)                      
      cat("}\n", file=fileWithAllCommands,append=T)
      if(nbOfPC>1){
        #####
        #Plot PCs 2 per 2:
        if(exists("PCA2D")){
          PCA2D<-checkedPCA2D(PCA2D,factorizedSP)
          param<-getStringFromListAndSP(PCA2D,colnames(factorizedSP),possibleValues=c("fill","alpha","color","shape"))
          addParam<-gsub("proposedColors","fixedColors",getAddPara(PCA2D,fixedColors))
          if(nchar(param)>0){
            param<-paste0("aes(",param,")")
            if(!grepl("shape",param) && grepl("fill",param)){
              param<-paste0(param,", shape=21, stroke=2")
            } else if(grepl("fill",param)){
              param<-paste0(param,", stroke=2")
              addParam<-paste0(addParam,"+\n scale_shape_manual(values=c(",paste(rep(21:25,length.out=length(levels(factorizedSP[,PCA2D$shape]))),collapse = ","),")) +
guides(fill=guide_legend(override.aes = list(shape = 21)),alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))")
            }
          }
        } else {
          param<-""
          addParam<-""
        }
        pdfSize<-max(7,6.7+0.035*nSamples)
        pngSize<-max(480,460+2.4*nSamples)
        for(i in 1:(nbOfPC-1)){
          for(j in (i+1):nbOfPC){
            if(usePng){
              png(paste0(outputFolder,"/PC",i,"-PC",j,".png"),width=pngSize,height=pngSize)
            } else {
              pdf(paste0(outputFolder,"/PC",i,"-PC",j,".pdf"),title=paste0("PC",i,"-PC",j),width=pdfSize,height=pdfSize)
            }
            cmd<-paste0("ggplot(new.df, aes(PC",i,",PC",j,")) +
                          geom_point(",param,",size=3) +
                          theme_grey(base_size = 20) +
                          xlab(paste0(\"PC",i,": \",var[",i,"],\"% variance\"))+
                          ylab(paste0(\"PC",j,": \",var[",j,"],\"% variance\"))",addParam)
            print(eval(parse(text = cmd)))
            dev.off()
          }
        }
        cat("for(i in 1:(nbOfPC-1)){\n", file=fileWithAllCommands,append=T)
        cat("for(j in (i+1):nbOfPC){\n", file=fileWithAllCommands,append=T)
        if(usePng){
          cat("png(paste0(outputFolder,\"/PC\",i,\"-PC\",j,\".png\"),width=",pngSize,",height=",pngSize,")\n", file=fileWithAllCommands,append=T)
        } else {
          cat("pdf(paste0(outputFolder,\"/PC\",i,\"-PC\",j,\".pdf\"),title=paste0(\"PC\",i,\"-PC\",j),width=",pdfSize,",height=",pdfSize,")\n", file=fileWithAllCommands,append=T)
        }
        cat("cmd<-paste0(\"ggplot(new.df, aes(PC\",i,\",PC\",j,\")) +
geom_point(",gsub("\\\"","\\\\\"",param),",size=3) +
theme_grey(base_size = 20) +
xlab(paste0(\\\"PC\",i,\": \\\",var[\",i,\"],\\\"% variance\\\"))+
ylab(paste0(\\\"PC\",j,\": \\\",var[\",j,\"],\\\"% variance\\\"))",gsub("\\\"","\\\\\"",addParam),"\")\n", file=fileWithAllCommands,append=T)
        cat("print(eval(parse(text = cmd)))\n", file=fileWithAllCommands,append=T)
        cat("dev.off()\n", file=fileWithAllCommands,append=T)                      
        cat("}\n", file=fileWithAllCommands,append=T)                    
        cat("}\n", file=fileWithAllCommands,append=T)
      }
      if(exists("getGeneContributionToPCA")){
        if(is.logical(getGeneContributionToPCA)){
          if(getGeneContributionToPCA){
            pcaDF<-cbind(expressionDF[rownames(sample.pca$rotation),metaCols],sample.pca$rotation*sample.pca$rotation)
            write.table(pcaDF,paste0(outputFolder,"/geneContributionToPCA.txt"),sep = "\t",row.names = F,quote=F)
            cat("The table with gene contribution is written in ",outputFolder,"/geneContributionToPCA.txt.\n",sep = "")
            
            cat("pcaDF<-cbind(expressionDF[rownames(sample.pca$rotation),metaCols],sample.pca$rotation*sample.pca$rotation)
write.table(pcaDF,paste0(outputFolder,\"/geneContributionToPCA.txt\"),sep = \"\\t\",row.names = F,quote=F)\n", file=fileWithAllCommands,append=T)
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
      sampleDists <- dist(t(rldata))
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- samplesToPlot
      colnames(sampleDistMatrix) <- samplesToPlot
      colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      cat("sampleDists <- dist(t(rldata),method=\"euclidean\")
library(pheatmap)
library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samplesToPlot
colnames(sampleDistMatrix) <- samplesToPlot
colors <- colorRampPalette(rev(brewer.pal(9, \"Blues\")))(255)\n", file=fileWithAllCommands,append=T)
      #an annot data frame is made only with the factors used in the pca
      cols<-NULL
      if(exists("PCA1D")){
        cols <- unique(c(cols,which(colnames(factorizedSP) %in% PCA1D)))
      }
      if(exists("PCA2D")){
        cols<-unique(c(cols,which(colnames(factorizedSP) %in% PCA2D)))
      }
      if(!is.null(cols)){
        cat("cols<-c(",paste0(cols,collapse=","),")\n", file=fileWithAllCommands,append=T)
      } else {
        cat("cols<-NULL\n", file=fileWithAllCommands,append=T)
      }
      annot <- NULL
      cat("annot <- NULL\n", file=fileWithAllCommands,append=T)
      if (length(cols) == 1) {
        annot <- data.frame(factorizedSP[,cols])
        colnames(annot) <- colnames(factorizedSP)[cols]
        rownames(annot) <- samplesToPlot
        cat("annot <- data.frame(factorizedSP[,cols])
colnames(annot) <- colnames(factorizedSP)[cols]
rownames(annot) <- samplesToPlot\n", file=fileWithAllCommands,append=T)
      } else if(length(cols)>1) {
        annot<-factorizedSP[,cols]
        cat("annot<-factorizedSP[,cols]\n", file=fileWithAllCommands,append=T)
      }
      pdfSize<-max(7,5+0.15*nSamples)
      pngSize<-max(480,237+12*nSamples)
      if(usePng){
        png(paste0(outputFolder,"/CorrelationMatrix_EuclComp.png"),width=pngSize,height=pngSize)
        cat("png(paste0(outputFolder,\"/CorrelationMatrix_EuclComp.png\"),width=",pngSize,",height=",pngSize,")\n", file=fileWithAllCommands,append=T)
      } else {
        pdf(paste0(outputFolder,"/CorrelationMatrix.pdf"),title="CorrelationMatrix",onefile = TRUE,width=pdfSize,height=pdfSize)
        cat("pdf(paste0(outputFolder,\"/CorrelationMatrix.pdf\"),title=\"CorrelationMatrix\",onefile = TRUE,width=",pdfSize,",height=",pdfSize,")\n", file=fileWithAllCommands,append=T)
      }  
      pheatmap(
        sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        cellwidth = 10,
        cellheight = 10,
        annotation = annot,
        annotation_colors = fixedColors,
        main="Euclidean distance - complete clustering",
        clustering_method="complete",
        col = colors
      )
      cat("pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  cellwidth = 10,
  cellheight = 10,
  annotation = annot,
  annotation_colors = fixedColors,
  main=\"Euclidean distance - complete clustering\",
  col = colors,
  clustering_method=\"complete\"
)\n", file=fileWithAllCommands,append=T)
      if(usePng){
        dev.off()
        png(paste0(outputFolder,"/CorrelationMatrix_EuclWard.png"),width=pngSize,height=pngSize)
        cat("dev.off()
png(paste0(outputFolder,\"/CorrelationMatrix_EuclWard.png\"),width=",pngSize,",height=",pngSize,")\n", file=fileWithAllCommands,append=T)
      }
      pheatmap(
        sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        cellwidth = 10,
        cellheight = 10,
        annotation = annot,
        annotation_colors = fixedColors,
        main="Euclidean distance - ward clustering",
        clustering_method="ward.D2",
        col = colors
      )
      cat("pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  cellwidth = 10,
  cellheight = 10,
  annotation = annot,
  annotation_colors = fixedColors,
  main=\"Euclidean distance - ward clustering\",
  clustering_method=\"ward.D2\",
  col = colors
)\n", file=fileWithAllCommands,append=T)
      if(usePng){
        dev.off()
        png(paste0(outputFolder,"/CorrelationMatrix_SpearWard.png"),width=pngSize,height=pngSize)
        cat("dev.off()
png(paste0(outputFolder,\"/CorrelationMatrix_SpearWard.png\"),width=",pngSize,",height=",pngSize,")\n", file=fileWithAllCommands,append=T)
      }
      correlationMatrix <- cor(rldata, method="spearman")
      newSampleDist <- as.dist(1 - correlationMatrix)
      pheatmap(
        correlationMatrix,
        clustering_distance_rows = newSampleDist,
        clustering_distance_cols = newSampleDist,
        cellwidth = 10,
        cellheight = 10,
        annotation = annot,
        annotation_colors = fixedColors,
        main="spearmanCor - ward clustering",
        clustering_method="ward.D2",
        col = rev(colors)
      )
      dev.off()
      cat("correlationMatrix <- cor(rldata, method=\"spearman\")
newSampleDist <- as.dist(1 - correlationMatrix)
pheatmap(
    correlationMatrix,
    clustering_distance_rows = newSampleDist,
    clustering_distance_cols = newSampleDist,
    cellwidth = 10,
    cellheight = 10,
    annotation = annot,
    annotation_colors = fixedColors,
    main=\"spearmanCor - ward clustering\",
    clustering_method=\"ward.D2\",
    col = rev(colors)
)
dev.off()\n", file=fileWithAllCommands,append=T)
    }
  } else {
    cat("plotMatrixAndClustering is not logical. No clustering was performed.\n")
  }
}

if(exists("fileWithGenes")){
  if(file.exists(fileWithGenes)){
    cat(paste0("fileWithGenes <- \"",fileWithGenes,"\"\n"),file=fileWithAllCommands,append=T)
    dfGene<-read.delim(fileWithGenes)
    cat("dfGene<-read.delim(fileWithGenes)\n", file=fileWithAllCommands,append=T)
    #colOfGeneID indicates in which column the gene list should be look for
    colOfGeneID <- colnames(dfGene)[1]
    cat("colOfGeneID <- colnames(dfGene)[1]\n", file=fileWithAllCommands,append=T)
    if (!(colOfGeneID %in% colnames(expressionDF))) {
      stop("The first line of the gene file does not correspond to a column in the expression file.")
    }
    if(exists("geneIDToAdd")){
      if(geneIDToAdd%in%colnames(expressionDF)){
        cat(paste0("geneIDToAdd <- \"",geneIDToAdd,"\"\n"),file=fileWithAllCommands,append=T)
      } else {
        cat("The geneIDToAdd provided is not part of the columns of tableWithNormalizedExpression. It will be omitted.\n")
        rm(geneIDToAdd)
      }
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
          cat("dataToPlot<-log2(1+data)\n", file=fileWithAllCommands,append=T)
        } else {
          dataToPlot<-data
          cat("dataToPlot<-data\n", file=fileWithAllCommands,append=T)
        }
      } else {
        cat("useLogExpression is not logicial. Will use log values.\n")
        useLogExpression<-T
        ylab<-paste0("log2(1+",ylab,")")
        dataToPlot<-log2(1+data)
        cat("dataToPlot<-log2(1+data)\n", file=fileWithAllCommands,append=T)
      }
    } else {
      cat("Will use log values.\n")
      useLogExpression<-T
      ylab<-paste0("log2(1+",ylab,")")
      dataToPlot<-log2(1+data)
      cat("dataToPlot<-log2(1+data)\n", file=fileWithAllCommands,append=T)
    }
    if(exists("doNotPlotGeneByGene")){
      if(!is.logical(doNotPlotGeneByGene)){
        cat("doNotPlotGeneByGene is not logical. Each gene will be plotted.\n")
        doNotPlotGeneByGene<-F
      }
    } else {
      doNotPlotGeneByGene<-F
    }
    if(!doNotPlotGeneByGene){
      #Set the parameters:
      if(exists("plotGenesPara")){
        param<-getStringFromListAndSP(plotGenesPara,colnames(factorizedSP),possibleValues=c("fill","alpha","color","shape"))
        additionalParam<-gsub("proposedColors","fixedColors",getAddPara(plotGenesPara,fixedColors))
        if(nchar(param)>0){
          param<-paste0("aes(",param,")")
          if(!grepl("shape",param) && grepl("fill",param)){
            param<-paste0(param,", shape=21, stroke=2")
          } else if(grepl("fill",param)){
            param<-paste0(param,", stroke=2")
            additionalParam<-paste0(additionalParam,"+ scale_shape_manual(values=c(",paste(rep(21:25,length.out=length(levels(factorizedSP[,plotGenesPara$shape]))),collapse = ","),")) +
          guides(fill=guide_legend(override.aes = list(shape = 21)),alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))")
          }
        }
      } else {
        param<-""
        additionalParam<-""
        plotGenesPara<-list()
      }    
      #Set all additional parameters:
      #y lim
      if(exists("useSameYmaxForAllGenes")){
        if(is.logical(useSameYmaxForAllGenes)){
          if(useSameYmaxForAllGenes){
            additionalParam<-paste0(additionalParam,"+ expand_limits(y=c(0,",max(dataToPlot[expressionDF[,colOfGeneID]%in%dfGene[,1],]),"))")
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
      nbXvalues<-length(levels(factorizedSP[,xaxisForGenes]))
      if(nbXvalues>4){
        additionalParam<-paste0(additionalParam,"+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))")
      }
      cat("library(ggplot2)\n", file=fileWithAllCommands,append=T)
      pdfSize<-max(7,6.5+0.12*nbXvalues)
      pngSize<-max(480,445+5*nbXvalues)
      cat("plotting ")
      for (i in 1:nrow(dfGene)) {
        # for each gene in the fg file the ggplot command line is applied and the plot is stored in a file in the directory of the gene list.
        geneID <- dfGene[i,1]
        if (geneID %in% expressionDF[,colOfGeneID]) {
          cat(geneID,",")
          subdf <- dataToPlot[expressionDF[,colOfGeneID] == geneID,]
          if(nrow(subdf)>1){
            cat("\nWarning: there is more than one line with the ID:",geneID,"in the file with expression values.\nThe first one will be plotted.\n")
            print(expressionDF[expressionDF[,colOfGeneID] == geneID,metaCols])
          }
          new.df <- cbind(factorizedSP,t(subdf[1,]))
          colnames(new.df)[ncol(new.df)]<-"y"
          if(exists("geneIDToAdd")){
            geneLabel<-expressionDF[expressionDF[,colOfGeneID]==geneID,geneIDToAdd][1]
            titleToPut<-paste(geneID,geneLabel,sep="-")
          } else{
            titleToPut<-geneID
          }
          if(usePng){
            png(paste0(outputFolder,"/",xaxisForGenes,"-",titleToPut,".png"),width=pngSize,height=pngSize)
          } else {
            pdf(paste0(outputFolder,"/",xaxisForGenes,"-",titleToPut,".pdf"),title=paste0(xaxisForGenes,"-",titleToPut),width=pdfSize,height=pdfSize)
          }
          cmd<-paste0("ggplot(new.df, aes(",xaxisForGenes,",y)) +
                      geom_point(",param,",size=3)  +
                      labs(title = titleToPut) +
                      theme_grey(base_size = 20) +
                      ylab(\"",ylab,"\")",additionalParam)
          print(eval(parse(text = cmd)))
          dev.off()
        } else{
          cat(paste(geneID,"is not in the file with expression values.\n"))
        }
      }
      cat("\n")
      cat("for (i in 1:nrow(dfGene)) {
geneID <- dfGene[i,1]
if (geneID %in% expressionDF[,colOfGeneID]) {
subdf <- dataToPlot[expressionDF[,colOfGeneID] == geneID,]
new.df <- cbind(factorizedSP,t(subdf[1,]))
colnames(new.df)[ncol(new.df)]<-\"y\"\n", file=fileWithAllCommands,append=T)
      if(exists("geneIDToAdd")){
        cat("geneLabel<-expressionDF[expressionDF[,colOfGeneID]==geneID,geneIDToAdd][1]
          titleToPut<-paste(geneID,geneLabel,sep=\"-\")\n", file=fileWithAllCommands,append=T)
      } else{
        cat("titleToPut<-geneID\n", file=fileWithAllCommands,append=T)
      }
      if(usePng){
        cat(paste0("png(paste0(outputFolder,\"/",xaxisForGenes,"-\",titleToPut,\".png\"),width=",pngSize,",height=",pngSize,")\n"), file=fileWithAllCommands,append=T)
      } else {
        cat(paste0("pdf(paste0(outputFolder,\"/",xaxisForGenes,"-\",titleToPut,\".pdf\"),title=paste0(\"",xaxisForGenes,"-\",titleToPut),width=",pdfSize,",height=",pdfSize,")\n"), file=fileWithAllCommands,append=T)
      }
      cat("cmd<-\"ggplot(new.df, aes(",xaxisForGenes,",y)) +
geom_point(",gsub("\\\"","\\\\\"",param),",size=3)  +
labs(title = titleToPut) +
theme_grey(base_size = 20) +
",paste0("ylab(\\\"",ylab,"\\\")"),gsub("\\\"","\\\\\"",additionalParam),"\"\n", file=fileWithAllCommands,append=T)
      cat("print(eval(parse(text = cmd)))
dev.off()
}
}\n", file=fileWithAllCommands,append=T)  
  }
  if(exists("addGlobalHeatmap")){
    if(is.logical(addGlobalHeatmap)){
      if(addGlobalHeatmap){
        if(exists("keepGeneOrder")){
          if(!is.logical(keepGeneOrder)){
            cat("keepGeneOrder is not logical. It will be set to F.\n")
            keepGeneOrder<-F
          }
        } else {
          keepGeneOrder<-F
        }
        if(exists("clusterSamples")){
          if(!is.logical(clusterSamples)){
            cat("clusterSamples is not logical. It will be set to F.\n")
            clusterSamples<-F
          }
        } else {
          clusterSamples<-F
        }
        cols<-unique(which(colnames(factorizedSP)%in%c(xaxisForGenes,unlist(plotGenesPara))))
        cat("library(pheatmap)
        cols<-c(",paste(cols,collapse=","),")\n", file=fileWithAllCommands,append=T)
        if (length(cols) == 1) {
          annotForHM <- data.frame(factorizedSP[,cols])
          colnames(annotForHM) <- colnames(factorizedSP)[cols]
          rownames(annotForHM) <- samplesToPlot
          cat("annotForHM <- data.frame(factorizedSP[,cols])
colnames(annotForHM) <- colnames(factorizedSP)[cols]
rownames(annotForHM) <- samplesToPlot\n", file=fileWithAllCommands,append=T)
        } else if(length(cols)>1) {
          annotForHM<-factorizedSP[,cols]
          cat("annotForHM<-factorizedSP[,cols]\n", file=fileWithAllCommands,append=T)
        }
        if(exists("geneIDToAdd")){
          colWithName<-geneIDToAdd
          cat("colWithName<-geneIDToAdd\n", file=fileWithAllCommands,append=T)
        } else{
          colWithName<-colOfGeneID
          cat("colWithName<-colOfGeneID\n", file=fileWithAllCommands,append=T)
        }
        pdfSize<-max(7,5+0.15*nSamples)
        pngSize<-max(480,237+12*nSamples)
        if(usePng){
          png(paste0(outputFolder,"/HeatmapWithGenes.png"),width=pngSize,height=500+17*nrow(dfGene))
          cat("png(paste0(outputFolder,\"/HeatmapWithGenes.png\"),width=",pngSize,",height=",500+17*nrow(dfGene),")\n", file=fileWithAllCommands,append=T)
        } else {
          pdf(paste0(outputFolder,"/HeatmapWithGenes.pdf"),title="HeatmapWithGenes",onefile = TRUE,width=pdfSize,height=5+0.22*nrow(dfGene))
          cat("pdf(paste0(outputFolder,\"/HeatmapWithGenes.pdf\"),title=\"HeatmapWithGenes\",onefile = TRUE,width=",pdfSize,",height=",5+0.22*nrow(dfGene),")\n", file=fileWithAllCommands,append=T)
        }
        epsilon=0.0000001
        df.sub<-dataToPlot[match(dfGene[,1],expressionDF[,colOfGeneID]),]
        if(anyDuplicated(expressionDF[expressionDF[,colOfGeneID]%in%dfGene[,1],colOfGeneID])>0){
          cat("Warning: at least one gene id in the file correspond to multiple rows in the expression file.\n Only the first row with the good gene id will be used.\n")
        }
        if(useLogExpression){
          breaksListAbs<-c(seq(0,12,length.out = 100),max(max(df.sub),12)+epsilon)
        } else {
          breaksListAbs<-c(seq(0,max(df.sub),length.out = 100),max(df.sub)+epsilon)
        }
        pheatmap(df.sub,
                 labels_row=expressionDF[match(dfGene[,1],expressionDF[,colOfGeneID]),colWithName],
                 cluster_rows=!keepGeneOrder,
                 cluster_cols=clusterSamples,
                 annotation=annotForHM,
                 annotation_colors = fixedColors,
                 breaks = breaksListAbs,
                 main=ylab,
                 cellheight = 16)
        dev.off()
        cat("epsilon=0.0000001
          df.sub<-dataToPlot[match(dfGene[,1],expressionDF[,colOfGeneID]),]\n", file=fileWithAllCommands,append=T)
        if(useLogExpression){
          cat("breaksListAbs<-c(seq(0,12,length.out = 100),max(max(df.sub),12)+epsilon)\n", file=fileWithAllCommands,append=T)
        } else {
          cat("breaksListAbs<-c(seq(0,max(df.sub),length.out = 100),max(df.sub)+epsilon)\n", file=fileWithAllCommands,append=T)
        }
        cat("pheatmap(df.sub,
                   labels_row=expressionDF[match(dfGene[,1],expressionDF[,colOfGeneID]),colWithName],
                   cluster_rows=",!keepGeneOrder,",
                   cluster_cols=",clusterSamples,",
                   annotation=annotForHM,
                   annotation_colors = fixedColors,
                   breaks = breaksListAbs,
                   main=\"",ylab,"\",
                   cellheight = 16)
          dev.off()\n", file=fileWithAllCommands,append=T)
        
      }
    } else {
      cat("addGlobalHeatmap is not logical. No heatmap will be plotted.\n")
    }
  }
  } else {
    cat("The file provided as fileWithGenes does not exists.\n")
  }
}
