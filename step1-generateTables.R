options(scipen=999)
options(stringsAsFactors=F)
rm(list=ls())

library(tools)


if(length(commandArgs(TRUE))>0){
  f<-commandArgs(TRUE)[1]
} else{
  cat("Select the config file.\n")
  #Ask for the config file
  f <- file.choose()
}
#Check the necessary options
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

checkFunctionPath<-function(){
  if(!exists("RNAseqFunctionPath")){
    stop("The RNAseqFunctionPath is not provided.")
  } else {
    if(!file.exists(RNAseqFunctionPath)){
      stop("The file provided in RNAseqFunctionPath:",RNAseqFunctionPath," does not exists.")
    }
  }
}

if(exists("mergeCounts")){
  if(!is.logical(mergeCounts)){
    cat("mergeCounts is not boolean so the counts will not be merged.\n")
  } else {
    if(mergeCounts){
      checkFunctionPath()
      source(RNAseqFunctionPath)
      htseqCounts<-mergeCounts_function(samplesPlanDF)
      if(!exists("outputFolderForStep1")){
        outputFolderForStep1<-paste0(getwd(),"/")
      } else {
        dir.create(outputFolderForStep1,recursive=T,showWarnings=F)
      }
      write.table(htseqCounts, file=paste0(outputFolderForStep1,"/AllHTSeqCounts.txt"),sep='\t',row.names = F,quote=F)
      cat(paste0("The counts have been merged in ",outputFolderForStep1,"/AllHTSeqCounts.txt.\n"))
    }
  }
}

if(exists("subsetCounts")){
  if(!is.logical(subsetCounts)){
    cat("subsetCounts is not boolean so the counts will not be subset.\n")
  } else {
    if(subsetCounts){
      inputOK<-F
      if(exists("htseqCounts")){
        cat("Will subset the counts generated at the previous step.\n")
        inputOK<-T
        name<-"AllHTSeqCounts"
      } else {
        if(exists("initialTableWithCount")){
          if(file.exists(initialTableWithCount)){
            htseqCounts<-read.delim(initialTableWithCount)
            if(!"Ens_ID"%in%colnames(htseqCounts)){
              if(exists("geneIDColInInitialTable")){
                if(geneIDColInInitialTable%in%colnames(htseqCounts)){
                  colnames(htseqCounts)[colnames(htseqCounts)==geneIDColInInitialTable]<-"Ens_ID"
                  cat("Will subset the counts provided in the initialTableWithCount :",initialTableWithCount,".\n")
                  name<-basename(file_path_sans_ext(initialTableWithCount))
                  inputOK<-T
                } else {
                  cat("The table provided as initialTableWithCount : ",initialTableWithCount," does not contain a column called Ens_ID and the geneIDColInInitialTable specified is not part of the column names. The counts will not be subset.\n")
                }
              } else {
                cat("The table provided as initialTableWithCount : ",initialTableWithCount," does not contain a column called Ens_ID and the geneIDColInInitialTable (the column name with the ensembl ids in this table) has not been specified. The counts will not be subset.\n")           
              }
            } else {
              cat("Will subset the counts provided in the initialTableWithCount :",initialTableWithCount,".\n")
              name<-basename(file_path_sans_ext(initialTableWithCount))
              inputOK<-T
            }
          } else {
            cat("The file provided in the initialTableWithCount :",initialTableWithCount," does not exists. The counts will not be subset.\n")
          }
        } else {
          cat("subsetCounts was put as T but there is no table generated at previous step and there is no file provided as initialTableWithCount.The counts will not be subset.\n")
        }
      }
      if(inputOK){
        if(exists("genesToRmFromCounts")){
          if(file.exists(genesToRmFromCounts)){
            temp.df<-read.delim(genesToRmFromCounts,h=F)
            htseqCounts<-subset(htseqCounts,!Ens_ID%in%temp.df$V1)
            if(!exists("outputFolderForStep1")){
              outputFolderForStep1<-paste0(getwd(),"/")
            } else {
              dir.create(outputFolderForStep1,recursive=T,showWarnings=F)
            }
            write.table(htseqCounts, file=paste0(outputFolderForStep1,"/",name,"_subset.txt"),sep='\t',row.names = F,quote=F)
            cat(paste0("The counts have been subset in ",outputFolderForStep1,"/",name,"_subset.txt.\n"))
          } else {
            cat("The file specified as genesToRmFromCounts:",genesToRmFromCounts," does not exists. No subset will be done.\n")
          }
        } else {
          cat("There is no file specified as genesToRmFromCounts. No subset will be done.\n")
        }
      }
    }
  }
}


if(exists("mergeFPKM")){
  if(!is.logical(mergeFPKM)){
    cat("mergeFPKM is not boolean so the FPKM will not be merged.\n")
  } else {
    if(mergeFPKM){
      if(exists("oneLinePerEnsemblID")){
        if(!is.logical(oneLinePerEnsemblID)){
          cat("oneLinePerEnsemblID is not boolean so no simplification will be done.\n")
          oneLinePerEnsemblID<-F
        }
      } else {
        cat("oneLinePerEnsemblID is not specified so no simplification will be done.\n")
        oneLinePerEnsemblID<-F
      }
      if(!exists("mergeFPKM_function")){
        checkFunctionPath()
        source(RNAseqFunctionPath)
      }
      FPKMCuff<-mergeFPKM_function(samplesPlanDF,sumDup=oneLinePerEnsemblID)
      if(!exists("outputFolderForStep1")){
        outputFolderForStep1<-paste0(getwd(),"/")
      } else {
        dir.create(outputFolderForStep1,recursive=T,showWarnings=F)
      }
      if(oneLinePerEnsemblID){
        name<-"AllCufflinks_Simplified"
      } else {
        name<-"AllCufflinks"
      }
      write.table(FPKMCuff, file=paste0(outputFolderForStep1,"/",name,".txt"),sep='\t',row.names = F,quote=F)
      cat(paste0("The FPKM have been merged in ",outputFolderForStep1,"/",name,".txt.\n"))
    }
  }
}

if(exists("normFPKMWithAnoukMethod")){
  if(!is.logical(normFPKMWithAnoukMethod)){
    cat("normFPKMWithAnoukMethod is not boolean so the FPKM table will not be normalized.\n")
  } else {
    if(normFPKMWithAnoukMethod){
      inputOK<-F
      if(exists("FPKMCuff")){
        cat("Will normalize the FPKM generated at the previous step.\n")
        inputOK<-T
      } else {
        if(exists("initialTableWithFPKM")){
          if(file.exists(initialTableWithFPKM)){
            FPKMCuff<-read.delim(initialTableWithFPKM)
            if(length(grep("FPKM_",colnames(FPKMCuff)))==0){
              cat("The table provided in initialTableWithFPKM:",initialTableWithFPKM," does not have any column name begining by FPKM_. No normalization will be performed.\n")
            } else {
              cat("Will normalize the FPKM provided in the initialTableWithFPKM :",initialTableWithFPKM,".\n")
              name<-basename(file_path_sans_ext(initialTableWithFPKM))
              inputOK<-T
            }
          } else {
            cat("The file provided in the initialTableWithFPKM :",initialTableWithFPKM," does not exists. No normalization will be performed.\n")
          }
        } else {
          cat("normFPKMWithAnoukMethod was put as T but there is no table generated at previous step and there is no file provided as initialTableWithFPKM.  No normalization will be performed.\n")
        }
      }
      if(inputOK){
        if(!exists("chrToRemoveBeforeNormWithAnoukMethod")){
          chrToRemoveBeforeNormWithAnoukMethod<-NA
        }
        if(!exists("nbOfGenesWithAnoukMethod")){
          cat("No nbOfGenesWithAnoukMethod have been specified. 1000 wil be used.")
          nbOfGenesWithAnoukMethod<-1000
        } else {
          if(!is.numeric(nbOfGenesWithAnoukMethod)){
            cat("The value put in nbOfGenesWithAnoukMethod is not a number. 1000 will be used.")
            nbOfGenesWithAnoukMethod<-1000
          }
        }
        if(!exists("normalizeHKRank")){
          checkFunctionPath()
          source(RNAseqFunctionPath)
        }
        resl<-normalizeHKRank(FPKMCuff,nbOfGenesWithAnoukMethod,chrToRemoveBeforeNormWithAnoukMethod)
        if(!exists("outputFolderForStep1")){
          outputFolderForStep1<-paste0(getwd(),"/")
        } else {
          dir.create(outputFolderForStep1,recursive=T,showWarnings=F)
        }
        write.table(resl[["normData"]],file=paste0(outputFolderForStep1,"/",name,"_norm.txt"),sep='\t',row.names = F,quote=F)
        cat(paste0("The FPKM have been normalized in ",outputFolderForStep1,"/",name,"_norm.txt.\n"))
        write.table(data.frame(sample=names(resl$normcoeff),coef=resl$normcoeff),file=paste0(outputFolderForStep1,"/",name,"_coeffsUsedForNorm.txt"),sep='\t',row.names = F,quote=F)
        cat(paste0("The coefficients used to normalize are in ",outputFolderForStep1,"/",name,"_coeffsUsedForNorm.txt.\n"))
        if(exists("keepGenesUsedForNorm")){
          if(is.logical(keepGenesUsedForNorm)){
            if(keepGenesUsedForNorm){
              write.table(data.frame(Ens_ID=resl[["hkgenesENSID"]],gene_name=resl[["hkgenesName"]]),file=paste0(outputFolderForStep1,"/",name,"_HKGenesUsedForNorm.txt"),sep='\t',row.names = F,quote=F)
              cat(paste0("The genes used to normalize are written in ", outputFolderForStep1,"/",name,"_HKGenesUsedForNorm.txt.\n"))
            }
          }
        }
      }
    }
  }
}
