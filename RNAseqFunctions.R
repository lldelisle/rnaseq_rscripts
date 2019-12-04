mergeCounts_function<-function(samplesPlanDF){
  #This function will generate one table with all the count values
  #The files must be htseqCount_shortName.txt in the working directory or specified in the samplesplan
  if(!"htseq_count_file"%in%colnames(samplesPlanDF)){
    samplesPlanDF$htseq_count_file<-sapply(samplesPlanDF$sample,function(s){paste0(getwd(),"/htseqCount_",s,".txt")})
  }

  samplesPlanDF$htseq_count_file_exists<-sapply(samplesPlanDF$htseq_count_file,file.exists)
  
  for(i in which(!samplesPlanDF$htseq_count_file_exists)){
    cat("The file with FPKM values for",samplesPlanDF$sample[i],":",samplesPlanDF$htseq_count_file[i],"does not exist. The sample is omitted.\n")
  }

  samplesPlanDF<-subset(samplesPlanDF,htseq_count_file_exists)
  
  #They are all merged
  for( i in 1:nrow(samplesPlanDF)){
    NexthtseqCounts<-read.delim(samplesPlanDF$htseq_count_file[i],h=T,stringsAsFactors=F)
    colnames(NexthtseqCounts)<-c("Ens_ID",samplesPlanDF$sample[i])
    if(!exists("htseqCounts")){
      htseqCounts<-NexthtseqCounts
    } else {
      htseqCounts<-merge(htseqCounts,NexthtseqCounts,all=T)
    }
  }
  #I remove the lines with stats:
  #Begining by __ from htseq-count
  #Begining by N_ from STAR
  htseqCounts<-htseqCounts[!grepl("^[N_]_",htseqCounts$Ens_ID),]
  return(htseqCounts)
}

mergeFPKM_function<-function(samplesPlanDF,sumDup=F,colToKeep=c("gene_id","gene_short_name","locus")){
  #This function will generate one table with all the FPKM values and the name of the gene and the locus
  #If sumDup=F you can have different lines for one gene_id
  #If sumDup=T all the values with the same gene_id will be summed up
  #The files must be FPKM_shortName.txt in the working directory or specified in the samplesplan
  
  if(!"cufflinks_file"%in%colnames(samplesPlanDF)){
    samplesPlanDF$cufflinks_file<-sapply(samplesPlanDF$sample,function(s){paste0(getwd(),"/FPKM_",s,".txt")})
  }

  samplesPlanDF$cufflinks_file_exists<-sapply(samplesPlanDF$cufflinks_file,file.exists)
  
  for(i in which(!samplesPlanDF$cufflinks_file_exists)){
    cat("The file with FPKM values for",samplesPlanDF$sample[i],":",samplesPlanDF$cufflinks_file[i],"does not exist. The sample is omitted.\n")
  }

  samplesPlanDF<-subset(samplesPlanDF,cufflinks_file_exists)
    
  #They are all merged
  for( i in 1:nrow(samplesPlanDF)){
    NextFPKMCuff<-read.delim(samplesPlanDF$cufflinks_file[i],h=T,stringsAsFactors=F)[,c(colToKeep,"FPKM")]
    colnames(NextFPKMCuff)<-c(colToKeep,paste0("FPKM_",samplesPlanDF$sample[i]))
    if(!exists("FPKMCuff")){
      FPKMCuff<-NextFPKMCuff
    } else {
      FPKMCuff<-merge(FPKMCuff,NextFPKMCuff,all=T)
    }
  }
  
  if(sumDup){
    #If sumDup is asked the lines with the same gene_id will be summed
    simplified<-aggregate(FPKMCuff[,4:ncol(FPKMCuff)],by=list(FPKMCuff$gene_id),FUN=sum)
    simplified<-data.frame(gene_id=simplified$Group.1,gene_short_name=FPKMCuff$gene_short_name[match(simplified$Group.1,FPKMCuff$gene_id)],
                           locus=FPKMCuff$locus[match(simplified$Group.1,FPKMCuff$gene_id)],simplified[,2:ncol(simplified)], stringsAsFactors=FALSE)
    #For the gene_id that were duplicated the locus is modified to take into account all transcripts
    dup<-FPKMCuff$gene_id[duplicated((FPKMCuff$gene_id))]
    dupFPKM<-FPKMCuff[FPKMCuff$gene_id%in%dup,c("gene_id","locus")]
    ldup<-matrix(unlist(strsplit(unlist(strsplit(as.character(dupFPKM$locus),":")),"-")),ncol=3,byrow=T)
    minLoc<-aggregate(as.numeric(ldup[,2]),by=list(dupFPKM$gene_id),FUN=min)
    maxLoc<-aggregate(as.numeric(ldup[,3]),by=list(dupFPKM$gene_id),FUN=max)
    chrLoc<-ldup[match(minLoc$Group.1,dupFPKM$gene_id),1]
    locusToPutBack<-data.frame(ensID=minLoc$Group.1,paste(chrLoc,paste(minLoc$x,maxLoc$x[match(minLoc$Group.1,maxLoc$Group.1)],sep="-"),sep=":"), stringsAsFactors=FALSE)
    simplified$locus[match(locusToPutBack$ensID,simplified$gene_id)]<-locusToPutBack[,2]
    return(simplified)
  } else{
    return(FPKMCuff)
  }
}

##From Anouk
normalizeHKRank<-function(FPKMCuff,nbgenes,chrToRm=c("chrM")){
  expdata<-FPKMCuff[,grep("FPKM_",colnames(FPKMCuff))]
  haveRowNames<-F
  if(all(c("gene_id","locus")%in%colnames(FPKMCuff))){
    rownames(expdata)<-paste(FPKMCuff$gene_id,FPKMCuff$locus,sep='__')
    haveRowNames<-T
  } else {
    rownames(expdata)<-paste(1:nrow(FPKMCuff),rep("NA",nrow(FPKMCuff)),sep='__')
  }
  
  if(!all(is.na(chrToRm))){
    #I check that the rownames have info:
    if(!haveRowNames){
      cat("The cufflinks table did not have the gene_id and locus info. It is impossible to remove the specified chr.\n")
    } else {
      #I remove the chrToRm genes
      linesWithChrToRm<-unlist(sapply(chrToRm,function(c){grep(paste0("__",c,":"),rownames(expdata))}))
      if(length(linesWithChrToRm)>0){
        cat(length(linesWithChrToRm),"lines removed.\n")
        expdata<-expdata[-linesWithChrToRm,]
      } else{
        cat("There is no gene in the chromosomes",chrToRm,"\n")
      }    
    }
  }
  
  ### get genes that are expressed in all samples
  expdata.nonzero=expdata[which(apply(expdata,1,min)>0 ),]
  
  ## transform the expression levels into ranks
  
  expdata.ranks=apply(expdata.nonzero,2,rank)
  
  ## compute the variance of the ranks for each gene
  
  expdata.ranks.var=apply(expdata.ranks,1,var,na.rm=T)
  
  ## rank the genes according to their variance
  
  expdata.nonzero.consrank=rank(expdata.ranks.var)
  
  ## compute the median rank over all samples, for each gene
  
  median.rank=apply(expdata.ranks,1,median,na.rm=T)
  
  ## get the genes which have a median rank in the 25%-75% range
  
  interquartile=median.rank>(0.25*length(median.rank)) & median.rank<(0.75*length(median.rank))
  
  ## get the house-keeping genes: the nbgenes genes with the most conserved ranks, among those that fall in the interquartile range
  
  hkgenes=names(sort(expdata.nonzero.consrank[interquartile])[1:nbgenes])
  
  ## compute the normalization coefficient for each sample =  the median of the RPKM of the hkgenes in that sample
  
  normcoeff=apply(expdata.nonzero[hkgenes,],2,median,na.rm=T)
  
  ## we want to bring all of the medians at the average of the medians (computed on all expressed genes)
  
  ##  normcoeff=normcoeff/mean(apply(expdata.nonzero,2,median))
  
  normcoeff=normcoeff/mean(normcoeff)
  
  ## finally, normalize the data
  
  hkgenesEnsID<-unlist(strsplit(hkgenes,"__"))[seq(1,2*length(hkgenes),2)]
    
  if(haveRowNames && "gene_short_name"%in%colnames(FPKMCuff)){
    hkgenesName<-FPKMCuff$gene_short_name[match(hkgenesEnsID,FPKMCuff$gene_id)]
  } else {
    hkgenesName<-NA
  }
  
  FPKMCuff.norm<-data.frame(FPKMCuff[,grep("FPKM_",colnames(FPKMCuff),invert=T)],t(t(FPKMCuff[,grep("FPKM_",colnames(FPKMCuff))])/normcoeff))
  
  results=list("hkgenesName"=hkgenesName,"hkgenesENSID"=hkgenesEnsID,"normcoeff"=normcoeff,"normData"=FPKMCuff.norm)
  
}

simplifyDF<-function(df,samples,unusedColnames=c("htseq_count_file","cufflinks_file"),keepFactorsWithOneValue=F){
  newdf<-df
  rownames(newdf)<-newdf$sample
  newdf<-newdf[samples,-which(colnames(newdf)%in%unusedColnames)]
  for(cn in colnames(newdf)){
    uniqVal<-unique(newdf[,cn])
    if(length(uniqVal)==1){
      if(keepFactorsWithOneValue){
        newdf[,cn]<-factor(newdf[,cn],levels=uniqVal)
      } else {
        newdf[,cn]<-NULL
      }
    } else {
      newdf[,cn]<-factor(newdf[,cn],levels=uniqVal)
    }
  }
  return(newdf)
}

getStringFromListAndSP<-function(listOfPara,possibleFactors,possibleValues){
  param<-""
  for(v in possibleValues){
    if(v%in%names(listOfPara)){
      val<-listOfPara[[v]]
      if(val%in%possibleFactors){
        param<-paste0(param,", ",v,"=",val)
      } else {
        cat(val,"will not be used as",v,"because it is not part of the samplesPan table.\n")
      }
    }
  }
  if(length(setdiff(names(listOfPara),possibleValues))>0){
    cat("The parameters ",setdiff(names(listOfPara),possibleValues)," are not usable in this plot.\n")
  }
  return(param)
}

isValidColor <- function(colorname){
  if(is.numeric(colorname)){
    return(TRUE)
  }
  if(colorname %in% colors()){
    return(TRUE)
  }
  if(is.character(colorname)){
    if(nchar(colorname)==7 || nchar(colorname)==9){
      if(substr(colorname,1,1)=="#"){
        #I should do other checks
        return(TRUE)
      }
    }
  }
  return(FALSE)
}
