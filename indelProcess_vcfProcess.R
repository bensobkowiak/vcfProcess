#' Process INDEL file from vcfProcess
#' @param indels INDEL file from vcfProcess
#' @param names Sample names
#' @param output Genotype file
#' @param LowCov Minimum read depth at site to accurately call allele frequency
#' @param BICoutput Option to produce output file with probability of each sample being mixed
#' @return vcf, names and output genotype files with mixed infections removed
#' @export

indelProcess = function(indels,indel_header,names,hetProp,DP_low,outputfile,repeatfile,misPercent){
  
  if (!require(stringr)){
    install.packages("stringr",repos = "http://cran.us.r-project.org")
    library(stringr)
  }
  if (!require(ape)){
    install.packages("ape",repos = "http://cran.us.r-project.org")
    library(ape)
  }
  
  indelvcf<-indels
  indel_header<-indel_header
  names_indel<-names
  indel_alleles<-data.frame(indelvcf[,4:5])
  st<-unlist(str_split(indelvcf[1,format], ":"))
  GT<-which(st=='GT')
  DP<-which(st=='DP')
  AD<-which(st=='AD')
  hetProp<-hetProp
  DP_low<-DP_low
  outputfile<-outputfile
  misPercent<-misPercent
  
  indel_genotype<-data.frame(indelvcf[,4])
  indel_AD_comp<-data.frame(indelvcf[,4])
  indel_read<-data.frame(indelvcf[,4])
  for (i in 10:ncol(indelvcf)){
    split_fieldsGT<-sapply(strsplit(indelvcf[,i],split = ":"), function(x) x[GT])
    split_fieldsAD<-sapply(strsplit(indelvcf[,i],split = ":"), function(x) x[AD])
    split_fieldsDP<-sapply(strsplit(indelvcf[,i],split = ":"), function(x) x[DP])
    indel_genotype<-cbind(indel_genotype,split_fieldsGT)
    indel_AD_comp<-cbind(indel_AD_comp,split_fieldsAD)
    indel_read<-cbind(indel_read,split_fieldsDP)
  }
  
  indel_genotype<-as.matrix(indel_genotype[,2:ncol(indel_genotype)])
  indel_AD_comp<-as.matrix(indel_AD_comp[,2:ncol(indel_AD_comp)])
  output=matrix("N",nrow(indel_genotype),ncol(indel_genotype))
  rownames(output)=rownames(indel_genotype)
  colnames(output)=colnames(indel_genotype)
  max_length <- max(unlist(lapply(strsplit(indelvcf[,5],split = ","), length)))
  
  for (i in 1:nrow(indel_genotype)){
    output[i,as.numeric(which(unlist(indel_genotype[i,])=="./."))]<-"?"     # No call
    output[i,as.numeric(which(unlist(indel_genotype[i,])==".|."))]<-"?"     # No call
    output[i,as.numeric(which(unlist(indel_genotype[i,])=="0/0"))]<-0 # ref
    output[i,as.numeric(which(unlist(indel_genotype[i,])=="0|0"))]<-0 # ref
    for (j in 1:max_length){
      output[i,as.numeric(which(unlist(indel_genotype[i,])==paste0(j,"/",j) | unlist(indel_genotype[i,])==paste0(j,"|",j)))]<-j
    } 
    # het call
    for (j in 0:(max_length-1)){
      for (k in (j+1):max_length){
        mixed<-as.integer(which(indel_genotype[i,]==paste0(j,"/",k) | indel_genotype[i,]==paste0(j,"|",k)))
        if (length(mixed)>0){
          for (mix in 1:length(mixed)){
            if (hetProp<1){
              AD_vec<-as.numeric(unlist(strsplit(indel_AD_comp[i,mixed[mix]],split = ",")))
              if (sum(AD_vec)>0 && AD_vec[(j+1)]/sum(AD_vec)>hetProp){
                output[i,mixed[mix]]<-j
              } else if (sum(AD_vec)>0 && AD_vec[k+1]/sum(AD_vec)>hetProp){
                output[i,mixed[mix]]<-k
              } 
            }
          }
        }
      }
    }
  }
  colnames(output)<-names_indel
  
  ######## low read for all mark as '?' 
  
  indel_read<-as.matrix(indel_read[,2:ncol(indel_read)])
  suppressWarnings(class(indel_read)<-"numeric")
  lowread<-which((indel_read<DP_low) == 'TRUE') 
  output[lowread]<-'?'
  
  
  ###########REMOVE indels that have no variant (may have been included due to low read site)
  
  novariant<-character()
  for (i in 1:nrow(output)){
    novariant[i]<-any(output[i,]==1)
  }
  output<-output[which(novariant==TRUE),] 
  indelvcf<-indelvcf[which(novariant==TRUE),]
  
  ## Remove indels > missing percent in pop.
  percent_missing<-integer(nrow(output))
  for (i in 1:nrow(output)){
    percent_missing[i]<-(length(which(output[i,]=="?"))/ncol(output))*100
  }
  output<-output[which(percent_missing<=misPercent),]
  indelvcf<-indelvcf[which(percent_missing<=misPercent),]
  
  ####### REMOVE REPEAT AND MOBILE ELEMENT REGIONS
  if (!is.null(repeatfile)){
    rep<-as.data.frame(read.csv(repeatfile,header=TRUE))
    res=vector()
    for (i in 1:nrow(rep)){
      re<-c(rep[i,1]:rep[i,2])
      res<-c(res,re)   
    }
    q<-is.element(as.numeric(indelvcf[,2]),res)
    rep<-cbind(indelvcf[which(q == 'TRUE'),2],indelvcf[which(q == 'TRUE'),4],output[which(q=='TRUE'),])
    colnames(rep)<-c("Position","Reference",names)
    output<-output[which(q == 'FALSE'),]
    indelvcf<-indelvcf[which(q == 'FALSE'),]
    if (nrow(rep)!=0){
      write.csv(rep, file=paste0(outputfile,"_indels_removedRegions.csv"),row.names = F)
    } else {print("No repeat or mobile element SNPs")
    }
  }
  
  vfile<-as.data.frame(apply(indelvcf, 1, paste, collapse="\t"))
  rownames(vfile)<-NULL
  headrow<-paste0(c(indel_header,names_indel),collapse = "\t")
  names(headrow)<-names(indel_header)
  names(vfile)<-names(indel_header)
  outvcf<-rbind(indel_header,headrow,vfile)
  write.table(outvcf, file= paste0(outputfile,"_indelprocessed.vcf"),row.names =FALSE,sep="\t",quote=FALSE,col.names=FALSE)
  
  ## outfiles
  CSVoutput<-cbind(indelvcf[,2],output)
  colnames(CSVoutput)<-c("POS",names_indel)
  write.csv(CSVoutput,paste0(outputfile,"_indels.csv"),row.names = F)
  output_fast<-t(output)
  forfastaref<-as.list(apply(output_fast, 1, paste, collapse=""))
  write.fasta(forfastaref,names,nbchar=60,file.out=paste0(outputfile,"_indels.fasta"),open="w")
  indel_fasta<-read.fasta(paste0(outputfile,"_indels.fasta"))
  write.nexus.data(indel_fasta,paste0(outputfile,"_indels.nexus"),format = "standard")
}
