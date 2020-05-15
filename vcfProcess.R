#' Filter VCF files and convert to FASTA and CSV files
#' @param inputfile Single or multisample VCF file
#' @param outputfile Prefix for output files
#' @param caller Variant calling software used ("SAMtools" or "GATK", default= "SAMtools")
#' @param no.Cores Number of CPU cores to use - if >1, script will run some parallelization 
#' @param samples2remove Removes named samples from multisample inputs (default=NULL)
#' @param indelProcess Process INDEL variants (requires indelProcess_vcfProcess.R)
#' @param DP_low Minimum read depth to consider call (default=5)
#' @param lowqual Minimum variant quality to consider call (defauly=20)
#' @param hetProp Proportion allele frequency at hSNPs to assign call (default=0.9)
#' @param misPercent Percentage of sites missing across samples to remove variant (default=80)
#' @param repeatfile .csv file with start and stop coordinates for regions to remove variants
#' @param excludeMix If true, will remove samples with a high likelihood of mixed infection (requires MixInfect_vcfProcess.R)
#' @return filtered .vcf and .csv variant files, optional INDEL file and mixed infection files
#' @export

vcfProcess = function(inputfile,outputfile,caller="SAMtools",no.Cores=1,samples2remove=NULL,indelProcess=FALSE,DP_low=5,lowqual=20,hetProp=0.9,misPercent=80,repeatfile=NULL,excludeMix=FALSE){
  
  if (!require(stringr)){
    install.packages("stringr",repos = "http://cran.us.r-project.org")
    library(stringr)
  }
  if (!require(seqinr)){
    install.packages("seqinr",repos = "http://cran.us.r-project.org")
    library(seqinr)
  }
  options(stringsAsFactors = F)
  
  vcf<-read.table(inputfile)
  header_input<-as.matrix(read.table(inputfile,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",header_input)==TRUE)
  header<-as.data.frame(header_input[1:end_head-1])
  names<-unlist(strsplit(header_input[end_head],"\t"))
  head_start<-names[1:9]
  format<-which(names=='FORMAT')
  names<-names[10:length(names)]
  
  ###### Remove named samples
  if (!is.null(samples2remove)){
    remove<-samples2remove
    vcf<-vcf[,-(which(is.element(names,remove))+9)]
    names<-names[-which(is.element(names,remove))]
  }
  
  #### Order of format fields
  st<-unlist(str_split(vcf[1,format], ":"))
  GT<-which(st=='GT')
  DP<-which(st=='DP')
  PL<-which(st=='PL')
  AD<-which(st=='AD')
  if (length(AD)==0){
    AD<-which(st=='DP')
    hetProp<-1
  } # if no AD field, hetProp set to 1
  
  ####### Remove low quality variants
  lowqual_mat<-vcf[which((as.integer(vcf[,6]) < lowqual) ==T),]
  colnames(lowqual_mat)<-c(head_start,names)
  vcf<-vcf[which((as.integer(vcf[,6]) < lowqual) == F),]
  if (nrow(lowqual_mat)!=0){
    write.csv(lowqual_mat,file=paste0(outputfile,"_lowqualVariants.csv"),row.names = F)
  } else {print("No low quality variants")}
  
  ####### Remove overlapping variants from GATK calls
  if (length(grep("*",vcf[,5],fixed = T))>1 & caller=="GATK"){
    overlap_mat<-vcf[grep("*",vcf[,5],fixed = T),]
    colnames(overlap_mat)<-c(head_start,names)
    vcf<-vcf[-grep("*",vcf[,5],fixed = T),]
    if (nrow(overlap_mat)!=0){
      write.csv(overlap_mat,file=paste0(outputfile,"_overlapVariants.csv"),row.names = F)
    } else {print("No overlap variants")}
  }
  
  ######## Process INDELs
  ind<- lapply(1:nrow(vcf), function(i){
    length(unlist(strsplit(vcf[i,4],"")))>1 || length(unlist(strsplit(unlist(strsplit(vcf[i,5],","))[1],"")))>1
  }
  )
  indels<-vcf[which(ind==TRUE),]
  indelPos<-indels[,2]
  colnames(indels)<-c(head_start,names)
  vcf<-vcf[which(ind==FALSE),]
  write.csv(indels,file=paste0(outputfile,"_InDels.csv"),row.names = F)
  if (indelProcess){
    if (nrow(indels)!=0){
      indelfile<-indelProcess(indels,names,GT,AD,DP,hetProp,DP_low,outputfile,repeatfile.present,repeatfile)
    } else {print("No indels")
    }
  }
  
  #### Assign alleles
  genotype<-data.frame(vcf[,4:5])
  AD_comp<-data.frame(vcf[,4:5])
  read<-data.frame(vcf[,4:5])
  for (i in 10:ncol(vcf)){
    split_fieldsGT<-sapply(strsplit(vcf[,i],split = ":"), function(x) x[GT])
    split_fieldsAD<-sapply(strsplit(vcf[,i],split = ":"), function(x) x[AD])
    split_fieldsDP<-sapply(strsplit(vcf[,i],split = ":"), function(x) x[DP])
    genotype<-cbind(genotype,split_fieldsGT)
    AD_comp<-cbind(AD_comp,split_fieldsAD)
    read<-cbind(read,split_fieldsDP)
  }
  mixed_sites<-matrix(c("R","GA","AG","M","CA","AC","W","TA","AT","Y","TC","CT","S","GC","CG","K","TG","GT"),ncol=3,byrow = T)
  alt_allele<-strsplit(genotype[,2],split = ",")
  max_length <- max(unlist(lapply(strsplit(genotype[,2],split = ","), length)))
  
  if (no.Cores>1){
    if (!require(foreach)){
      install.packages("foreach",repos = "http://cran.us.r-project.org")
      library(foreach)
    }
    if (!require(doMC)){
      install.packages("doMC",repos = "http://cran.us.r-project.org")
      library(doMC)
    }
    registerDoMC(no.Cores)
    output<-foreach(i=1:nrow(genotype)) %dopar% {
      df<-data.frame(matrix("N",ncol=ncol(genotype)))
      df[1,as.numeric(which(unlist(genotype[i,])=="./."))]<-"N" # Missing
      df[1,as.numeric(which(unlist(genotype[i,])=="0/0"))]<-genotype[i,1] # Ref call
      for (j in 1:max_length){
        df[1,as.numeric(which(unlist(genotype[i,])==paste0(j,"/",j)))]<-unlist(strsplit(genotype[i,2],","))[j]
      } # Alt call
      alleles<-c(genotype[i,1],unlist(strsplit(genotype[i,2],",")))
      for (j in 0:(max_length-1)){
        for (k in (j+1):max_length){
          mixed<-which(genotype[i,]==paste0(j,"/",k))
          if (length(mixed)>0){
            for (mix in 1:length(mixed)){
              if (hetProp<1){
                AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mixed[mix]],split = ",")))
                if (sum(AD_vec)>0 && AD_vec[(j+1)]/sum(AD_vec)>hetProp){
                  df[1,mixed[mix]]<-alleles[j+1]
                } else if (sum(AD_vec)>0 && AD_vec[k+1]/sum(AD_vec)>hetProp){
                  df[1,mixed[mix]]<-alleles[k+1]
                } else {
                  df[1,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]}
              } else {
                df[1,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]}
            } 
          }
        }
      }
      df
    }
    output<- matrix(unlist(output), nrow=length(output), byrow=T)
    colnames(output)=colnames(genotype)
  } else {
    output=matrix("N",nrow(genotype),ncol(genotype))
    rownames(output)=rownames(genotype)
    colnames(output)=colnames(genotype)
    for (i in 1:nrow(genotype)){
      output[i,as.numeric(which(unlist(genotype[i,])=="./."))]<-"N" # Missing
      output[i,as.numeric(which(unlist(genotype[i,])=="0/0"))]<-genotype[i,1] # Ref call
      for (j in 1:max_length){
        output[i,as.numeric(which(unlist(genotype[i,])==paste0(j,"/",j)))]<-unlist(strsplit(genotype[i,2],","))[j]
      } # Alt call
      alleles<-c(genotype[i,1],unlist(strsplit(genotype[i,2],",")))
      for (j in 0:(max_length-1)){
        for (k in (j+1):max_length){
          mixed<-which(genotype[i,]==paste0(j,"/",k))
          if (length(mixed)>0){
            for (mix in 1:length(mixed)){
              if (hetProp<1){
                AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mixed[mix]],split = ",")))
                if (sum(AD_vec)>0 && AD_vec[(j+1)]/sum(AD_vec)>hetProp){
                  output[i,mixed[mix]]<-alleles[j+1]
                } else if (sum(AD_vec)>0 && AD_vec[k+1]/sum(AD_vec)>hetProp){
                  output[i,mixed[mix]]<-alleles[k+1]
                } else {
                  output[i,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]}
              } else {
                output[i,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]}
            } 
          }
        }
      }
    }
  }
  output<-output[,3:ncol(output)] 
  
  ######## Mark low read positions as N
  read<-as.matrix(read[,3:ncol(read)])
  class(read)<-"numeric"
  snprd<-read<DP_low
  lowread<-which(snprd == 'TRUE') 
  output[lowread]<-'N'
  
  ##### Remove invariant sites 
  new21<-as.matrix(cbind(vcf[,4],output))
  new22<-matrix(0,nrow(new21),ncol(new21))
  for (i in 1:nrow(new21)){
    for (j in 2:ncol(new21)){
      if (new21[i,j]==new21[i,1]){
        new22[i,j]='TRUE'}
      else  {
        new22[i,j]='FALSE'
      }
    }
  }
  
  new23<-new22[,2:ncol(new22)]
  rd<- lapply(1:nrow(new23), function(i){
    all(as.logical(new23[i,]))
  }
  )
  called<-do.call(rbind,rd)
  output<-output[which(called=='FALSE'),] 
  vcf<-vcf[which(called=='FALSE'),]
  
  ####### Remove SNPs in repeat/specified regions
  if (!is.null(repeatfile)){
    rep<-as.data.frame(read.csv(repeatfile,header=TRUE))
    res1=vector()
    for (i in 1:nrow(rep)){
      re<-c(rep[i,1]:rep[i,2])
      res1<-c(res1,re)   
    }
    res<-is.element(as.numeric(vcf[,2]),res1)
    rep<-cbind(vcf[which(res == 'TRUE'),2],vcf[which(res == 'TRUE'),4],output[which(res=='TRUE'),])
    colnames(rep)<-c("Position","Reference",names)
    output<-output[which(res == 'FALSE'),]
    vcf<-vcf[which(res == 'FALSE'),]
    if (nrow(rep)!=0){
      write.csv(rep, file=paste0(outputfile,"_SNPsinRepRegions.csv"),row.names = F)
    } else {print("No repeat or mobile element SNPs")
    }
  }
                           
  #### Test for mixed infection, remove mixed samples and reassign mixed calls 
  if (excludeMix){
    result<-MixInfect(vcf,names,output,hetProp,outputfile,format)
    vcf<-result$vcf
    names<-result$names
    output<-result$output
  }
  
  ##### Remove SNPs > misPercent missing data
  percent_missing<-integer(nrow(output))
  for (i in 1:nrow(output)){
    percent_missing[i]<-((length(which(output[i,]=="?"))+length(which(output[i,]=="N")))/ncol(output))*100
  }
  output<-output[which(percent_missing<=misPercent),]
  vcf<-vcf[which(percent_missing<=misPercent),]
  
  ########## Write VCF file of processed SNPs
  vfile<-as.data.frame(apply(vcf, 1, paste, collapse="\t"))
  rownames(vfile)<-NULL
  headrow<-paste0(c(head_start,names),collapse = "\t")
  names(headrow)<-names(header)
  names(vfile)<-names(header)
  outvcf<-rbind(header,headrow,vfile)
  write.table(outvcf, file= paste0(outputfile,"_processed.vcf"),row.names =FALSE,sep="\t",quote=FALSE,col.names=FALSE)
  
  ########### CSV file with SNP genotypes only with reference sequence and SNP position 
  colnames(output)<-names
  Reference<-as.character(vcf[,4])
  Position<-vcf[,2]
  forcsv<-cbind(Position,Reference,output)
  write.csv(forcsv,file=paste0(outputfile,"_with_ref.csv"),row.names=FALSE)
  
  ########### Fasta file of all isolates with reference
  namesfasta<-c("Reference",names)
  output_fast<-t(output)
  forfastaref<-rbind(Reference,output_fast)
  forfastaref<-as.list(apply(forfastaref, 1, paste, collapse=""))
  write.fasta(forfastaref,namesfasta,nbchar=60,file.out=paste0(outputfile,"_with_ref.fasta"),open="w")
  
  ########### FASTA file of samples only
  mat<-cbind(as.character(vcf[,4]),output)
  new22<-matrix(0,nrow(mat),ncol(mat))
  for (i in 1:nrow(mat)){
    for (j in 2:ncol(mat)){
      if (mat[i,j]==mat[i,1] | mat[i,j]=="N"){
        new22[i,j]='FALSE'}
      else  {
        new22[i,j]='TRUE'
      }
    }
  }
  new22<-new22[,2:ncol(new22)]
  df<- lapply(1:nrow(new22), function(i){
    all(as.logical(new22[i,]))
  }
  )
  out<-do.call(rbind,df)
  output<-output[which(out=='FALSE'),] 
  vcf<-vcf[which(out=='FALSE'),]
  Reference<-as.character(vcf[,4])
  Position<-vcf[,2]
  forcsv<-cbind(Position,Reference,output)
  write.csv(forcsv,file=paste0(outputfile,".csv"),row.names=FALSE)
  output_fast<-t(output)
  forfasta<-as.list(apply(output_fast, 1, paste, collapse=""))
  write.fasta(forfasta,names,nbchar=60,file.out=paste0(outputfile,".fasta"),open="w")
}
