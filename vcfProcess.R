#' Filter VCF files and convert to FASTA and CSV files
#' @param inputfile Single or multisample SNP or SNP/INDEL VCF file
#' @param outputfile Prefix for output files
#' @param indelfile Name of file containing INDELs file if separate to SNPs
#' @param no.Cores Number of CPU cores to use - if >1, script will run some parallelization 
#' @param samples2remove Samples to remove from multisample inputs, can be string or .txt file in single column (default=NULL)
#' @param samples2include Samples to include from multisample inputs, can be string or .txt file in single column (default=NULL (keep all samples))
#' @param processIndel Process INDEL variants (requires indelProcess_vcfProcess.R)
#' @param DP_low Minimum read depth to consider call (default=5)
#' @param lowqual Minimum variant quality to consider call (defauly=20)
#' @param hetProp Proportion allele frequency at hSNPs to assign call (default=0.9)
#' @param hetasN Code hSNPs as 'N' (TRUE) or by FASTA nucleic acid code (FALSE)
#' @param misPercent Percentage of sites missing across samples to remove variant (default=90)
#' @param repeatfile .csv file with start and stop coordinates for regions to remove variants
#' @param disINDEL Remove SNPs within int distance of INDELs 
#' @param MixInfect If TRUE, run MixInfect to test for mixed infections, if excludeMix = TRUE, will remove samples with a high likelihood of mixed infection (requires MixInfect_vcfProcess.R)
#' @return filtered .vcf and .csv variant files, optional INDEL and mixed infection files
#' @export

if (!require(stringr)){
  install.packages("stringr",repos = "http://cran.us.r-project.org")
  library(stringr)
}
if (!require(seqinr)){
  install.packages("seqinr",repos = "http://cran.us.r-project.org")
  library(seqinr)
}
if (!require(foreach)){
  install.packages("foreach",repos = "http://cran.us.r-project.org")
  library(foreach)
}
if (!require(doMC)){
  install.packages("doMC",repos = "http://cran.us.r-project.org")
  library(doMC)
}

options(stringsAsFactors = F)

vcfProcess = function(inputfile,outputfile="output",
                      indelfile=NULL,no.Cores=1,samples2remove=NULL,samples2include=NULL,filter=TRUE,
                      processIndel=FALSE,DP_low=5,lowqual=20,hetProp=0.9,hetasN=TRUE,misPercent=90,
                      repeatfile=NULL,disINDEL=NULL,MixInfect=TRUE,excludeMix=FALSE){
  

  source("~/Documents/Scripts/MixInfect_vcfProcess.R")
  source("~/Documents/Scripts/indelProcess_vcfProcess.R")
  
  vcf<-read.table(inputfile)
  if (!is.null(indelfile)){
    indelvcf<-read.table(indelfile)
    indel_header_input<-as.matrix(read.table(indelfile,comment.char=" ",sep="\n"))
    end_head_indel<-which(grepl("#CHROM",indel_header_input)==TRUE)
    indel_header<-as.data.frame(indel_header_input[1:end_head_indel-1])
    vcf<-rbind(vcf,indelvcf)
  }
  header_input<-as.matrix(read.table(inputfile,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",header_input)==TRUE)
  header<-as.data.frame(header_input[1:end_head-1])
  names<-unlist(strsplit(header_input[end_head],"\t"))
  head_start<-names[1:9]
  format<-which(names=='FORMAT')
  names<-names[10:length(names)]
  
  ###### Remove or only include named samples
  if (!is.null(samples2remove)){
    remove<-read.table(samples2remove,check.names = F)[,1]
    vcf<-vcf[,-(which(is.element(names,remove))+9)]
    names<-names[-which(is.element(names,remove))]
  }
  if (!is.null(samples2include)){
    include<-read.table(samples2include,check.names = F)[,1]
    vcf<-cbind(vcf[,1:9],vcf[,(which(is.element(names,include))+9)])
    names<-names[which(is.element(names,include))]
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
  
  ###### Remove variants if failed filter
  if (filter){
    filterCol<-which(head_start=="FILTER")
    failedQC_mat<-vcf[which(vcf[,filterCol] != "PASS"),]
    colnames(failedQC_mat)<-c(head_start,names)
    vcf<-vcf[which(vcf[,filterCol] == "PASS"),]
    if (nrow(failedQC_mat)!=0){
      write.csv(failedQC_mat,file=paste0(outputfile,"_filteredVariants.csv"),row.names = F)
    } else {print("No filtered variants")}
  }
  
  ####### Remove low quality variants
  lowqual_mat<-vcf[which((as.integer(vcf[,6]) < lowqual) ==T),]
  colnames(lowqual_mat)<-c(head_start,names)
  vcf<-vcf[which((as.integer(vcf[,6]) < lowqual) == F),]
  if (nrow(lowqual_mat)!=0){
    write.csv(lowqual_mat,file=paste0(outputfile,"_lowqualVariants.csv"),row.names = F)
  } else {print("No low quality variants")}
  
  ####### Remove spanning deletion variants
  if (length(grep("*",vcf[,5],fixed = T))>1){
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
  write.csv(indels,file=paste0(outputfile,"_InDels_prefiltered.csv"),row.names = F)
  if (processIndel){
    if (nrow(indels)!=0){
      if (is.null(indelfile)){
        indel_header<-header
      }
      indelfile<-indelProcess(indels,indel_header,names,format,hetProp,DP_low,outputfile,repeatfile,misPercent)
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
  max_length <- max(unlist(lapply(strsplit(genotype[,2],split = ","), length)))
  
  if (no.Cores>1){
    registerDoMC(no.Cores)
    output<-foreach(i=1:nrow(genotype)) %dopar% {
      df<-data.frame(matrix("N",ncol=ncol(genotype)))
      df[1,as.numeric(which(unlist(genotype[i,])=="./."))]<-"N" # Missing
      df[1,as.numeric(which(unlist(genotype[i,])==".|."))]<-"N" # Missing
      df[1,as.numeric(which(unlist(genotype[i,])=="0/0"))]<-genotype[i,1] # Ref call
      df[1,as.numeric(which(unlist(genotype[i,])=="0|0"))]<-genotype[i,1] # Ref call
      for (j in 1:max_length){
        df[1,as.numeric(which(unlist(genotype[i,])==paste0(j,"/",j) | unlist(genotype[i,])==paste0(j,"|",j)))]<-unlist(strsplit(genotype[i,2],","))[j]
      } # Alt call
      alleles<-c(genotype[i,1],unlist(strsplit(genotype[i,2],",")))
      for (j in 0:(max_length-1)){
        for (k in (j+1):max_length){
          mixed<-which(genotype[i,]==paste0(j,"/",k) | genotype[i,]==paste0(j,"|",k))
          if (length(mixed)>0){
            for (mix in 1:length(mixed)){
              if (hetProp<1){
                AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mixed[mix]],split = ",")))
                if (sum(AD_vec)>0 && AD_vec[(j+1)]/sum(AD_vec)>hetProp){
                  df[1,mixed[mix]]<-alleles[j+1]
                } else if (sum(AD_vec)>0 && AD_vec[k+1]/sum(AD_vec)>hetProp){
                  df[1,mixed[mix]]<-alleles[k+1]
                } else if (hetasN==FALSE){
                  df[1,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]
                }
              } else if (hetasN==FALSE) {
                df[1,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]
              }
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
      output[i,as.numeric(which(unlist(genotype[i,])==".|."))]<-"N" # Missing
      output[i,as.numeric(which(unlist(genotype[i,])=="0/0"))]<-genotype[i,1] # Ref call
      output[i,as.numeric(which(unlist(genotype[i,])=="0|0"))]<-genotype[i,1] # Ref call
      for (j in 1:max_length){
        output[i,as.numeric(which(unlist(genotype[i,])==paste0(j,"/",j) | unlist(genotype[i,])==paste0(j,"|",j)))]<-unlist(strsplit(genotype[i,2],","))[j]
      } # Alt call
      alleles<-c(genotype[i,1],unlist(strsplit(genotype[i,2],",")))
      for (j in 0:(max_length-1)){
        for (k in (j+1):max_length){
          mixed<-which(genotype[i,]==paste0(j,"/",k) | genotype[i,]==paste0(j,"|",k))
          if (length(mixed)>0){
            for (mix in 1:length(mixed)){
              if (hetProp<1){
                AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mixed[mix]],split = ",")))
                if (sum(AD_vec)>0 && AD_vec[(j+1)]/sum(AD_vec)>hetProp){
                  output[i,mixed[mix]]<-alleles[j+1]
                } else if (sum(AD_vec)>0 && AD_vec[k+1]/sum(AD_vec)>hetProp){
                  output[i,mixed[mix]]<-alleles[k+1]
                } else if (hetasN==FALSE) {
                  output[i,mixed[mix]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alleles[j+1],alleles[k+1])) | is.element(mixed_sites[,3],paste0(alleles[j+1],alleles[k+1]))),1]}
              } else if (hetasN==FALSE) {
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
  suppressWarnings(class(read)<-"numeric")
  snprd<-read<DP_low
  output[which(snprd | is.na(snprd))]<-'?'
  
  ##### Remove invariant sites 
  new21<-as.matrix(cbind(vcf[,4],output))
  if (ncol(new21)>2){
  rows_to_keep <- which(rowSums(new21[, -1] != new21[, 1] & 
                                  new21[, -1] != "?", na.rm = TRUE) > 0)
  output<-output[rows_to_keep,] 
  vcf<-vcf[rows_to_keep,]
  } else {
    rows_to_keep <- which(new21[,2]!= new21[, 1] & 
                            new21[, 2] != "?")
    output<-as.data.frame(output[rows_to_keep])
    vcf<-vcf[rows_to_keep,]
  }

  
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
    output<-as.data.frame(output[which(res == 'FALSE'),])
    vcf<-vcf[which(res == 'FALSE'),]
    if (nrow(rep)!=0){
      write.csv(rep, file=paste0(outputfile,"_SNPsinRepRegions.csv"),row.names = F)
    } else {print("No repeat or mobile element SNPs")
    }
  }
  
  ####### Remove SNPs within int of INDEL
  if (!is.null(disINDEL) & length(indelPos)>0 & processIndel==TRUE){
    res1=vector()
    for (i in 1:length(indelPos)){
      re<-(indelPos[i]-disINDEL):(indelPos[i]+disINDEL)
      res1<-c(res1,re)   
    }
    res<-is.element(as.numeric(vcf[,2]),res1)
    inds<-cbind(vcf[which(res == 'TRUE'),2],vcf[which(res == 'TRUE'),4],output[which(res=='TRUE'),])
    colnames(inds)<-c("Position","Reference",names)
    output<-output[which(res == 'FALSE'),]
    vcf<-vcf[which(res == 'FALSE'),]
    if (nrow(inds)!=0){
      write.csv(inds, file=paste0(outputfile,"_SNPswithin",disINDEL,"bpINDELs.csv"),row.names = F)
    } else {print("No SNPS within distance of INDELs")
    }
  }
  
  #### Test for mixed infection, remove mixed samples and reassign mixed calls 
  if (MixInfect){
    result<-MixInfect(vcf,names,output,hetProp,outputfile,format,excludeMix)
    vcf<-result$vcf
    names<-result$names
    output<-result$output
  }
  
  ##### Remove SNPs > misPercent missing data
  percent_missing<-integer(nrow(output))
  for (i in 1:nrow(output)){
    percent_missing[i]<-((length(which(output[i,]=="?"))+length(which(output[i,]=="N")))/ncol(output))*100
  }
  output<-as.data.frame(output[which(percent_missing<=misPercent),])
  vcf<-vcf[which(percent_missing<=misPercent),]
  
  ########## Write VCF file of processed SNPs
  vfile<-as.data.frame(apply(vcf, 1, paste, collapse="\t"))
  rownames(vfile)<-NULL
  headrow<-paste0(c(head_start,names),collapse = "\t")
  names(headrow)<-names(header)
  names(vfile)<-names(header)
  outvcf<-rbind(header,headrow,vfile)
  write.table(outvcf, file= paste0(outputfile,"_SNPs_processed.vcf"),row.names =FALSE,sep="\t",quote=FALSE,col.names=FALSE)
  
  ########### CSV file with SNP genotypes only with reference sequence and SNP position 
  colnames(output)<-names
  Reference<-as.character(vcf[,4])
  Position<-vcf[,2]
  forcsv<-cbind(Position,Reference,output)
  write.csv(forcsv,file=paste0(outputfile,"_SNPs_with_ref.csv"),row.names=FALSE)
  
  ########### Fasta file of all isolates with reference
  namesfasta<-c("Reference",names)
  output_fast<-t(output)
  forfastaref<-rbind(Reference,output_fast)
  forfastaref<-as.list(apply(forfastaref, 1, paste, collapse=""))
  write.fasta(forfastaref,namesfasta,nbchar=60,file.out=paste0(outputfile,"_SNPs_with_ref.fasta"),open="w")
  
  ########### FASTA file of samples only
  new21<-as.matrix(cbind(vcf[,4],output))
  if (ncol(new21)>2){
    rows_to_keep <- which(rowSums(new21[, -1] != new21[, 1] & 
                                    new21[, -1] != "?", na.rm = TRUE) > 0)
    output<-output[rows_to_keep,] 
    vcf<-vcf[rows_to_keep,]
  } else {
    rows_to_keep <- which(new21[,2]!= new21[, 1] & 
                            new21[, 2] != "?")
    output<-as.data.frame(output[rows_to_keep,])
    vcf<-vcf[rows_to_keep,]
  }

  Reference<-as.character(vcf[,4])
  Position<-vcf[,2]
  forcsv<-cbind(Position,Reference,output)
  names(forcsv)<-c("Position",namesfasta)
  write.csv(forcsv,file=paste0(outputfile,"_SNPs.csv"),row.names=FALSE)
  output_fast<-t(output)
  forfasta<-as.list(apply(output_fast, 1, paste, collapse=""))
  write.fasta(forfasta,names,nbchar=60,file.out=paste0(outputfile,"_SNPs.fasta"),open="w")
}

# Define command-line options
option_list <- list(
  make_option(c("-i", "--inputfile"), type="character", default=NULL, help="Input VCF file", metavar="character"),
  make_option(c("-o", "--outputfile"), type="character", default="output", help="Prefix for output files", metavar="character"),
  make_option(c("--indelfile"), type="character", default=NULL, help="Name of file containing INDELs file if separate to SNPs", metavar="character"),
  make_option(c("-c", "--no.Cores"), type="integer", default=1, help="Number of CPU cores to use", metavar="integer"),
  make_option(c("--samples2remove"), type="character", default=NULL, help="Samples to remove from multisample inputs", metavar="character"),
  make_option(c("--samples2include"), type="character", default=NULL, help="Samples to include from multisample inputs", metavar="character"),
  make_option(c("--filter"), type="logical", default=TRUE, help="Use the 'FILTER' column in VCF file to filter SNPs", metavar="logical"),
  make_option(c("--processIndel"), type="logical", default=FALSE, help="Process INDEL variants", metavar="logical"),
  make_option(c("--DP_low"), type="integer", default=5, help="Minimum read depth to consider call", metavar="integer"),
  make_option(c("--lowqual"), type="integer", default=20, help="Minimum variant quality to consider call", metavar="integer"),
  make_option(c("--hetProp"), type="numeric", default=0.9, help="Proportion allele frequency at hSNPs to assign call", metavar="numeric"),
  make_option(c("--hetasN"), type="logical", default=TRUE, help="Code hSNPs as 'N' or by FASTA nucleic acid code", metavar="logical"),
  make_option(c("--misPercent"), type="integer", default=90, help="Percentage of sites missing across samples to remove variant", metavar="integer"),
  make_option(c("--repeatfile"), type="character", default=NULL, help=".csv file with start and stop coordinates for regions to remove variants", metavar="character"),
  make_option(c("--disINDEL"), type="integer", default=NULL, help="Remove SNPs within int distance of INDELs", metavar="integer"),
  make_option(c("--MixInfect"), type="logical", default=TRUE, help="Run MixInfect to test for mixed infections", metavar="logical"),
  make_option(c("--excludeMix"), type="logical", default=FALSE, help="Remove samples with a high likelihood of mixed infection", metavar="logical")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$inputfile)) {
  print_help(opt_parser)
  stop("Input file must be provided.", call.=FALSE)
}

# Run the function with parsed options
vcfProcess(opt$inputfile, opt$outputfile, opt$indelfile, opt$no.Cores, opt$samples2remove, opt$samples2include, 
           opt$filter, opt$processIndel, opt$DP_low, opt$lowqual, opt$hetProp, opt$hetasN, opt$misPercent, 
           opt$repeatfile, opt$disINDEL, opt$MixInfect, opt$excludeMix)

                           
