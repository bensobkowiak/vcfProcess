#' Filter VCF files and convert to FASTA and CSV files
#' @param inputfile Single or multisample VCF file
#' @param outputfile Prefix for output files
#' @param samples2remove If is not NULL, removes named samples from multisample inputs
#' @param indelProcess Process INDEL variants (requires indelProcess_vcfProcess.R)
#' @param DP_low Minimum read depth to consider call
#' @param lowqual Minimum variant quality to consider call
#' @param hetProp Proportion allele frequency at hSNPs to assign call
#' @param misPercent Percentage of sites missing across samples to remove variant
#' @param repeatfile .csv file with start and stop coordinates for regions to remove variants
#' @param excludeMix If true, will remove samples with a high likelihood of mixed infection (requires MixInfect_vcfProcess.R)
#' @return filtered .vcf and .csv variant files, optional INDEL file and mixed infection files
#' @export

vcfProcess = function(inputfile,outputfile,samples2remove=NULL,indelProcess=FALSE,DP_low=5,lowqual=20,hetProp=0.9,misPercent=80,repeatfile=NULL,excludeMix=FALSE){
  
  ######################################   VCF FILE POST PROCESSING  ##################################
  
  if (!require(stringr)){
    install.packages("stringr",repos = "http://cran.us.r-project.org")
    library(stringr)
  }
  if (!require(seqinr)){
    install.packages("seqinr",repos = "http://cran.us.r-project.org")
    library(seqinr)
  }
  options(stringsAsFactors = F)
  
  ######READ IN VCF FILE, REPEAT REGION FILE AND RETAIN HEADER############################################################
  
  vcf<-read.table(inputfile)
  header_input<-as.matrix(read.table(inputfile,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",header_input)==TRUE)
  header<-as.data.frame(header_input[1:end_head-1])
  names<-unlist(strsplit(header_input[end_head],"\t"))
  head_start<-names[1:9]
  format<-which(names=='FORMAT')
  names<-names[10:length(names)]
  
  ###### Remove samples
  
  if (is.null(samples2remove)==FALSE){
    remove<-samples2remove
    vcf<-vcf[,-(which(is.element(names,remove))+9)]
    names<-names[-which(is.element(names,remove))]
  }
  
  
  ##################REMOVE UNWANTED SNPS ############################################################
  
  #### DETERMINE ORDER OF FORMAT FIELDS ####
  
  st<-unlist(str_split(vcf[1,format], ":"))
  GT<-which(st=='GT')
  DP<-which(st=='DP')
  PL<-which(st=='PL')
  AD<-which(st=='AD')
  
  # if no AD field, hetProp set to 1
  if (length(AD)==0){
    AD<-which(st=='DP')
    hetProp<-1
  }
  
  ######### REMOVE LOW QUAL VARIANTS
  
  lowqual_mat<-vcf[which((as.integer(vcf[,6]) < lowqual) ==T),]
  colnames(lowqual_mat)<-c(head_start,names)
  vcf<-vcf[which((as.integer(vcf[,6]) < lowqual) == F),]
  if (nrow(lowqual_mat)!=0){
    write.csv(lowqual_mat,file=paste0(outputfile,"_lowqualsnps.csv"),row.names = F)
  } else {print("No low quality SNPs")}
  
  ######## MAKE INDEL FILE (if applicable) 
  
  ind<- lapply(1:nrow(vcf), function(i){
    length(unlist(strsplit(vcf[i,4],"")))>1 || length(unlist(strsplit(unlist(strsplit(vcf[i,5],","))[1],"")))>1
  }
  )
  indels<-vcf[which(ind==TRUE),]
  colnames(indels)<-c(head_start,names)
  vcf<-vcf[which(ind==FALSE),]
  write.csv(indels,file=paste0(outputfile,"_InDels.csv"),row.names = F)
  if (indelProcess==TRUE){
    if (nrow(indels)!=0){
      indelfile<-indelProcess(indels,names,GT,AD,DP,hetProp,DP_low,outputfile,repeatfile.present,repeatfile)
    } else {print("No indels")
    }
  }
  
  #######SPLIT FIELDS AND MAKE NEW MATRICES ###########################################################
  
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
  
  ####### ALLELE ASSIGNING #######################################################################
  
  output=matrix("N",nrow(genotype),ncol(genotype))
  rownames(output)=rownames(genotype)
  colnames(output)=colnames(genotype)
  mixed_sites<-matrix(c("R","GA","AG","M","CA","AC","W","TA","AT","Y","TC","CT","S","GC","CG","K","TG","GT"),ncol=3,byrow = T)
  
  for (i in 1:nrow(genotype)){
    alt<-unlist(strsplit(genotype[i,2],split = ","))
    output[i,as.numeric(which(unlist(genotype[i,])=="./."))]<-"N"
    output[i,as.numeric(which(unlist(genotype[i,])=="0/0"))]<-genotype[i,1]
    output[i,as.numeric(which(unlist(genotype[i,])=="1/1"))]<-alt[1]
    output[i,as.numeric(which(unlist(genotype[i,])=="2/2"))]<-alt[2]
    output[i,as.numeric(which(unlist(genotype[i,])=="3/3"))]<-alt[3]
    mix01<-which(genotype[i,]=="0/1")
    mix02<-which(genotype[i,]=="0/2")
    mix03<-which(genotype[i,]=="0/3")
    mix12<-which(genotype[i,]=="1/2")
    mix13<-which(genotype[i,]=="1/3")
    mix23<-which(genotype[i,]=="2/3")
    if (length(mix01)>0){
      for (j in 1:length(mix01)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix01[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[1]/sum(AD_vec)>hetProp){
            output[i,mix01[j]]<-genotype[i,1]
          } else if (sum(AD_vec)>0 & AD_vec[2]/sum(AD_vec)>hetProp){
            output[i,mix01[j]]<-alt[1]
          } else {
            output[i,mix01[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[1])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[1]))),1]}
        } else {
          output[i,mix01[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[1])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[1]))),1]}
      }
    }
    if (length(mix02)>0){
      for (j in 1:length(mix02)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix02[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[1]/sum(AD_vec)>hetProp){
            output[i,mix02[j]]<-genotype[i,1]
          } else if (sum(AD_vec)>0 & AD_vec[3]/sum(AD_vec)>hetProp){
            output[i,mix02[j]]<-alt[2]
          } else {
            output[i,mix02[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[2])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[2]))),1]}
        } else {
          output[i,mix02[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[2])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[2]))),1]}
      } 
    }
    if (length(mix03)>0){
      for (j in 1:length(mix03)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix03[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[1]/sum(AD_vec)>hetProp){
            output[i,mix03[j]]<-genotype[i,1]
          } else if (sum(AD_vec)>0 & AD_vec[4]/sum(AD_vec)>hetProp){
            output[i,mix03[j]]<-alt[3]
          } else {
            output[i,mix03[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[3])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[3]))),1]}
        } else {
          output[i,mix03[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(genotype[i,1],alt[3])) | is.element(mixed_sites[,3],paste0(genotype[i,1],alt[3]))),1]}
      }
    }
    if (length(mix12)>0){
      for (j in 1:length(mix12)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix12[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[2]/sum(AD_vec)>hetProp){
            output[i,mix12[j]]<-alt[1]
          } else if (sum(AD_vec)>0 & AD_vec[3]/sum(AD_vec)>hetProp){
            output[i,mix12[j]]<-alt[2]
          } else {
            output[i,mix12[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[2],alt[1])) | is.element(mixed_sites[,3],paste0(alt[2],alt[1]))),1]}
        } else {
          output[i,mix12[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[2],alt[1])) | is.element(mixed_sites[,3],paste0(alt[2],alt[1]))),1]}
      }
    }
    if (length(mix13)>0){
      for (j in 1:length(mix13)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix13[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[2]/sum(AD_vec)>hetProp){
            output[i,mix13[j]]<-alt[1]
          } else if (sum(AD_vec)>0 & AD_vec[4]/sum(AD_vec)>hetProp){
            output[i,mix13[j]]<-alt[3]
          } else {
            output[i,mix13[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[3],alt[1])) | is.element(mixed_sites[,3],paste0(alt[3],alt[1]))),1]}
        } else {
          output[i,mix13[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[3],alt[1])) | is.element(mixed_sites[,3],paste0(alt[3],alt[1]))),1]}
      }
    }
    if (length(mix23)>0){
      for (j in 1:length(mix23)){
        if (hetProp<1){
          AD_vec<-as.numeric(unlist(strsplit(AD_comp[i,mix23[j]],split = ",")))
          if (sum(AD_vec)>0 & AD_vec[3]/sum(AD_vec)>hetProp){
            output[i,mix23[j]]<-alt[2]
          } else if (sum(AD_vec)>0 & AD_vec[4]/sum(AD_vec)>hetProp){
            output[i,mix23[j]]<-alt[3]
          } else {
            output[i,mix23[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[2],alt[3])) | is.element(mixed_sites[,3],paste0(alt[2],alt[3]))),1]}
        } else {
          output[i,mix23[j]]=mixed_sites[which(is.element(mixed_sites[,2],paste0(alt[2],alt[3])) | is.element(mixed_sites[,3],paste0(alt[2],alt[3]))),1]}
      }
    }
  }
  
  output<-output[,3:ncol(output)] 
  
  ######## low read for all mark as 'N' 
  
  read<-as.matrix(read[,3:ncol(read)])
  class(read)<-"numeric"
  snprd<-read<DP_low
  lowread<-which(snprd == 'TRUE') 
  output[lowread]<-'N'
  
  ########### Remove snps that have no variant
  
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
  
  ####### REMOVE REPEAT AND MOBILE ELEMENT REGIONS
  
  if (is.null(repeatfile)==FALSE){
    rep<-as.data.frame(read.csv(repeatfile,header=TRUE))
    pos<-as.numeric(vcf[,2])
    res=vector()
    for (i in 1:nrow(rep)){
      re<-c(rep[i,1]:rep[i,2])
      res<-c(re,y)   
    }
    res<-is.element(pos,res)
    rep<-cbind(vcf[which(res == 'TRUE'),2],vcf[which(res == 'TRUE'),4],output[which(res=='TRUE'),])
    colnames(rep)<-c("Position","Reference",names)
    output<-output[which(res == 'FALSE'),]
    vcf<-vcf[which(res == 'FALSE'),]
    if (nrow(rep)!=0){
      write.csv(rep, file=paste0(outputfile,"_SNPsinRepRegions.csv"),row.names = F)
    } else {print("No repeat or mobile element SNPs")
    }
  }
  
  ### Check for mixed infection, remove mixed samples and reassign mixed calls #####
  
  if (excludeMix==TRUE){
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
  
  ################# writing different file types##############################################
  
  ########## Write VCF file of processed SNPs
  
  vfile<-as.data.frame(apply(vcf, 1, paste, collapse="\t"))
  rownames(vfile)<-NULL
  headrow<-paste0(c(head_start,names),collapse = "\t")
  names(headrow)<-names(header)
  names(vfile)<-names(header)
  outvcf<-rbind(header,headrow,v)
  write.table(outvcf, file= paste0(outputfile,".vcf"),row.names =FALSE,sep="\t",quote=FALSE,col.names=FALSE)
  
  ########### CSV file with SNP genotypes only with reference sequence and SNP position 
  
  colnames(output)<-names
  Reference<-as.character(vcf[,4])
  Position<-vcf[,2]
  forcsv<-cbind(Position,Reference,output)
  write.csv(forcsv,file=paste0(outputfile,"_with_ref.csv"),row.names=FALSE)
  
  
  ########### Fasta file of all isolates (with ref)
  
  namesfasta<-c("Reference",names)
  output_fast<-t(output)
  forfastaref<-rbind(Reference,output_fast)
  forfastaref<-as.list(apply(forfastaref, 1, paste, collapse=""))
  write.fasta(forfastaref,namesfasta,nbchar=60,file.out=paste0(outputfile,"_with_ref.fasta"),open="w")
  
  ########### FASTA of samples only (removing monomorphic sites in samples)
  
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
