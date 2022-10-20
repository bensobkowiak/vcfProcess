#' Identify samples with mixed infection and remove from vcfProcess
#' @param vcf Processed VCF file
#' @param names Sample names
#' @param output Genotype file
#' @param LowCov Minimum read depth at site to accurately call allele frequency
#' @param BICoutput Option to produce output file with probability of each sample being mixed
#' @return vcf, names and output genotype files with mixed infections removed
#' @export

MixInfect = function(vcf,names,output,hetProp,outputfile,format,excludeMix,LowCov=10){
  
  if (!require(mclust)){
    install.packages("mclust")
    library("mclust")
  }
  
  vcf<-vcf
  names<-names
  output<-output
  hetProp<-hetProp
  outputfile<-outputfile
  format<-format
  
  #### Determine AD field and create matrix of AD and GT
  AD<-which(unlist(str_split(vcf[1,format], ":"))=='AD')
  GT<-which(unlist(str_split(vcf[1,format], ":"))=='GT')
  
  #### Make new matrices of separared GT and AD fields
  vcf_mat<-vcf[,10:ncol(vcf)]
  GT_mat<-matrix(0,nrow(vcf_mat),ncol(vcf_mat))
  AD_mat<-matrix(0,nrow(vcf_mat),ncol(vcf_mat))
  for (i in 1:nrow(vcf_mat)){
    for (j in 1:ncol(vcf_mat)){
      s<-unlist(str_split(vcf_mat[i,j], ":"))
      GT_mat[i,j]<-s[GT]
      AD_mat[i,j]<-s[AD]
    } 
  }
  GT_mat[which(output=="?")]<-"?"
  
  ### Create output file
  
  outfile<-as.data.frame(matrix(NA,nrow=length(names),ncol=7))
  colnames(outfile)<-c("Sample name","Mix or Non-mix","Mixed SNPs","Total SNPs","Proportion het/total SNPs","Number of strains in mix","Major strain proportion")
  outfile[,1]<-names
  
  #### No of het SNPs and total and proportions
  mixes<-matrix(0,ncol=ncol(GT_mat),nrow=4)
  for (i in 1:ncol(GT_mat)){
    mixes[1,i]<-length(which(GT_mat[,i]=="0/1"))
  }
  for (i in 1:ncol(GT_mat)){
    mixes[2,i]<-length(which(GT_mat[,i]=="1/1"))
  }
  for (i in 1:ncol(GT_mat)){
    mixes[3,i]<-mixes[1,i]+mixes[2,i]
  }
  for (i in 1:ncol(GT_mat)){
    mixes[4,i]<-(mixes[1,i]/mixes[3,i])*100
  }
  
  outfile[,3]<-mixes[1,]
  outfile[,4]<-mixes[3,]
  outfile[,5]<-mixes[4,]
  outfile[,2]<-'Non-mix'
  outfile[,6]<-1
  
  #################### ESTIMATE PROPORTIONS OF MIXED SAMPLES (up to 3 mixes in sample) #######################
  
  mixgeno<-as.data.frame(GT_mat[,which(outfile[,5]>0.5 & outfile[,3]>10)])
  mixAD<-as.data.frame(AD_mat[,which(outfile[,5]>0.5 & outfile[,3]>10)])
  mixnames<-names[which(outfile[,5]>0.5 & outfile[,3]>10)]
  
  if (length(mixnames)>0){
    
    #### Remove any SNPs that do not have a hetero call at any sample
    new22<-matrix(0,nrow(mixgeno),ncol(mixgeno))
    for (i in 1:nrow(mixgeno)){
      for (j in 1:ncol(mixgeno)){
        if (mixgeno[i,j]=="0/1"){
          new22[i,j]='FALSE'}
        else  {new22[i,j]='TRUE'
        }
      }
    }
    d<- lapply(1:nrow(new22), function(i){
      all(as.logical(new22[i,]))
    }
    )
    g<-do.call(rbind,d)
    mix_GT<-cbind(vcf[which(g=='FALSE'),2],mixgeno[which(g=='FALSE'),])
    mix_AD<-cbind(vcf[which(g=='FALSE'),2],mixAD[which(g=='FALSE'),])
    
    ## Create proportion matrix
    
    output_prop=matrix(0,nrow(mix_GT),ncol(mix_GT))
    output_prop[,1]<-mix_GT[,1]
    colnames(output_prop)<-c("Position",mixnames)
    for (i in 1:nrow(output_prop)){
      for (j in 2:ncol(output_prop)){
        if (mix_GT[i,j]=="0/1"){
          props<-as.numeric(unlist(strsplit(mix_AD[i,j],split = ",")))
          if (length(props)==2 && sum(props)>=LowCov){
            output_prop[i,j]<-props[2]/sum(props)
          } else if (length(props)>2 &&sum(max(props),max(props[props!=max(props)]))>=LowCov){
            output_prop[i,j]<-max(props)/sum(max(props),max(props[props!=max(props)]))
          } else{
            output_prop[i,j]<-0
          }
        } else {
          output_prop[i,j]<-0
        }
      }
    }
    
    ####### Run Guassian Mclust and idenitify 2 or 3 mixes
    
    BICvalues<-as.data.frame(matrix(NA,ncol=4,nrow = ncol(output_prop)-1))
    BICvalues[,1]<-mixnames
    colnames(BICvalues)<-c("Sample identifier","G=2","G=4","G=6")
    
    for (i in 1:nrow(BICvalues)){
      b<-as.numeric(output_prop[,(i+1)])
      b<-b[which(b!=0)]
      b<-c(b,1-b)
      if (length(b)>0){
        a<-mclustBIC(b,G=c(2,4,6),verbose = F)[,2]
        if (length(a)==3){
          BICvalues[i,2:4]<-a
        } else { BICvalues[i,2:4]<-c(a,rep(NA,3-length(a)))
        }
      }
      ind<-which(is.element(outfile$`Sample name`,BICvalues$`Sample identifier`[i]))
      if (BICvalues[i,2]>=20 && is.na(BICvalues[i,2])==FALSE){
        d<-Mclust(b,G=2,verbose = F)
        outfile[ind,2]<-'Mix'
        outfile[ind,6]<-2
        outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][1]
      } else if (BICvalues[i,4]>=20 && is.na(BICvalues[i,4])==FALSE){
        d<-Mclust(b,G=6,verbose = F)
        outfile[ind,2]<-'Mix'
        outfile[ind,6]<-3
        means<-d$parameters$mean[order(d$parameters$mean,decreasing = T)]
        if (sum(means[5:6])<0.5){
          outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][3]
        } else {
          outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][4]
        }
      }
    }
    
    ### remove mixed infections 
    
    mixed<-which(outfile$`Mix or Non-mix`=="Mix")
    mixoutfile<-outfile[mixed,]
    write.csv(mixoutfile,paste0(outputfile,"_MixSampleSummary.csv"),row.names = F)
    if (excludeMix){
    newoutput<-output[,-mixed]
    newnames<-names[-mixed]
    newvcf<-cbind(vcf[,1:9],vcf_mat[,-mixed])
    mixed_SNPs<-cbind(vcf[,2],vcf[,4],output[,mixed])
    colnames(mixed_SNPs)<-c("Position","Reference",names[mixed])
    write.csv(mixed_SNPs,paste0(outputfile,"_MixedSNPs.csv"),row.names = F)
    
    ## Remove invariant sites in remaining samples
    ref<-as.character(newvcf[,4])
    mat<-cbind(ref,newoutput)
    new22<-matrix(0,nrow(mat),ncol(mat))
    for (i in 1:nrow(mat)){
      for (j in 2:ncol(mat)){
        if (mat[i,j]==mat[i,1]){
          new22[i,j]='FALSE'}
        else  {
          new22[i,j]='TRUE'
        }
      }
    }
    new22<-new22[,2:ncol(new22)]
    d<- lapply(1:nrow(new22), function(i){
      all(as.logical(new22[i,]))
    }
    )
    g<-do.call(rbind,d)
    newoutput<-newoutput[which(g=='FALSE'),] 
    removedSNPs<-newoutput[which(g=='TRUE'),] 
    newvcf<-newvcf[which(g=='FALSE'),]
    
    # return new values
    vcf=newvcf
    names=newnames
    output=newoutput
    }
  } else {
    print("No Mixed infection")
  }
  return(list(vcf=vcf,names=names,output=output))
}

