#' Filter VCF files and convert to FASTA and CSV files
#' @param inputfiles Single or multisample SNP or SNP/INDEL VCF file
#' @param outputfile Prefix for output files
#' @param indelfiles Name of file containing INDELs file if separate to SNPs
#' @param no.Cores Number of CPU cores to use - if >1, script will run some parallelization 
#' @param samples2remove Samples to remove from multisample inputs, can be string or .txt file in single column (default=NULL)
#' @param samples2include Samples to include from multisample inputs, can be string or .txt file in single column (default=NULL (keep all samples))
#' @param processIndel Process INDEL variants (requires indelProcess_vcfProcess.R)
#' @param DP_low Minimum read depth to consider call (default=5)
#' @param lowqual Minimum variant quality to consider call (defauly=20)
#' @param hetProp Proportion allele frequency at hSNPs to assign call (default=0.9)
#' @param hetasN Code hSNPs as 'N' (TRUE) or by FASTA nucleic acid code (FALSE)
#' @param misPercent Percentage of sites missing across samples to remove variant (default=90)
#' @param repeatfile .csv file with start and stop coordinates for regions to mask
#' @param MixRepeatFile .csv file with start and stop coordinates for regions to mask when running MixInfect (can be NULL if repeatfile specified, otherwise will generate warning)
#' @param disINDEL Remove SNPs within int distance of INDELs 
#' @param MixInfect2 If TRUE, run MixInfect2 to test for mixed infections, if excludeMix = TRUE, will remove samples with a high likelihood of mixed infection (requires MixInfect2_vcfProcess.R)
#' @param LowCov Lowest read depth at mixed sites for MixInfect2 (default = 10)
#' @return filtered .vcf and .csv variant files, optional INDEL and mixed infection files
#' @export

if (!require("stringr")) { install.packages("stringr") }
if (!require("seqinr")) { install.packages("seqinr") }
if (!require("optparse")) { install.packages("optparse") }
if (!require("foreach")) { install.packages("foreach") }
if (!require("doMC")) { install.packages("doMC") }

library(seqinr)
library(stringr)
library(optparse)
library(foreach)
library(doMC)


vcfProcess <- function(inputfiles, outputfile = "output", indelfiles = NULL, no.Cores = 1, samples2remove = NULL, samples2include = NULL, filter = TRUE,
                       processIndel = FALSE, DP_low = 5, lowqual = 20, hetProp = 0.9, hetasN = TRUE, misPercent = 90,
                       repeatfile = NULL, MixRepeatFile = NULL, disINDEL = NULL, MixInfect2 = TRUE, LowCov = 10, excludeMix = FALSE) {
  
  source("~/Documents/Scripts/vcfProcess_functions.R")
  source("~/Documents/Scripts/MixInfect2_vcfProcess.R")
  source("~/Documents/Scripts/indelProcess_vcfProcess.R")
  
  header_info <- get_vcf_header(inputfiles[1])
  header <- header_info$header
  names <- header_info$names
  head_start <- names[1:9]
  format <- which(names == 'FORMAT')
  names <- names[10:length(names)]
  names_all<-names
  filterCol <- which(head_start == "FILTER")
  
  if (is.null(MixRepeatFile)){
    if (!is.null(repeatfile)){
      MixRepeatFile <- repeatfile
    } else {
      sprintf("Warning! No masking file specified for MixInfect2, number of mixed infections likely overestimated. Consider re-running with a masking file")
    }
  }
  
  for (vf in 1:length(inputfiles)) {
    vcf_n <- read_vcf(inputfiles[vf])
    if (!is.null(indelfiles)) {
      indelvcf_n <- read_vcf(indelfiles[vf])
      indel_header_info <- get_vcf_header(indelfiles[vf])
      indel_header <- indel_header_info$header
      vcf_n <- rbind(vcf_n, indelvcf_n)
    }
    
    output_vf<-paste0(outputfile, "_", vf)
    
    sample_info <- remove_or_include_samples(vcf_n, names, samples2remove, samples2include)
    vcf_n <- sample_info$vcf
    names_final <- sample_info$names
    
    vcf_n <- filter_variants(vcf_n, output_vf, filterCol, filter,head_start,names_final)
    vcf_n <- remove_low_quality_variants(vcf_n, output_vf, lowqual, head_start, names_final)
    vcf_n <- remove_spanning_deletions(vcf_n, output_vf, head_start, names_final)
    
    indel_info_n <- process_indels(vcf_n, head_start, names_final)
    vcf_n <- indel_info_n$vcf
    indels_n <- indel_info_n$indels
    indelPos_n <- indel_info_n$indelPos
    
    allele_info_n <- assign_alleles(vcf_n, names_final, head_start, no.Cores, hetProp, hetasN, DP_low)
    output_n <- allele_info_n$output
    read <- allele_info_n$read
    
    output_n <- mark_low_read_positions(output_n, read, DP_low)
    
    invariant_info <- remove_invariant_sites(output_n, vcf_n)
    output_n <- invariant_info$output
    vcf_n <- invariant_info$vcf
    
    repeat_info <- remove_variants_in_repeat_regions(output_n, vcf_n, repeatfile, output_vf, names_final)
    output_n <- repeat_info$output
    vcf_n <- repeat_info$vcf
    
    indel_info <- remove_snps_within_distance_of_indels(output_n, vcf_n, indelPos_n, disINDEL, output_vf, names_final)
    output_n <- indel_info$output
    vcf_n <- indel_info$vcf
    
    if (exists("vcf")){
      vcf <- rbind(vcf, vcf_n)
      indels <- rbind(indels, indels_n)
      output<-rbind(output,output_n)
    } else {
      vcf <- vcf_n
      indels <- indels_n
      output<-output_n
    }
  }
  print("finished allele calling")
  
  mixed_infection_info <- if (MixInfect2) {
    MixInfect2_vcfProcess(vcf, names_final, output, hetProp, outputfile, MixRepeatFile, format, excludeMix, no.Cores, LowCov)
  } else {
    list(vcf = vcf, names = names_final, output = output)
  }
  if (excludeMix){
    vcf <- mixed_infection_info$vcf
    names_final <- mixed_infection_info$names
    output <- mixed_infection_info$output
  }
  missing_data_info <- remove_snps_with_high_missing_data(output, vcf, misPercent)
  output <- missing_data_info$output
  vcf <- missing_data_info$vcf
  
  write_output_files(output, vcf, outputfile, names_final, head_start, header)
  
  if (processIndel && nrow(indels) != 0) {
    if (is.null(indelfiles)) {
      indel_header <- header
    }
    indelfile <- indelProcess(indels, indel_header, names_final, format, hetProp, DP_low, outputfile, repeatfile, misPercent)
  } else {
    print("No indels")
  }
}

option_list <- list(
  make_option(c("-i","--inputfiles"), type = "character", action = "store", help = "Input VCF files, comma-separated", metavar = "character"),
  make_option(c("-o","--outputfile"), type = "character", default = "output", help = "Prefix for output files", metavar = "character"),
  make_option(c("--indelfiles"), type = "character", action = "store", help = "Names of files containing INDELs if separate to SNPs, comma-separated", metavar = "character"),
  make_option(c("-c", "--no.Cores"), type = "integer", default = 1, help = "Number of CPU cores to use", metavar = "integer"),
  make_option(c("--samples2remove"), type = "character", default = NULL, help = "Samples to remove from multisample inputs", metavar = "character"),
  make_option(c("--samples2include"), type = "character", default = NULL, help = "Samples to include from multisample inputs", metavar = "character"),
  make_option(c("--filter"), type = "logical", default = TRUE, help = "Use the 'FILTER' column in VCF file to filter SNPs", metavar = "logical"),
  make_option(c("--processIndel"), type = "logical", default = FALSE, help = "Process INDEL variants", metavar = "logical"),
  make_option(c("--DP_low"), type = "integer", default = 5, help = "Minimum read depth to consider call", metavar = "integer"),
  make_option(c("--lowqual"), type = "integer", default = 20, help = "Minimum variant quality to consider call", metavar = "integer"),
  make_option(c("--hetProp"), type = "numeric", default = 0.9, help = "Proportion allele frequency at hSNPs to assign call", metavar = "numeric"),
  make_option(c("--hetasN"), type = "logical", default = TRUE, help = "Code hSNPs as 'N' or by FASTA nucleic acid code", metavar = "logical"),
  make_option(c("--misPercent"), type = "integer", default = 90, help = "Percentage of sites missing across samples to remove variant", metavar = "integer"),
  make_option(c("--repeatfile"), type = "character", default = NULL, help = ".csv file with start and stop coordinates for regions to remove variants", metavar = "character"),
  make_option(c("--MixRepeatFile"), type = "character", default = NULL, help = ".csv file with start and stop coordinates for regions to mask when running MixInfect (can be NULL if repeatfile specified, otherwise will generate warning)", metavar = "character"),
  make_option(c("--disINDEL"), type = "integer", default = NULL, help = "Remove SNPs within int distance of INDELs", metavar = "integer"),
  make_option(c("--MixInfect2"), type = "logical", default = TRUE, help = "Run MixInfect2 to test for mixed infections", metavar = "logical"),
  make_option(c("--excludeMix"), type = "logical", default = FALSE, help = "Remove samples with a high likelihood of mixed infection", metavar = "logical")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Split comma-separated inputs into lists
if (!is.null(opt$inputfiles)) {
  opt$inputfiles <- strsplit(opt$inputfiles, ",")[[1]]
}
if (!is.null(opt$indelfiles)) {
  opt$indelfiles <- strsplit(opt$indelfiles, ",")[[1]]
}


# Check if input files are provided
if (is.null(opt$inputfiles) || length(opt$inputfiles) == 0) {
  print_help(opt_parser)
  stop("At least one input file must be provided.", call. = FALSE)
}


# Run the function with parsed options
vcfProcess(opt$inputfiles, opt$outputfile, opt$indelfiles, opt$no.Cores, opt$samples2remove, opt$samples2include,
           opt$filter, opt$processIndel, opt$DP_low, opt$lowqual, opt$hetProp, opt$hetasN, opt$misPercent,
           opt$repeatfile, opt$MixRepeatFile, opt$disINDEL, opt$MixInfect2, opt$excludeMix)

