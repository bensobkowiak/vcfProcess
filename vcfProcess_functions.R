## Functions for vcfProcess.R

options(stringsAsFactors = FALSE)
get_vcf_header <- function(inputfile) {
  header_input <- as.matrix(read.table(inputfile, comment.char = " ", sep = "\n"))
  end_head <- which(grepl("#CHROM", header_input) == TRUE)
  header <- as.data.frame(header_input[1:end_head - 1])
  names <- unlist(strsplit(header_input[end_head], "\t"))
  list(header = header, names = names)
}

read_vcf <- function(inputfile) {
  read.table(inputfile)
}

remove_or_include_samples <- function(vcf, names, samples2remove, samples2include) {
  if (!is.null(samples2remove)) {
    remove <- read.table(samples2remove, check.names = FALSE)[, 1]
    vcf <- vcf[, which(!names %in% remove) + 9]
    names <- names[!(names %in% remove)]
  }
  if (!is.null(samples2include)) {
    include <- read.table(samples2include, check.names = FALSE)[, 1]
    vcf <- cbind(vcf[, 1:9], vcf[, which(names %in% include) + 9])
    names <- names[names %in% include]
  }
  list(vcf = vcf, names = names)
}

filter_variants <- function(vcf, outputfile, filterCol, filter, head_start, names) {
  if (filter) {
    failedQC_mat <- vcf[vcf[, filterCol] != "PASS", ]
    colnames(failedQC_mat) <- c(head_start, names)
    vcf <- vcf[vcf[, filterCol] == "PASS", ]
    if (nrow(failedQC_mat) != 0) {
      write.csv(failedQC_mat, file = paste0(outputfile, "_filteredVariants.csv"), row.names = FALSE)
    } else {
      print("No filtered variants")
    }
  }
  vcf
}

remove_low_quality_variants <- function(vcf, outputfile, lowqual, head_start, names) {
  lowqual_mat <- vcf[as.integer(vcf[, 6]) < lowqual, ]
  colnames(lowqual_mat) <- c(head_start, names)
  vcf <- vcf[as.integer(vcf[, 6]) >= lowqual, ]
  if (nrow(lowqual_mat) != 0) {
    write.csv(lowqual_mat, file = paste0(outputfile, "_lowqualVariants.csv"), row.names = FALSE)
  } else {
    print("No low quality variants")
  }
  vcf
}

remove_spanning_deletions <- function(vcf, outputfile, head_start, names) {
  if (length(grep("*", vcf[, 5], fixed = TRUE)) > 1) {
    overlap_mat <- vcf[grep("*", vcf[, 5], fixed = TRUE), ]
    colnames(overlap_mat) <- c(head_start, names)
    vcf <- vcf[-grep("*", vcf[, 5], fixed = TRUE), ]
    if (nrow(overlap_mat) != 0) {
      write.csv(overlap_mat, file = paste0(outputfile, "_overlapVariants.csv"), row.names = FALSE)
    } else {
      print("No overlap variants")
    }
  }
  vcf
}

process_indels <- function(vcf, head_start, names) {
  ind <- lapply(1:nrow(vcf), function(i) {
    length(unlist(strsplit(vcf[i, 4], ""))) > 1 || length(unlist(strsplit(unlist(strsplit(vcf[i, 5], ","))[1], ""))) > 1
  })
  indels <- vcf[which(ind == TRUE), ]
  indelPos <- indels[, 2]
  colnames(indels) <- c(head_start, names)
  vcf <- vcf[which(ind == FALSE), ]
  list(vcf = vcf, indels = indels, indelPos = indelPos)
}

assign_alleles <- function(vcf, names, head_start, no.Cores, hetProp, hetasN, DP_low) {
  genotype <- data.frame(vcf[, 4:5])
  AD_comp <- data.frame(vcf[, 4:5])
  read <- data.frame(vcf[, 4:5])
  for (i in 10:ncol(vcf)) {
    split_fieldsGT <- sapply(strsplit(vcf[, i], split = ":"), function(x) x[1])
    split_fieldsAD <- sapply(strsplit(vcf[, i], split = ":"), function(x) x[2])
    split_fieldsDP <- sapply(strsplit(vcf[, i], split = ":"), function(x) x[3])
    genotype <- cbind(genotype, split_fieldsGT)
    AD_comp <- cbind(AD_comp, split_fieldsAD)
    read <- cbind(read, split_fieldsDP)
  }
  mixed_sites <- matrix(c("R", "GA", "AG", "M", "CA", "AC", "W", "TA", "AT", "Y", "TC", "CT", "S", "GC", "CG", "K", "TG", "GT"), ncol = 3, byrow = TRUE)
  max_length <- max(unlist(lapply(strsplit(genotype[, 2], split = ","), length)))
  
  if (no.Cores > 1) {
    registerDoMC(no.Cores)
    output <- foreach(i = 1:nrow(genotype)) %dopar% {
      df <- data.frame(matrix("N", ncol = ncol(genotype)))
      df[1, as.numeric(which(unlist(genotype[i, ]) == "./."))] <- "N" # Missing
      df[1, as.numeric(which(unlist(genotype[i, ]) == ".|."))] <- "N" # Missing
      df[1, as.numeric(which(unlist(genotype[i, ]) == "0/0"))] <- genotype[i, 1] # Ref call
      df[1, as.numeric(which(unlist(genotype[i, ]) == "0|0"))] <- genotype[i, 1] # Ref call
      for (j in 1:max_length) {
        df[1, as.numeric(which(unlist(genotype[i, ]) == paste0(j, "/", j) | unlist(genotype[i, ]) == paste0(j, "|", j)))] <- unlist(strsplit(genotype[i, 2], ","))[j]
      } # Alt call
      alleles <- c(genotype[i, 1], unlist(strsplit(genotype[i, 2], ",")))
      for (j in 0:(max_length - 1)) {
        for (k in (j + 1):max_length) {
          mixed <- which(genotype[i, ] == paste0(j, "/", k) | genotype[i, ] == paste0(j, "|", k))
          if (length(mixed) > 0) {
            for (mix in 1:length(mixed)) {
              if (hetProp < 1) {
                AD_vec <- as.numeric(unlist(strsplit(AD_comp[i, mixed[mix]], split = ",")))
                if (sum(AD_vec) > 0 && AD_vec[(j + 1)] / sum(AD_vec) > hetProp) {
                  df[1, mixed[mix]] <- alleles[j + 1]
                } else if (sum(AD_vec) > 0 && AD_vec[k + 1] / sum(AD_vec) > hetProp) {
                  df[1, mixed[mix]] <- alleles[k + 1]
                } else if (hetasN == FALSE) {
                  df[1, mixed[mix]] <- mixed_sites[which(is.element(mixed_sites[, 2], paste0(alleles[j + 1], alleles[k + 1])) | is.element(mixed_sites[, 3], paste0(alleles[j + 1], alleles[k + 1]))), 1]
                }
              } else if (hetasN == FALSE) {
                df[1, mixed[mix]] <- mixed_sites[which(is.element(mixed_sites[, 2], paste0(alleles[j + 1], alleles[k + 1])) | is.element(mixed_sites[, 3], paste0(alleles[j + 1], alleles[k + 1]))), 1]
              }
            }
          }
        }
      }
      df
    }
    output <- matrix(unlist(output), nrow = length(output), byrow = TRUE)
    colnames(output) <- colnames(genotype)
  } else {
    output <- matrix("N", nrow(genotype), ncol(genotype))
    rownames(output) <- rownames(genotype)
    colnames(output) <- colnames(genotype)
    for (i in 1:nrow(genotype)) {
      output[i, as.numeric(which(unlist(genotype[i, ]) == "./."))] <- "N" # Missing
      output[i, as.numeric(which(unlist(genotype[i, ]) == ".|."))] <- "N" # Missing
      output[i, as.numeric(which(unlist(genotype[i, ]) == "0/0"))] <- genotype[i, 1] # Ref call
      output[i, as.numeric(which(unlist(genotype[i, ]) == "0|0"))] <- genotype[i, 1] # Ref call
      for (j in 1:max_length) {
        output[i, as.numeric(which(unlist(genotype[i, ]) == paste0(j, "/", j) | unlist(genotype[i, ]) == paste0(j, "|", j)))] <- unlist(strsplit(genotype[i, 2], ","))[j]
      } # Alt call
      alleles <- c(genotype[i, 1], unlist(strsplit(genotype[i, 2], ",")))
      for (j in 0:(max_length - 1)) {
        for (k in (j + 1):max_length) {
          mixed <- which(genotype[i, ] == paste0(j, "/", k) | genotype[i, ] == paste0(j, "|", k))
          if (length(mixed) > 0) {
            for (mix in 1:length(mixed)) {
              if (hetProp < 1) {
                AD_vec <- as.numeric(unlist(strsplit(AD_comp[i, mixed[mix]], split = ",")))
                if (sum(AD_vec) > 0 && AD_vec[(j + 1)] / sum(AD_vec) > hetProp) {
                  output[i, mixed[mix]] <- alleles[j + 1]
                } else if (sum(AD_vec) > 0 && AD_vec[k + 1] / sum(AD_vec) > hetProp) {
                  output[i, mixed[mix]] <- alleles[k + 1]
                } else if (hetasN == FALSE) {
                  output[i, mixed[mix]] <- mixed_sites[which(is.element(mixed_sites[, 2], paste0(alleles[j + 1], alleles[k + 1])) | is.element(mixed_sites[, 3], paste0(alleles[j + 1], alleles[k + 1]))), 1]
                }
              } else if (hetasN == FALSE) {
                output[i, mixed[mix]] <- mixed_sites[which(is.element(mixed_sites[, 2], paste0(alleles[j + 1], alleles[k + 1])) | is.element(mixed_sites[, 3], paste0(alleles[j + 1], alleles[k + 1]))), 1]
              }
            }
          }
        }
      }
    }
  }
  output <- output[, 3:ncol(output)]
  list(output = output, read = read)
}

mark_low_read_positions <- function(output, read, DP_low) {
  read<-as.matrix(read[,3:ncol(read)])
  suppressWarnings(class(read)<-"numeric")
  snprd<-read<DP_low
  output[which(snprd | is.na(snprd))]<-'?'
  output
}

remove_invariant_sites <- function(output, vcf) {
  new21 <- as.matrix(cbind(vcf[, 4], output))
  if (ncol(new21) > 2) {
    rows_to_keep <- which(rowSums(new21[, -1] != new21[, 1] & new21[, -1] != "?", na.rm = TRUE) > 0)
    output <- output[rows_to_keep, ]
    vcf <- vcf[rows_to_keep, ]
  } else {
    rows_to_keep <- which(new21[, 2] != new21[, 1] & new21[, 2] != "?")
    output <- as.data.frame(output[rows_to_keep, ])
    vcf <- vcf[rows_to_keep, ]
  }
  list(output = output, vcf = vcf)
}

remove_variants_in_repeat_regions <- function(output, vcf, repeatfile, outputfile, names) {
  if (!is.null(repeatfile)) {
    rep <- as.data.frame(read.csv(repeatfile, header = TRUE))
    res1 <- vector()
    for (i in 1:nrow(rep)) {
      re <- c(rep[i, 1]:rep[i, 2])
      res1 <- c(res1, re)
    }
    res <- is.element(as.numeric(vcf[, 2]), res1)
    rep <- cbind(vcf[which(res == 'TRUE'), 2], vcf[which(res == 'TRUE'), 4], output[which(res == 'TRUE'), ])
    colnames(rep) <- c("Position", "Reference", names)
    output <- as.data.frame(output[which(res == 'FALSE'), ])
    vcf <- vcf[which(res == 'FALSE'), ]
    if (nrow(rep) != 0) {
      write.csv(rep, file = paste0(outputfile, "_SNPsinRepRegions.csv"), row.names = FALSE)
    } else {
      print("No repeat or mobile element SNPs")
    }
  }
  list(output = output, vcf = vcf)
}

remove_snps_within_distance_of_indels <- function(output, vcf, indelPos, disINDEL, outputfile, names) {
  if (!is.null(disINDEL) & length(indelPos) > 0) {
    res1 <- vector()
    for (i in 1:length(indelPos)) {
      re <- (indelPos[i] - disINDEL):(indelPos[i] + disINDEL)
      res1 <- c(res1, re)
    }
    res <- is.element(as.numeric(vcf[, 2]), res1)
    inds <- cbind(vcf[which(res == 'TRUE'), 2], vcf[which(res == 'TRUE'), 4], output[which(res == 'TRUE'), ])
    colnames(inds) <- c("Position", "Reference", names)
    output <- output[which(res == 'FALSE'), ]
    vcf <- vcf[which(res == 'FALSE'), ]
    if (nrow(inds) != 0) {
      write.csv(inds, file = paste0(outputfile, "_SNPswithin", disINDEL, "bpINDELs.csv"), row.names = FALSE)
    } else {
      print("No SNPS within distance of INDELs")
    }
  }
  list(output = output, vcf = vcf)
}

remove_snps_with_high_missing_data <- function(output, vcf, misPercent) {
  percent_missing <- integer(nrow(output))
  for (i in 1:nrow(output)) {
    percent_missing[i] <- ((length(which(output[i, ] == "?")) + length(which(output[i, ] == "N"))) / ncol(output)) * 100
  }
  output <- as.data.frame(output[which(percent_missing <= misPercent), ])
  vcf <- vcf[which(percent_missing <= misPercent), ]
  list(output = output, vcf = vcf)
}

write_output_files <- function(output, vcf, outputfile, names, head_start, header) {
  vfile <- as.data.frame(apply(vcf, 1, paste, collapse = "\t"))
  rownames(vfile) <- NULL
  headrow <- paste0(c(head_start, names), collapse = "\t")
  names(headrow) <- names(header)
  names(vfile) <- names(header)
  outvcf <- rbind(header, headrow, vfile)
  write.table(outvcf, file = paste0(outputfile, "_SNPs_processed.vcf"), row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
  
  colnames(output) <- names
  Reference <- as.character(vcf[, 4])
  Position <- vcf[, 2]
  write.table(data.frame(SNP=1:length(Position),Position=Position), file = paste0(outputfile, "_SNPsIndex_with_ref.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  forcsv <- cbind(Position, Reference, output)
  write.csv(forcsv, file = paste0(outputfile, "_SNPs_with_ref.csv"), row.names = FALSE)
  
  namesfasta <- c("Reference", names)
  output_fast <- t(output)
  forfastaref <- rbind(Reference, output_fast)
  forfastaref <- as.list(apply(forfastaref, 1, paste, collapse = ""))
  write.fasta(forfastaref, namesfasta, nbchar = 60, file.out = paste0(outputfile, "_SNPs_with_ref.fasta"), open = "w")
  
  new21 <- as.matrix(cbind(vcf[, 4], output))
  if (ncol(new21) > 2) {
    rows_to_keep <- which(rowSums(new21[, -1] != new21[, 1] & new21[, -1] != "?", na.rm = TRUE) > 0)
    output <- output[rows_to_keep, ]
    vcf <- vcf[rows_to_keep, ]
  } else {
    rows_to_keep <- which(new21[, 2] != new21[, 1] & new21[, 2] != "?")
    output <- as.data.frame(output[rows_to_keep, ])
    vcf <- vcf[rows_to_keep, ]
  }
  
  Reference <- as.character(vcf[, 4])
  Position <- vcf[, 2]
  write.table(data.frame(SNP=1:length(Position),Position=Position), file = paste0(outputfile, "_SNPsIndex.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  forcsv <- cbind(Position, Reference, output)
  names(forcsv) <- c("Position", namesfasta)
  write.csv(forcsv, file = paste0(outputfile, "_SNPs.csv"), row.names = FALSE)
  output_fast <- t(output)
  forfasta <- as.list(apply(output_fast, 1, paste, collapse = ""))
  write.fasta(forfasta, names, nbchar = 60, file.out = paste0(outputfile, "_SNPs.fasta"), open = "w")
}
