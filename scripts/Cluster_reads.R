#
# Copyright 2022 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@iit.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

#load BioStrings package
suppressMessages(library(Biostrings))
#load stats package
suppressMessages(library(stats))
#load factoextra package
suppressMessages(library(factoextra))

Cluster_reads <- function(fastq_file, first_allele_preliminary, IQR_outliers_coef_precl, IQR_outliers_coef, min_clipped_len, sd_noise_score, num_alleles) {
  IQR_outliers_coef <- as.numeric(IQR_outliers_coef)
  IQR_outliers_coef_precl <- as.numeric(IQR_outliers_coef_precl)
  min_clipped_len <- as.numeric(min_clipped_len)
  sd_noise_score <- as.numeric(sd_noise_score)
  sample_name <- gsub(pattern = "_trimmed", replacement = "", x = basename(gsub(pattern = "\\.fastq", replacement = "", x = fastq_file)))
  sample_dir <- gsub(pattern = "_trimmed", replacement = "", x = gsub(pattern = "inSilicoPCR", replacement = "readsAssignment", x = gsub(pattern = "\\.fastq", replacement = "", x = fastq_file)))
  dir.create(sample_dir)
  logfile <- paste0(dirname(dirname(sample_dir)), "/logfile.txt")
  num_reads_sample <- as.double(system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk seq -A ", fastq_file, " | grep \"^>\" | wc -l"), intern=TRUE))
  DNAStringSet_obj_preliminary <- readDNAStringSet(first_allele_preliminary, "fasta")
  reference_length <- width(DNAStringSet_obj_preliminary)
  #map all reads to the preliminary Allele in order to identify reads coming from the two alleles
  sam_file_reads_to_first_allele <- paste0(sample_dir, "/", sample_name, "_reads_to_allele_preliminary.sam")
  bam_file_reads_to_first_allele <- paste0(sample_dir, "/", sample_name, "_reads_to_allele_preliminary.bam")
  print(sam_file_reads_to_first_allele)
  system(paste0("/opt/conda/envs/CharONT_env/bin/minimap2 -ax map-ont ", first_allele_preliminary, " ", fastq_file, " | /opt/conda/envs/CharONT_env/bin/samtools view -h -F 2308 -o " , sam_file_reads_to_first_allele))
  system(paste0("/opt/conda/envs/CharONT_env/bin/samtools view -hSb ", sam_file_reads_to_first_allele, " | /opt/conda/envs/CharONT_env/bin/samtools sort -o ", bam_file_reads_to_first_allele))
  sam_reads_to_first_allele_tmp <- read.table(file = sam_file_reads_to_first_allele, fill = TRUE, comment.char = "@", stringsAsFactors = FALSE, sep = "\t", quote = "")
  ind_rm <- grep(pattern = "NM", x =  sam_reads_to_first_allele_tmp[, 1])
  if (length(ind_rm) > 0) {
    sam_reads_to_first_allele <- sam_reads_to_first_allele_tmp[-ind_rm, ]
  } else {
    sam_reads_to_first_allele <- sam_reads_to_first_allele_tmp
  }
  #extract reads names and cigar strings
  reads_names <- sam_reads_to_first_allele[, 1]
  start_mapping_coord <- as.numeric(sam_reads_to_first_allele[, 4])
  cigar_strings <- sam_reads_to_first_allele[, 6]
  score <- matrix(0, nrow = length(cigar_strings), ncol = 2)
  del_block_lengths <- c()
  max_del_block_lengths <- c()
  ins_block_lengths <- c()
  max_ins_block_lengths <- c()
  max_clipping_block_lengths <- c()
  max_clipping_block_lengths_signed <- c()
  
  #cycle over reads and extract for each read the longest portion non matching the reference
  for (k in 1:length(cigar_strings)) {
    cigar_string <- cigar_strings[k]
    nonref_symbols_coords <- gregexpr(pattern= "[DIS]", text = cigar_string)[[1]]
    ref_symbols_coords <- gregexpr(pattern= "M", text = cigar_string)[[1]]
    num_starting_coords <- gregexpr(pattern = "[0-9]+", text = cigar_string)[[1]]
    nonref_starting_coords <- gregexpr(pattern = "[0-9]+[DIS]+", text = cigar_string)[[1]]
    match_starting_coords <- gregexpr(pattern = "[0-9]+M+", text = cigar_string)[[1]]
    del_starting_coords <- gregexpr(pattern= "[0-9]+D", text = cigar_string)[[1]]
    ins_starting_coords <- gregexpr(pattern= "[0-9]+I", text = cigar_string)[[1]]
    clipping_starting_coords <- gregexpr(pattern = "[0-9]+[S]+", text = cigar_string)[[1]]
    mapped_length_curr_read <- 0
    mapped_length_curr_read_tmp <- c()
    for (n in 1:length(ref_symbols_coords)) {
      mapped_length_curr <- as.numeric(substr(x = cigar_string, start = (match_starting_coords[n]), stop = (match_starting_coords[n] + attr(match_starting_coords, "match.length")[n] - 2)))
      mapped_length_curr_read_tmp <- c(mapped_length_curr_read_tmp, mapped_length_curr)
    }
    mapped_length_curr_read <- sum(mapped_length_curr_read_tmp)
    del_block_lengths_curr_read <- c()
    ins_block_lengths_curr_read <- c()
    clipping_block_lengths_curr_read <- c()
    clipping_block_lengths_curr_read_signed <- c()
    #if there are non matching bases
    if (nonref_starting_coords[1] != -1) {
      #if there are deletions
      if (del_starting_coords[1] != -1) {
        #extract lengths of deletions
        for (j in 1:length(del_starting_coords)) {
          del_block_curr <- as.numeric(substr(x = cigar_string, start = (del_starting_coords[j]), stop = (del_starting_coords[j] + attr(del_starting_coords, "match.length")[j] - 2)))
          del_block_lengths_curr_read <- c(del_block_lengths_curr_read, del_block_curr)
        }
        del_block_lengths <- c(del_block_lengths, sum(del_block_lengths_curr_read))
        max_del_block_lengths <- c(max_del_block_lengths, max(del_block_lengths_curr_read))
        mapped_length_curr_read <- mapped_length_curr_read + sum(del_block_lengths_curr_read)
      } else {
        del_block_lengths <- c(del_block_lengths, 0)
        max_del_block_lengths <- c(max_del_block_lengths, 0)
      }
      #if there are insertions
      if (ins_starting_coords[1] != -1) {
        #extract lengths of insertions
        for (m in 1:length(ins_starting_coords)) {
          ins_block_curr <- as.numeric(substr(x = cigar_string, start = (ins_starting_coords[m]), stop = (ins_starting_coords[m] + attr(ins_starting_coords, "match.length")[m] - 2)))
          ins_block_lengths_curr_read <- c(ins_block_lengths_curr_read, ins_block_curr)
        }
        ins_block_lengths <- c(ins_block_lengths, sum(ins_block_lengths_curr_read))
        max_ins_block_lengths <- c(max_ins_block_lengths, max(ins_block_lengths_curr_read))
        mapped_length_curr_read <- mapped_length_curr_read + sum(ins_block_lengths_curr_read)
      } else {
        ins_block_lengths <- c(ins_block_lengths, 0)
        max_ins_block_lengths <- c(max_ins_block_lengths, 0)
      }
      #if there are soft-clipped bases
      if (clipping_starting_coords[1] != -1) {
        #extract length of soft-clipped bases if they are longer than min_clipped_len
        for (l in 1:length(clipping_starting_coords)) {
          #if the soft-clipping is at the 5' of the read, subtract the starting coordinate
          if (clipping_starting_coords[l] == 1) {
            if (as.numeric(substr(x = cigar_string, start = (clipping_starting_coords[l]), stop = (clipping_starting_coords[l] + attr(clipping_starting_coords, "match.length")[l] - 2))) > min_clipped_len) {
              clipping_block_curr_signed <- as.numeric(substr(x = cigar_string, start = (clipping_starting_coords[l]), stop = (clipping_starting_coords[l] + attr(clipping_starting_coords, "match.length")[l] - 2))) - start_mapping_coord[k]
              clipping_block_curr <- abs(clipping_block_curr_signed)
              clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, clipping_block_curr)
              clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, clipping_block_curr_signed)
            } else {
              clipping_block_curr_signed <- 0
              clipping_block_curr <- 0
              clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, 0)
              clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, 0)
            }
            #if the soft-clipping is at the 3' of the read, subtract the reference length, the length of the mapping and the starting coordinate
          } else {
            if (as.numeric(substr(x = cigar_string, start = (clipping_starting_coords[l]), stop = (clipping_starting_coords[l] + attr(clipping_starting_coords, "match.length")[l] - 2))) > min_clipped_len) {
              clipping_block_curr_signed <- as.numeric(substr(x = cigar_string, start = (clipping_starting_coords[l]), stop = (clipping_starting_coords[l] + attr(clipping_starting_coords, "match.length")[l] - 2))) - (reference_length - mapped_length_curr_read  - start_mapping_coord[k])
              clipping_block_curr <- abs(clipping_block_curr_signed)
              clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, clipping_block_curr)
              clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, clipping_block_curr_signed)
            } else {
              clipping_block_curr_signed <- 0
              clipping_block_curr <- 0
              clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, 0)
              clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, 0)
            }
          }
          if (clipping_block_curr < min_clipped_len) {
            clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, 0)
            clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, 0)
            
          } else {
            clipping_block_lengths_curr_read <- c(clipping_block_lengths_curr_read, clipping_block_curr)
            clipping_block_lengths_curr_read_signed <- c(clipping_block_lengths_curr_read_signed, clipping_block_curr_signed)
          }
        }
        max_clipping_block_lengths <- c(max_clipping_block_lengths, max(abs(clipping_block_lengths_curr_read)))
        max_clipping_block_lengths_signed <- c(max_clipping_block_lengths_signed, clipping_block_lengths_curr_read_signed[which(abs(clipping_block_lengths_curr_read_signed) == max(abs(clipping_block_lengths_curr_read_signed)))][1])
      } else {
        max_clipping_block_lengths <- c(max_clipping_block_lengths, 0)
        max_clipping_block_lengths_signed <- c(max_clipping_block_lengths_signed, 0)
      }
    } else {
      del_block_lengths_curr_read <- 0
      max_del_block_lengths <- c(max_del_block_lengths, 0)
      ins_block_lengths_curr_read <- 0
      max_ins_block_lengths <- c(max_ins_block_lengths, 0)
      clipping_block_lengths_curr_read <- 0
      max_clipping_block_lengths <- c(max_clipping_block_lengths, 0)
      max_clipping_block_lengths_signed <- c(max_clipping_block_lengths_signed, 0)
    }
  }
  #longest indel model
  #if allele #1 is the one with the shortest repeat soft-clipped portions are probably associated with insertions
  if (max(max_clipping_block_lengths_signed) > 0) {
    ylabplot <- "Max INS or Clipped INS (bp)"
  } else {
    ylabplot <- "Max INS (bp)"
  }
  #if allele #1 is the one with the longest repeat, soft-clipped portions are probably associated with deletions
  if (min(max_clipping_block_lengths_signed) < 0) {
    xlabplot <- "Max DEL or Clipped DEL (bp)"
  } else {
    xlabplot <- "Max DEL (bp)"
  }
  #cycle over reads
  for (k in 1:length(cigar_strings)) {
    #check if read has a soft-clipped portion associated with insertion bigger than insertion
    if (max_clipping_block_lengths_signed[k] > 0) {
      score[k, 2] <- max(max_clipping_block_lengths_signed[k], max_ins_block_lengths[k])
    } else {
      score[k, 2] <- max_ins_block_lengths[k]
    }
    #check if read has a soft-clipped portion associated with deletion bigger than deletion
    if (max_clipping_block_lengths_signed[k] < 0) {
      score[k, 1] <- c(max(-max_clipping_block_lengths_signed[k], max_del_block_lengths[k]))
    } else {
      score[k, 1] <- max_del_block_lengths[k]
    }
  }
  #add some random noise to avoid IQR = 0 with low coverage
  score <- score + rnorm(n = nrow(score)*2, mean = 0, sd = sd_noise_score)
  #check there are at least 2 different score values if clustering has to be performed  
  preclustering_outliers_score_dels <- boxplot.stats(score[, 1], coef = IQR_outliers_coef_precl)$out
  preclustering_outliers_score_ins <- boxplot.stats(score[, 2], coef = IQR_outliers_coef_precl)$out
  ind_preclustering_outliers <- which(score[, 1] %in% preclustering_outliers_score_dels | score[, 2] %in% preclustering_outliers_score_ins)
  num_preclustering_outliers <- length(ind_preclustering_outliers)
  if (num_preclustering_outliers > 0) {
    preclustering_outliers_score <- matrix(score[ind_preclustering_outliers, ], ncol = 2)
    preclustering_score_no_outliers <- matrix(score[-ind_preclustering_outliers, ], ncol = 2)
  } else {
    preclustering_outliers_score <- c()
    preclustering_score_no_outliers <- matrix(score, ncol = 2)
  }
  #skip clustering and assign all reads to one allele if studying haploid chromosome or if there are not 2 different maximum non-matching lengths across all reads
  if (num_alleles == "1" || nrow(unique(preclustering_score_no_outliers)) < 2) {
    if (nrow(unique(preclustering_score_no_outliers)) < 2) {
      cat(text = paste0("WARNING: not possible to identify multiple alleles for sample ", sample_name, ", haploid analysis is performed"), sep = "\n")
      cat(text = paste0("WARNING: not possible to identify multiple alleles for sample ", sample_name, ", haploid analysis is performed"), file = logfile, sep = "\n", append = TRUE)
    }
    first_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_allele_1.fastq")
    first_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_allele_1.fasta")
    first_allele_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_allele_1.txt")
    first_allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_allele_1.fastq")
    first_allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_allele_1.fasta")
    #remove outliers which may be associated with somatic mutations
    ind_outliers <- ind_preclustering_outliers
    num_outliers <- length(ind_outliers) 
    if (num_outliers > 0) {
      score_no_outliers <- score[-ind_outliers, ]
      reads_names_no_outliers <- reads_names[-ind_outliers]
    } else {
      score_no_outliers <- score
      reads_names_no_outliers <- reads_names
    }
    num_outliers_reference <- num_preclustering_outliers
    num_outliers_alternative <- 0
    ind_outliers_reference <- ind_preclustering_outliers
    ind_outliers_alternative <- c()
    outliers_reference_score <- preclustering_outliers_score
    outliers_alternative_score <- c()
    cluster_alternative_index <- c()
    cluster_alternative_index_no_outliers <- c()
    write.table(x=reads_names_no_outliers, quote = FALSE, file = first_allele_reads_names, row.names = FALSE, col.names = FALSE)
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fastq_file, " ", first_allele_reads_names, " > ", first_allele_reads_fq))
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk seq -A ", first_allele_reads_fq, " > ", first_allele_reads_fa))
    allele_reads_fq_all <- first_allele_reads_fq
    #cluster reads into num_alleles groups (ref/alt alleles) based on length of the longest portion non matching the reference (first allele)
  } else {
    #do clustering
    if (num_alleles == "auto") {
      max_num_cl <- 10
      opt <- fviz_nbclust(as.matrix(preclustering_score_no_outliers), kmeans, method = "silhouette", k.max = min(length(unique(preclustering_score_no_outliers)), max_num_cl)) + theme_classic()
      plot(opt)
      ggsave(paste0(sample_dir, "/", sample_name, "_reads_scores_silhouette.png"))
      tol <- 1
      ncl <- max(as.numeric(opt$data$clusters[which(opt$data$y >= tol*max(opt$data$y))]))
      cat(text = paste0("Optimal number of clusters based on Silhouette coefficient: ", ncl), sep = "\n")
      cat(text = paste0("Optimal number of clusters based on Silhouette coefficient: ", ncl), file = logfile, sep = "\n", append = TRUE)
      clusters <- kmeans(preclustering_score_no_outliers, ncl, iter.max = 1000, nstart = 5)
      num_alleles <- ncl
    } else {
      num_alleles <- as.numeric(num_alleles)
      clusters <- kmeans(preclustering_score_no_outliers, num_alleles, iter.max = 1000, nstart = 5)
    }
    #find index and score of reference reads
    cluster_reference_id <- unique(which(rowMeans(clusters$centers) == min(rowMeans(clusters$centers))))
    cluster_reference_index_tmp <- which(clusters$cluster == cluster_reference_id)
    #find index and score of alternative alleles
    cluster_alternative_id <- (1:num_alleles)[-cluster_reference_id]
    cluster_alternative_index_tmp <- lapply(cluster_alternative_id, function (x) which(clusters$cluster == x))
    names(cluster_alternative_index_tmp) <- paste0("Allele_", cluster_alternative_id)
    ref_id <- paste0("Allele_", cluster_reference_id)
    alt_id <- paste0("Allele_", cluster_alternative_id)

    #if there are any preclustering outliers
    if (num_preclustering_outliers > 0) {
      #for each preclustering outlier, calculate distance from the clusters' centers
      diff <- matrix(data = 0, nrow = num_preclustering_outliers, ncol = num_alleles)
      for (l in 1:num_alleles) {
        diff[, l] <- apply(t(abs(apply(matrix(preclustering_outliers_score, ncol = 2), 1, function(x) x - clusters$centers[l, ]))), 1, function(y) sqrt(sum(y^2)))
      }
      #determine if preclustering outliers should be assigned to reference or alternative clusters
      preclustering_outliers_cluster_id <- apply(diff, 1, FUN=which.min)
      #associate preclustering outliers' scores to alleles
      preclustering_outliers_score_split <- lapply(split(preclustering_outliers_score, preclustering_outliers_cluster_id), function(x) matrix(data = x, ncol = 2, byrow = FALSE))
      names(preclustering_outliers_score_split) <- paste0("Allele_", names(preclustering_outliers_score_split))
      preclustering_outliers_score_reference <- preclustering_outliers_score_split[ref_id]
      preclustering_outliers_score_alternative <- preclustering_outliers_score_split[names(preclustering_outliers_score_split) %in% alt_id]
      
      #if any, assign preclustering outliers to reference cluster
      if (length(unlist(preclustering_outliers_score_reference)) > 0) {
        cluster_reference_score <- matrix(rbind(preclustering_outliers_score_reference[[1]], preclustering_score_no_outliers[cluster_reference_index_tmp, ]), ncol = 2)
      } else {
        cluster_reference_score <- matrix(preclustering_score_no_outliers[cluster_reference_index_tmp, ], ncol = 2)
      }
      
      #assign preclustering outliers to alternative clusters
      cluster_alternative_score <- list()
      for (l in alt_id) {
        if (length(preclustering_outliers_score_alternative[[l]]) > 0) {
          cluster_alternative_score[[l]] <- matrix(rbind(preclustering_outliers_score_alternative[[l]], preclustering_score_no_outliers[cluster_alternative_index_tmp[[l]], ]), ncol = 2)
        } else {
          cluster_alternative_score[[l]] <- matrix(preclustering_score_no_outliers[cluster_alternative_index_tmp[[l]], ], ncol = 2)
        }
      }
      #if there are not any preclustering outliers  
    } else {
      cluster_reference_score <- matrix(preclustering_score_no_outliers[cluster_reference_index_tmp, ], ncol = 2)
      cluster_alternative_score <- list()
      for (l in names(cluster_alternative_index_tmp)) {
        cluster_alternative_score[[l]] <- matrix(preclustering_score_no_outliers[cluster_alternative_index_tmp[[l]], ], ncol = 2)
      }
    }
    #identify reads from the reference allele
    cluster_reference_index <- c()
    for (k in 1:length(cluster_reference_score[, 1])) {
      ind <- which(apply(score, 1, function(x) all(x == matrix(cluster_reference_score[k, ], ncol = 2))))
      cluster_reference_index <- unique(c(cluster_reference_index, ind))
    }
    #identify reads from the alternative alleles
    cluster_alternative_index <- list()
    #for each alternative allele
    for (l in 1:length(cluster_alternative_score)) {
      cai_curr <- c()
      #identify reads assigned to alternative allele l
      for (k in 1:length(cluster_alternative_score[[l]][, 1])) {
        ind <- which(apply(score, 1, function(x) all(x == matrix(cluster_alternative_score[[l]][k, ], ncol = 2))))
        cai_curr <- unique(c(cai_curr, ind))
      }
      cluster_alternative_index[[l]] <- cai_curr
    }
    #remove outliers from reference cluster which may be associated with somatic mutations
    outliers_reference_score_dels <- boxplot.stats(cluster_reference_score[, 1], coef = IQR_outliers_coef)$out
    outliers_reference_score_ins <- boxplot.stats(cluster_reference_score[, 2], coef = IQR_outliers_coef)$out
    ind_outliers_reference <- intersect(which(score[, 1] %in% outliers_reference_score_dels | score[, 2] %in% outliers_reference_score_ins), cluster_reference_index)
    outliers_reference_score <- score[ind_outliers_reference, ]
    score_reference_no_outliers <- matrix(score[setdiff(cluster_reference_index, ind_outliers_reference), ], ncol = 2)
    num_outliers_reference <- nrow(cluster_reference_score) - nrow(score_reference_no_outliers)
    if (num_outliers_reference > 0) {
      cluster_reference_index_no_outliers <- setdiff(cluster_reference_index, ind_outliers_reference)
      cluster_reference_readnames <- reads_names[cluster_reference_index_no_outliers]
    } else {
      cluster_reference_index_no_outliers <- cluster_reference_index
      cluster_reference_readnames <- reads_names[cluster_reference_index]
    }
    #identify outliers from alternative clusters which may be associated with somatic mutations
    outliers_alternative_score_dels <- lapply(cluster_alternative_score, function(x) boxplot.stats(x[, 1], coef = IQR_outliers_coef)$out)
    outliers_alternative_score_ins <- lapply(cluster_alternative_score, function(x) boxplot.stats(x[, 2], coef = IQR_outliers_coef)$out)
    ind_outliers_alternative <- list()
    for (j in 1:(num_alleles - 1)) {
      ind_outliers_alternative[[j]] <- intersect(which(score[, 1] %in% outliers_alternative_score_dels[[j]] | score[, 2] %in% outliers_alternative_score_ins[[j]]), cluster_alternative_index[[j]])
    }
    #remove alternative alleles' outliers from scores 
    outliers_alternative_score <- score[unlist(ind_outliers_alternative), ]
    score_alternative_no_outliers <- matrix(score[setdiff(unlist(cluster_alternative_index), unlist(ind_outliers_alternative)), ], ncol = 2)
    num_outliers_alternative <- sum(unlist(lapply(cluster_alternative_score, function(x) nrow(x))))- nrow(score_alternative_no_outliers)
    if (num_outliers_alternative > 0) {
      cluster_alternative_index_no_outliers <- mapply(setdiff, cluster_alternative_index, ind_outliers_alternative, SIMPLIFY = FALSE)
      cluster_alternative_readnames <- lapply(cluster_alternative_index_no_outliers, function(x) reads_names[x])
    } else {
      cluster_alternative_index_no_outliers <- cluster_alternative_index
      cluster_alternative_readnames <- lapply(cluster_alternative_index, function(x) reads_names[x])
    }
    #sort clusters' centers for median and identify reads names associated with each allele
    median_cluster_reference_nonreflen <- clusters$centers[cluster_reference_id, ]
    median_cluster_alternative_nonreflen <- clusters$centers[cluster_alternative_id, ]
    clusters$median <- rbind(median_cluster_reference_nonreflen, median_cluster_alternative_nonreflen, deparse.level = 0)
    cluster_readnames <- c(list(cluster_reference_readnames), cluster_alternative_readnames)
    
    allele_reads_fq_all <- c()
    for (k in 1:num_alleles) {
      allele_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_allele_", k, ".txt")
      allele_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_allele_", k, ".fastq")
      allele_reads_fq_all <- c(allele_reads_fq_all, allele_reads_fq)
      allele_reads_fa <- paste0(sample_dir, "/", sample_name, "_reads_allele_", k, ".fasta")
      write.table(x = cluster_readnames[[k]], quote = FALSE, file = allele_reads_names, row.names = FALSE, col.names = FALSE)
      system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fastq_file, " ", allele_reads_names, " > ", allele_reads_fq))
      system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk seq -A ", allele_reads_fq, " > ", allele_reads_fa))
    }
    outliers_readnames <- reads_names[unlist(c(ind_outliers_reference, ind_outliers_alternative))]
    outliers_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_outliers.txt")
    write.table(x = outliers_readnames, quote = FALSE, file = outliers_reads_names, row.names = FALSE, col.names = FALSE)
    outliers_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_outliers.fastq")
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fastq_file, " ", outliers_reads_names, " > ", outliers_reads_fq))
  }
  num_outliers <- num_outliers_reference + num_outliers_alternative
  ind_outliers <- c(ind_outliers_reference, unlist(ind_outliers_alternative))
  allelic_ratio_outliers <- num_outliers/num_reads_sample
  allelic_ratio_perc_outliers <- allelic_ratio_outliers*100
  #col_list <- colours()[sample(x = 1:length(colours()), size = 10, replace = FALSE)]
  #col_list <- rainbow(n = 10)
  col_list <- c('#3cb44b', '#4363d8', '#ffe119', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  
  if (num_outliers > 0) {
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_smoothScatter.png"))
    smoothScatter(score[-ind_outliers, ], xlab = xlabplot, ylab = ylabplot, main = "Reads scores")
    dev.off()
    png(paste0(sample_dir, "/", sample_name, "_reads_scores.png"))
    plot(matrix(score[-ind_outliers, ], ncol = 2), xlab = xlabplot, ylab = ylabplot, main = "Reads scores", col = "black", type = "p", pch = 19, cex = 2, xlim = c(0, max(score[, 1])*1.5),  ylim = c(0, max(score[, 2])*1.5))
    points(matrix(score[ind_outliers, ], ncol = 2), col = "red2", type = "p", pch = 15, cex = 2)
    if (num_alleles > 1) {
      for (i in 1:(num_alleles - 1)) {
        points(matrix(score[cluster_alternative_index_no_outliers[[i]], ], ncol = 2), col = col_list[i], type = "p", pch = 19, cex = 2)
      }
      legend(x = "topright", legend = c(ref_id, alt_id, "Outliers"), col = c("black", col_list[1:(num_alleles - 1)], "red2"), cex = 1.5, pch = 19)
    } else {
      legend(x = "topright", legend = c("Allele_1", "Outliers"), col = c("black", "red2"), cex = 1.5, pch = 19)
    }
    dev.off()
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_logscale.png"))
    plot(matrix(score[-ind_outliers, ], ncol = 2), xlab = xlabplot, ylab = ylabplot, log = "xy", main = "Reads scores", col = "black", type = "p", pch = 19, cex = 2)
    points(matrix(score[ind_outliers, ], ncol = 2), col = "red2", type = "p", pch = 15, cex = 2)
    if (num_alleles > 1) {
      for (i in 1:(num_alleles - 1)) {
        points(matrix(score[cluster_alternative_index_no_outliers[[i]], ], ncol = 2), col = col_list[i], type = "p", pch = 19, cex = 2)
      }
      legend(x = "bottomright", legend = c(ref_id, alt_id, "Outliers"), col = c("black", col_list[1:(num_alleles - 1)], "red2"), cex = 1.5, pch = 19)
    } else {
      legend(x = "bottomright", legend = c("Allele_1", "Outliers"), col = c("black", "red2"), cex = 1.5, pch = 19)
    }
    dev.off()
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_no_outliers.png"))
    plot(matrix(score[-ind_outliers, ], ncol = 2), xlab = xlabplot, ylab = ylabplot, main = "Reads scores", col = "black", type = "p", pch = 19, cex = 2, xlim = c(0, max(score[-ind_outliers, 1])*1.5), ylim = c(0, max(score[-ind_outliers, 2])*1.5))
    if (num_alleles > 1) {
      for (i in 1:(num_alleles - 1)) {
        points(matrix(score[cluster_alternative_index_no_outliers[[i]], ], ncol = 2), col = col_list[i], type = "p", pch = 19, cex = 2)
      }
      legend(x = "topright", legend = c(ref_id, alt_id), col = c("black", col_list[1:(num_alleles - 1)]), cex = 1.5, pch = c(19, 19))
    } else {
      legend(x = "topright", legend = "Allele_1", col = "black", cex = 1.5, pch = 19)
    }
    dev.off()
    cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads (", sprintf("%.2f", allelic_ratio_perc_outliers), "%), possibly associated with somatic mutations, have been discarded"), sep = "\n")
    cat(text = paste0("Sample ", sample_name, ": ", sprintf("%d", num_outliers), " reads (", sprintf("%.2f", allelic_ratio_perc_outliers), "%), possibly associated with somatic mutations, have been discarded"),  file = logfile, sep = "\n", append = TRUE)
    outliers_readnames <- reads_names[ind_outliers]
    outliers_reads_names <- paste0(sample_dir, "/", sample_name, "_reads_names_outliers.txt")
    write.table(x=outliers_readnames, quote = FALSE, file = outliers_reads_names, row.names = FALSE, col.names = FALSE)
    outliers_reads_fq <- paste0(sample_dir, "/", sample_name, "_reads_outliers.fastq")
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fastq_file, " ", outliers_reads_names, " > ", outliers_reads_fq))
    reads_names_no_outliers <- reads_names[-ind_outliers]
  } else {
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_smoothScatter.png"))
    smoothScatter(score, xlab = xlabplot, ylab = ylabplot, main = "Reads scores")
    dev.off()
    png(paste0(sample_dir, "/", sample_name, "_reads_scores.png"))
    plot(matrix(score, ncol = 2), xlab = xlabplot, ylab = ylabplot, main = "Reads scores", col = "black", type = "p", pch = 19, cex = 2, xlim = c(0, max(score[, 1])*1.5),  ylim = c(0, max(score[, 2])*1.5))
    if (num_alleles > 1) {
      for (i in 1:(num_alleles - 1)) {
        points(matrix(score[cluster_alternative_index_no_outliers[[i]], ], ncol = 2), col = col_list[i], type = "p", pch = 19, cex = 2)
      }
      legend(x = "topright", legend = c(ref_id, alt_id), col = c("black", col_list[1:(num_alleles - 1)]), cex = 1.5, pch = 19)
    } else {
      legend(x = "topright", legend = "Allele_1", col = "black", cex = 1.5, pch = 19)
    }
    dev.off()
    png(paste0(sample_dir, "/", sample_name, "_reads_scores_logscale.png"))
    plot(matrix(score, ncol = 2), xlab = xlabplot, ylab = ylabplot, log = "xy", main = "Reads scores", col = "black", type = "p", pch = 19, cex = 2)
    if (num_alleles > 1) {
      for (i in 1:(num_alleles - 1)) {
        points(matrix(score[cluster_alternative_index_no_outliers[[i]], ], ncol = 2), col = col_list[i], type = "p", pch = 19, cex = 2)
      }
      legend(x = "bottomright", legend = c(ref_id, alt_id), col = c("black", col_list[1:(num_alleles - 1)]), cex = 1.5, pch = 19)
    } else {
      legend(x = "bottomright", legend = "Allele_1", col = "black", cex = 1.5, pch = 19)
    }
    dev.off()
    reads_names_no_outliers <- reads_names
  }
  cat(text = paste0("Reads files for sample ", sample_name, ": ", allele_reads_fq_all), sep = "\n")
  return(allele_reads_fq_all)
}

Cluster_reads(fastq_file, first_allele_preliminary, IQR_outliers_coef_precl, IQR_outliers_coef, min_clipped_len, sd_noise_score, num_alleles)
