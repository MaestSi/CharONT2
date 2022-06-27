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

Obtain_preliminary_consensus <- function(fastq_file, TRC, max_num_reads_clustering, THR, PLUR, num_threads) {
  TRC <- as.numeric(TRC)
  THR <- as.numeric(THR)
  PLUR <- as.numeric(PLUR)
  num_threads <- as.numeric(num_threads)
  fasta_file <- gsub(pattern = "\\.fastq", replacement = ".fasta", x = fastq_file)
  system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk seq -A ", fastq_file , " > ", fasta_file))
  sample_name <- gsub(pattern = "_trimmed", replacement = "", x = basename(gsub(pattern = "\\.fasta", replacement = "", x = fasta_file)))
  sample_dir <- gsub(pattern = "_trimmed", replacement = "", x = gsub(pattern = "inSilicoPCR", replacement = "preliminaryConsensusCalling", x = gsub(pattern = "\\.fasta", replacement = "", x = fasta_file)))
  subset_reads_fa <- paste0(sample_dir, "/", sample_name, "_subset_reads.fasta")
  dir.create(sample_dir)
  vsearch_clustering_dir <- paste0(sample_dir, "/clustering_vsearch")
  dir.create(vsearch_clustering_dir)
  #perform vsearch clustering and create a preliminary version for one Allele
  seed <- 1
  system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk sample -s ", seed , " ", fasta_file, " ",  max_num_reads_clustering, " > ", subset_reads_fa))
  ids_mac_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_ids_mac.txt")
  mac_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_mac.fasta")
  mac_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_reads_mac.fastq")
  system(paste0("/opt/conda/envs/CharONT_env/bin/vsearch --cluster_smallmem ", subset_reads_fa, " --usersort --id ", THR, " --iddef 2 --clusterout_sort --fasta_width 0 --maxseqlength 300000000 --strand both --sizeout --consout ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta --clusters ", vsearch_clustering_dir, "/", sample_name, "_cluster"))
  centroid_mac <- system(paste0("head -n1 ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta"), intern = TRUE)
  id_centroid <-  system(paste0("echo ", "\"", centroid_mac, "\"", " | sed 's/centroid=//g' | sed 's/;seqs.*$//g'"), intern = TRUE)
  clusters_vsearch <- list.files(path = vsearch_clustering_dir, pattern = paste0(sample_name, "_cluster"), full.names = TRUE)
  for (j in 1:length(clusters_vsearch)) {
    cluster_match <-  grep(id_centroid, readLines(clusters_vsearch[j]))
    if (length(cluster_match) > 0) {
      mac_file <- clusters_vsearch[j]
      system(paste0("cat ", mac_file, " | grep ", "\"",  "^>", "\"", "  | sed 's/^>//' | sed 's/;.*$//' > ", ids_mac_first_preliminary))
      system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fasta_file, " ", ids_mac_first_preliminary, " > ", mac_fa_first_preliminary))
      system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk subseq ", fastq_file, " ", ids_mac_first_preliminary, " > ", mac_fq_first_preliminary))
      break
    }
  }
  num_reads_sample <- as.double(system(paste0("cat ", fasta_file, " | grep \"^>\" | wc -l"), intern=TRUE))
  num_reads_mac_first_preliminary <- as.double(system(paste0("cat ", mac_fa_first_preliminary, " | grep \"^>\" | wc -l"), intern=TRUE))
  target_reads_consensus <- TRC
  first_allele_preliminary_tmp1 <- paste0(sample_dir, "/", sample_name, "_preliminary_allele_tmp1.fasta")
  first_allele_preliminary_tmp2 <- paste0(sample_dir, "/", sample_name, "_preliminary_allele_tmp2.fasta")
  first_allele_preliminary <- paste0(sample_dir, "/", sample_name, "_preliminary_allele.fasta")
  if (num_reads_mac_first_preliminary < 3) {
    system(paste0("head -n2 ", vsearch_clustering_dir, "/", sample_name, "_consensus.fasta > ", first_allele_preliminary))
  } else {
    if (num_reads_mac_first_preliminary < target_reads_consensus) {
      target_reads_consensus <- num_reads_mac_first_preliminary
    }
    plurality_value <- PLUR*target_reads_consensus
    sequences <- readDNAStringSet(fasta_file, "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_preliminary_allele.fastq")
    draft_reads_fa_first_preliminary <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_consensus, "_reads_preliminary_allele.fasta")
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk sample -s ", seed , " ", mac_fq_first_preliminary, " ",  target_reads_consensus, " > ", draft_reads_fq_first_preliminary))
    system(paste0("/opt/conda/envs/CharONT_env/bin/seqtk seq -A ", draft_reads_fq_first_preliminary, " > ", draft_reads_fa_first_preliminary))
    mfa_file_first_preliminary <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_first_preliminary)
    system(paste0("/opt/conda/envs/CharONT_env/bin/mafft --auto --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_first_preliminary, " > ", mfa_file_first_preliminary))
    system(paste0("/opt/conda/envs/CharONT_env/bin/cons -sequence ", mfa_file_first_preliminary, " -plurality ", PLUR, " -outseq ", first_allele_preliminary_tmp1))
    system(paste0("sed 's/[nN]//g' ", first_allele_preliminary_tmp1, " > ", first_allele_preliminary_tmp2))
    DNAStringSet_obj <- readDNAStringSet(first_allele_preliminary_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    names(DNAStringSet_obj_renamed) <- "Allele_preliminary"
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = first_allele_preliminary, format = "fasta", width = 20000)
    return(list(first_allele_preliminary, num_reads_sample))
  }
}

Obtain_preliminary_consensus(fastq_file, TRC, max_num_reads_clustering, THR, PLUR, num_threads)