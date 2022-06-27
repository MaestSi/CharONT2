#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/CharONT2
========================================================================================
 MaestSi/CharONT2 analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/CharONT2
----------------------------------------------------------------------------------------
*/
def helpMessage() {
        log.info"""
    Usage:
    nextflow -c CharONT2.conf run CharONT2.nf -profile docker
    Mandatory argument:
    -profile                                                              Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the CharONT2.conf file
    --fastq_files                                                         Path to fastq files, use wildcards to select multiple samples
    --results_dir                                                         Path to a folder where to store results
    --num_alleles                                                         num_alleles represents the number of expected alleles
    --PCRThr                                                              PCRThr is the identity threshold for in-silico PCR in case inSilicoPCR=true
    --primerSeqOne                                                        primerSeqOne is a primer sequence used for in-silico PCR in case inSilicoPCR=true
    --primerSeqTwo                                                        primerSeqTwo is a primer sequence used for in-silico PCR in case inSilicoPCR=true
    --scripts_dir                                                         scripts_dir is the directory containing all scripts
    --target_reads_consensus                                              target_reads_consensus defines the maximum number of reads used for consensus calling
    --target_reads_polishing                                              target_reads_polishing defines the maximum number of reads used for consensus polishing
    --max_reads_preliminary                                               max_reads_preliminary defines the maximum number of reads used for preliminary clustering and consensus calling
    --clustering_id_threshold                                             identity threshold for clustering preliminary allele assembly
    --plurality                                                           MAFFT plurality value: minimum fraction of aligned reads supporting a basis for including it in the preliminary consensus
    --min_maf                                                             minimum minor allele frequency; if less than min_maf*100% of reads are assigned to Allele //2, the sample is assumed homozygous
    --IQR_outliers_coef_precl                                             label as candidate outliers reads with score > 3rd_QR + IQR_outliers_coef_precl*IQR or score < 1st_QR - IQR_outliers_coef_precl*IQR
    --IQR_outliers_coef                                                   label as outliers reads with score > 3rd_QR + IQR_outliers_coef*IQR or score < 1st_QR - IQR_outliers_coef*IQR; IQR is computed within each cluster
    --fast_alignment_flag                                                 set fast_alignment_flag=1 if you want to perform fast multiple sequence alignment; otherwise set fast_alignment_flag=0
    --min_clipped_len                                                     minimum number of soft-clipped bases to be considered as a DEL or INS
    --sd_noise_score                                                      sd_noise_score is the standard deviation of gaussian-distributed noise with zero mean added to score
    --primers_length                                                      primers_length defines how many bases are trimmed from consensus sequences
    --medaka_model                                                        medaka model for consensus polishing
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of sample names, conditions, and FAST5s path.
Channel
    .fromPath(params.fastq_files)
    .map {tuple( it.name.split('\\.')[0], it )}
    .set{inputFiles_inSilicoPCR}

    // in-silico PCR
process inSilicoPCR {
  input:
    tuple val(sample), val(fastq) from inputFiles_inSilicoPCR
  output:
    val(sample) into inSilicoPCR_preliminaryConsensusCalling
  script:
  if(params.inSilicoPCR)
  """
    mkdir -p ${params.results_dir}/inSilicoPCR
    sample=\$(echo \$(basename \$(realpath ${fastq})) | sed \'s/\\.fa.*//\' | sed \'s/\\.fq.*//\')
    sam_file_one=\$sample"_in_silico_pcr_one.sam"
    sam_file_two=\$sample"_in_silico_pcr_two.sam"
    trimmed_reads_fq=${params.results_dir}"/inSilicoPCR/"\$sample"_trimmed.fastq"
    export PATH=\$PATH:/opt/conda/envs/CharONT_env/bin/
    mem=\$(echo ${task.memory} | sed \'s/ //\' | sed \'s/B//\' | sed \'s/G/g/\' | sed \'s/M/m/\')
    /opt/conda/envs/CharONT_env/bin/msa.sh in=${fastq} out=\$sam_file_one literal=${params.primerSeqOne} qin=33 cutoff=${params.PCRThr} -Xmx\$mem
    /opt/conda/envs/CharONT_env/bin/msa.sh in=${fastq} out=\$sam_file_two literal=${params.primerSeqTwo} qin=33 cutoff=${params.PCRThr} -Xmx\$mem
    /opt/conda/envs/CharONT_env/bin/cutprimers.sh in=${fastq} out=\$trimmed_reads_fq sam1=\$sam_file_one sam2=\$sam_file_two qin=33 fake=f include=t fixjunk -Xmx\$mem
    ln -s \$trimmed_reads_fq ./trimmed.fastq
  
  """
  else
  """
    mkdir -p ${params.results_dir}/inSilicoPCR
    reads_full=\$(realpath ${fastq})
    sample=\$(echo \$(basename \$(realpath ${fastq})) | sed \'s/\\.fa.*//\' | sed \'s/\\.fq.*//\')
    trimmed_reads_fq=${params.results_dir}"/inSilicoPCR/"\$sample"_trimmed.fastq"
    cp ${fastq} \$trimmed_reads_fq

  """
}

// Preliminary consensus sequence generation
process preliminaryConsensusCalling {
    input:
        val(sample) from inSilicoPCR_preliminaryConsensusCalling

    output:
        val(sample) into preliminaryConsensusCalling_readsAssignment

    script:
    if(params.preliminaryConsensusCalling)
    """
        mkdir -p ${params.results_dir}/preliminaryConsensusCalling
        export PATH=\$PATH:/opt/conda/envs/CharONT_env/bin/
        /opt/conda/envs/CharONT_env/bin/Rscript ${params.scripts_dir}/Obtain_preliminary_consensus.R fastq_file=${params.results_dir}/inSilicoPCR/${sample}_trimmed.fastq TRC=${params.target_reads_consensus} max_num_reads_clustering=${params.max_reads_preliminary} THR=${params.clustering_id_threshold} PLUR=${params.plurality} num_threads=${task.cpus}
    """
    else
    """
       echo "Skipped."
    """
  }


process readsAssignment {
    input:
         val(sample) from preliminaryConsensusCalling_readsAssignment
    
    output:
        val(sample) into readsAssignment_draftConsensusCalling

    script:
    if(params.readsAssignment)
    """
        mkdir -p ${params.results_dir}/readsAssignment
        export PATH=\$PATH:/opt/conda/envs/CharONT_env/bin/
        /opt/conda/envs/CharONT_env/bin/Rscript ${params.scripts_dir}/Cluster_reads.R fastq_file=${params.results_dir}/inSilicoPCR/${sample}_trimmed.fastq first_allele_preliminary=${params.results_dir}/preliminaryConsensusCalling/${sample}_trimmed/${sample}_trimmed_preliminary_allele.fasta IQR_outliers_coef_precl=${params.IQR_outliers_coef_precl} IQR_outliers_coef=${params.IQR_outliers_coef} min_clipped_len=${params.min_clipped_len} sd_noise_score=${params.sd_noise_score} num_alleles=${params.num_alleles}
    """
    else
    """
        echo "Skipped."
    """
}

process draftConsensusCalling {
    input:
        val(sample) from readsAssignment_draftConsensusCalling

    output:
        val(sample) into draftConsensusCalling_consensusPolishing        

    script:
    if(params.draftConsensusCalling)
    """
        mkdir -p ${params.results_dir}/draftConsensusCalling
        export PATH=\$PATH:/opt/conda/envs/CharONT_env/bin/

        fastq_files=\$(find ${params.results_dir}/readsAssignment/${sample}_trimmed | grep "\\d.*\\.fastq" | grep -v "outliers")

        for f in \$fastq_files; do
           allele_number=\$(echo \$(basename \$f) | sed \'s/.*_allele_//\' | sed \'s/\\.fastq//\')
          /opt/conda/envs/CharONT_env/bin/Rscript ${params.scripts_dir}/Obtain_draft_consensus.R fastq_file=\$f allele_num=\$allele_number min_maf=${params.min_maf} TRC=${params.target_reads_consensus} PLUR=${params.plurality} num_threads=${task.cpus} fast_alignment_flag=${params.fast_alignment_flag}
        done
        
    """
    else
    """
        echo "Skipped."
    """
}

process consensusPolishing {
    input:
        val(sample) from draftConsensusCalling_consensusPolishing

    output:
        val(sample) into consensusPolishing_trfAnnotate
    
    script:
    if(params.consensusPolishing)
    """
        mkdir -p ${params.results_dir}/consensusPolishing
        export PATH=\$PATH:/opt/conda/envs/CharONT_env/bin/

        draft_consensus_files=\$(find ${params.results_dir}/draftConsensusCalling/${sample}_trimmed | grep ${sample}"_trimmed_draft_allele.*\\.fasta" | grep -v "tmp")

        for f in \$draft_consensus_files; do
           allele_number=\$(echo \$(basename \$f) | sed \'s/.*_allele_//\' | sed \'s/\\.fasta//\')
           fastq_file=\$(find ${params.results_dir}/readsAssignment/${sample}_trimmed | grep "\\d.*\\.fastq" | grep "allele_"\$allele_number)
           /opt/conda/envs/CharONT_env/bin/Rscript ${params.scripts_dir}/Polish_consensus.R draft_consensus=\$f fastq_file=\$fastq_file allele_num=\$allele_number TRP=${params.target_reads_polishing} num_threads=${task.cpus} primers_length=${params.primers_length} medaka_model=${params.medaka_model}
        done
        
    """
    else
    """
        mkdir -p ${params.results_dir}/consensusPolishing
        draft_consensus_files=\$(find ${params.results_dir}/draftConsensusCalling/${sample}_trimmed | grep ${sample}"_trimmed_draft_allele.*\\.fasta" | grep -v "tmp")
        for f in \$draft_consensus_files; do
            polished_consensus_file=\$(echo \$f | sed \'s/draft/polished/\' | sed \'s/draftConsensusCalling/consensusPolishing/\')
            /opt/conda/envs/CharONT_env/bin/seqtk trimfq \$f -b ${params.primers_length} -e ${params.primers_length}  > \$polished_consensus_file
        done
    """
}

process trfAnnotate {
    input:
          val(sample) from consensusPolishing_trfAnnotate

    output:
        

    script:
    if(params.trfAnnotate)
    """
        mkdir -p ${params.results_dir}/trfAnnotate/
        polished_consensus_files=\$(find ${params.results_dir}/consensusPolishing/${sample}_trimmed | grep ${sample}"_trimmed_allele.*\\.fasta" | grep -v "untrimmed")

        for f in \$polished_consensus_files; do
            sample_name=\$(echo \$(basename \$f) | sed \'s/_allele.*//\')
            mkdir -p ${params.results_dir}/trfAnnotate/\$sample_name
            cp \$f ${params.results_dir}/trfAnnotate/\$sample_name
            cd ${params.results_dir}/trfAnnotate/\$sample_name
            /opt/conda/envs/CharONT_env/bin/trf \$f 2 7 7 80 10 50 500
            /opt/conda/envs/CharONT_env/bin/trf \$f 2 3 5 80 10 50 500
            /opt/conda/envs/CharONT_env/bin/trf \$f 2 500 500 80 10 50 500
        done
    """
    else
    """
        mkdir -p ${params.results_dir}/trfAnnotate/
        for f in \$polished_consensus_files; do
            sample_name=\$(echo \$(basename \$f) | sed \'s/_allele.*//\')
            mkdir -p ${params.results_dir}/trfAnnotate/\$sample_name
            cp \$f ${params.results_dir}/trfAnnotate/\$sample_name
        done
    """
}

    
