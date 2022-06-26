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
    nextflow -c charont.conf run charont.nf --fastq_files = "/path/to/*.fastq" --results_dir = "/path/to/results_dir" 
-profile docker
    Mandatory argument:
    -profile                        Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the charont.conf file
      
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
        
    """
}

process consensusPolishing {
    input:


    output:
        

    script:
    if(params.consensusPolishing)
    """
        
    """
    else
    """
        
    """
}

process trfAnnotate {
    input:


    output:
        

    script:
    if(params.trfAnnotate)
    """
        
    """
    else
    """
        
    """
}

    
