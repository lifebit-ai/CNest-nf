#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf --project test --step 4

    Mandatory arguments:
      --project       [string] Name of the project
      --step          [int] step of the workflow to run. One of [1, 2, 3, 4]


    step 1 arguments:
      --ref           [file] Reference FASTA
      --indexb        [file] Index in BED format for fast counting
      --design        [file] A CSV file with header and three columns (name,cram,crai)

    step 2 arguments:
      --index         [file] Index in tab format (index_tab.txt)
      --bindir        [path] Path to the directory of all bin files
      --binlist       [file] A txt file with paths to all bin files (one per row)

    step 3 arguments:
      --index         [file] Index in tab format (index_tab.txt)
      --bindir        [path] Path to the directory of all bin files
      --gender        [file] Gender file from step 2
      --batch_size         [int]  Batch size for references
      --samples       [file] Samples to process

    step 4 arguments:
      --rbindir
      --cordir
      --index
      --gender
      --cov

    Optional arguments:
      --wgs           [int] indicate the memory factor for WGS
      --test          [flag/int] Run the pipeline in a test mode and subset total samples to n samples. (5 default if no number provided)
      --help          [flag] Show help messages

    """.stripIndent()
}

// Show help message

if (params.help) exit 0, helpMessage()

/*
================================================================================
                                Set parameters
================================================================================
*/

if (params.bed) ch_bed = Channel.value(file(params.bed))
if (params.ref) ch_ref = Channel.value(file(params.ref))
if (params.design) {
  Channel.fromPath(params.design)
    .splitCsv(sep: ',', skip: params.skip_lines)
    .map { name, file_path, index_path -> [ name, file(file_path), file(index_path) ] }
    .set { ch_files_sets }
}

// Helper files
if (params.index_tab) ch_index_tab = Channel.value(file(params.index_tab))
if (params.indexb) ch_index_bed = Channel.value(file(params.indexb))
if (params.gender) ch_gender = Channel.value(file(params.gender))
if (params.cov) ch_cov = Channel.value(file(params.cov))

// Test mode
if (params.test == true) {subset_samples = 5} else {subset_samples = params.test}
if (params.test && params.design) ch_files_sets = ch_files_sets.take(subset_samples)
if (params.test && (params.bindir || params.binlist || params.rbindir || params.samples)) ch_sample_names = ch_sample_names.take(subset_samples)

/*
================================================================================
                                Steps
================================================================================
*/
if (params.step =~ 1) {
  ch_bedgz = Channel.value(file(params.bedgz))

  process step_0_bedgz_uncompress {
    tag "${params.project}"
    //echo true

    input:
    file(bedgz) from ch_bedgz

    output:
    file("hg38.1kb.baits.bed") into ch_bed

    when:
    !params.bed

    script:
    if (params.test)
      """
      gzip -cd ${bedgz} | head -1000 > "hg38.1kb.baits.bed"
      """
    else
      """
      gzip -cd ${bedgz} > "hg38.1kb.baits.bed"
      """
  }

  // Step1 create work directory
  process step_1_index_gen {
    tag "${params.project}"
    publishDir "results/", mode: params.mode, pattern: "${params.project}/index*"
    //echo true

    input:
    file(bed) from ch_bed

    output:
    path "${params.project}/index_tab.txt" into ch_index_tab
    path "${params.project}/index.txt" into ch_index
    path "${params.project}/index.bed" into ch_index_bed

    script:
    """
    cnest.py step1 --project ${params.project} --bed ${bed}
    """
  }
}

if (params.step =~ 2) {
  process step_2_bin_gen {
    tag "sample:${name}"
    publishDir "results/", mode: params.mode
    //echo true

    input:
    set val(name), file(file_path), file(index_path) from ch_files_sets
    file("genome.fa") from ch_ref
    path "${params.project}/index.bed" from ch_index_bed

    output:
    path "${params.project}/bin/$name" into ch_bin_sample

    script:
    if (params.test)
      """
      mkdir -p ${params.project}/tmp/ ${params.project}/bin/
      cnest.py --debug step2 --project ${params.project} --sample ${name} --input ${file_path} --fasta genome.fa --fast
      """
    else
      """
      mkdir -p ${params.project}/tmp/ ${params.project}/bin/
      cnest.py step2 --project ${params.project} --sample ${name} --input ${file_path} --fasta genome.fa --fast
      """
  }
  process make_bin_dir {

    input:
    path sample_bin_files from ch_bin_sample.collect()

    output:
    path "${params.project}/bin" into ch_bin

    script:
    """
    mkdir -p ${params.project}/bin
    for file in $sample_bin_files; do
      cp -L "\$file" ${params.project}/bin
    done
    """
  }
}

if (params.step =~ 3) {

  process step_3_gender_qc {
    //echo true
    publishDir "results/", mode: params.mode
    time '10h'

    input:
    path bin_dir from ch_bin
    path index from ch_index_tab

    output:
    path "gender_qc.txt"
    path "gender_classification.txt" into ch_gender
    path "mean_coverage.txt" into ch_cov

    script:
    """
      cnest.py step3 \
      --indextab $index \
      --bindir $bin_dir \
      --qc gender_qc.txt \
      --gender gender_classification.txt \
      --cov mean_coverage.txt
    """
  }
}

if (params.step =~ 4 || params.step =~ 5 || params.step =~ 6) {
  if (params.bindir) {
    ch_input_files = Channel.fromPath("${params.bindir}/*")
  } else if (params.step =~ 2) {
    def bin_path = ch_bin.view().val
    ch_input_files = Channel.fromPath("${bin_path}/*")
  } else {
    log.error "Please provide --bindir or run step 2 first"
    exit 1
  }
  if (params.rbindir) ch_input_files = Channel.fromPath("${params.rbindir}/*")

  println "Total number of samples in bin directory - "
  number_of_input_files = ch_input_files.count().view().val
  number_of_batches = (int) Math.ceil(number_of_input_files/params.target_size)
  println "Number of batches to run - "
  println number_of_batches

  if(params.start_batch > number_of_batches){
    log.error "--start_batch $params.start_batch must be less than or equal to $number_of_batches for this run."
    exit 1
  }

  if(params.run_until_n_batches && params.run_until_n_batches > number_of_batches ){
    log.error "--run_until_n_batches $params.run_until_n_batches must be less than or equal to $number_of_batches for this run."
    exit 1
  }

  if(params.run_until_n_batches){
    Channel
    .of( params.start_batch-1..params.run_until_n_batches-1)
    .map {
      (it * params.target_size) + 1
    }
    .into { ch_start_pos_1; ch_start_pos_2; ch_start_pos_3 }
  }else{
    Channel
    .of( params.start_batch-1..number_of_batches-1)
    .map {
      (it * params.target_size) + 1
    }
    .into { ch_start_pos_1; ch_start_pos_2; ch_start_pos_3 }
    // this above channel will produce 1, 11, 21, 31 ... for starting position
  }
}

if (params.step =~ 4){
  if (params.bindir) ch_bin = Channel.value(file(params.bindir))

  process step_4_logR_ratio {
    tag "start_pos_${start_pos}"
    //echo true
    publishDir "results/", mode: params.mode

    input:
    path bin_dir from ch_bin
    path index from ch_index_tab
    each start_pos from ch_start_pos_1

    output:
    path "${params.project}/cor/*" into ch_cor_files

    script:
    // for odd number of samples in a batch change target_size and batch_size
    def num_samples_in_current_batch =  number_of_input_files - start_pos
    if (num_samples_in_current_batch < params.target_size){
      batch_size = num_samples_in_current_batch + 1
      target_size = num_samples_in_current_batch + 1
    }else{
      batch_size = params.batch_size
      target_size = params.target_size
    }
    """
      mkdir -p ${params.project}/cor/
      cnest_dev.py step4 \
        --bindir $bin_dir \
        --indextab $index \
        --batch ${params.batch_size} \
        --tlen ${params.target_size} \
        --spos ${start_pos} \
        --cordir ${params.project}/cor/

      cat ${params.project}/cor/*
    """
  }

  process make_cor_dir {
    input:
    path ch_cor_files from ch_cor_files.collect()

    output:
    path "${params.project}/cor" into ch_cor_dir

    script:
    """
    mkdir -p ${params.project}/cor
    for file in $ch_cor_files; do
      cp -L "\$file" ${params.project}/cor
    done
    rm ${params.project}/cor/NA
    """
  }
}

if (params.step =~ 5){
  if (params.bindir) ch_bin = Channel.value(file(params.bindir))
  if (params.cordir) ch_cor_dir = Channel.fromPath("${params.cordir}")

  process step_5_log2_rbin_gen {
    tag "start_pos_${start_pos}_target_size_${target_size}"
    //echo true
    publishDir "results/", mode: params.mode

    input:
    path bin_dir from ch_bin
    path cor_dir from ch_cor_dir
    path index from ch_index_tab
    each start_pos from ch_start_pos_2
    path gender from ch_gender

    output:
    path "${params.project}/rbin/*" into ch_rbin_dir_files

    script:
    // for odd number of samples in a batch change target_size and batch_size
    def num_samples_in_current_batch =  number_of_input_files - start_pos
    if (num_samples_in_current_batch < params.target_size){
      batch_size = num_samples_in_current_batch + 1
      target_size = num_samples_in_current_batch + 1
    }else{
      batch_size = params.batch_size
      target_size = params.target_size
    }
    """
      echo "CPU = $task.cpus"
      echo "Memory = $task.memory"
      mkdir -p ${params.project}/rbin/
      cnest_dev.py step5 \
        --bindir $bin_dir \
        --cordir $cor_dir \
        --rbindir ${params.project}/rbin \
        --gender $gender \
        --indextab $index \
        --cor ${params.cor} \
        --batch $batch_size \
        --tlen $target_size \
        --spos ${start_pos} \
        --skipem
    """
  }

  process make_rbin_dir {
    input:
    path ch_rbin_dir_files from ch_rbin_dir_files.collect()

    output:
    path "${params.project}/rbin" into ch_rbin_dir

    script:
    """
    mkdir -p ${params.project}/rbin
    for file in $ch_rbin_dir_files; do
      cp -L "\$file" ${params.project}/rbin
    done
    """
  }
}

if (params.step =~ 6){
  if (params.bindir) ch_bin = Channel.value(file(params.bindir))
  if (params.cordir) ch_cor_dir = Channel.fromPath("${params.cordir}")
  if (params.rbindir) ch_rbin_dir = Channel.fromPath("${params.rbindir}")

  process step_6_hmm_call {
    tag "start_pos_${start_pos}"
    //echo true
    publishDir "results/", mode: params.mode

    input:
    path rbin_dir from ch_rbin_dir
    path cor_dir from ch_cor_dir
    path index from ch_index_tab
    path gender_file from ch_gender
    path cov_file from ch_cov
    each start_pos from ch_start_pos_3

    output:
    path "${params.project}/cnv/*"

    script:
    // for odd number of samples in a batch change target_size and batch_size
    def num_samples_in_current_batch =  number_of_input_files - start_pos
    if (num_samples_in_current_batch < params.target_size){
      batch_size = num_samples_in_current_batch + 1
      target_size = num_samples_in_current_batch + 1
    }else{
      batch_size = params.batch_size
      target_size = params.target_size
    }
    """
    mkdir -p ${params.project}/cnv/
    cnest_dev.py step6 \
      --rbindir $rbin_dir \
      --cordir $cor_dir \
      --cnvdir ${params.project}/cnv/ \
      --gender $gender_file \
      --indextab $index \
      --cov $cov_file \
      --covc ${params.covc} \
      --cor ${params.cor} \
      --batch $batch_size \
      --tlen $target_size \
      --spos ${start_pos} \
      --skipem
    """
  }
}
