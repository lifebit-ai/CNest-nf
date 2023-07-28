#!/usr/bin/env nextflow

// Re-usable componext for adding a helpful help message in our Nextflow script
def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf --project test --part 4

    Mandatory arguments:
      --project       [string] Name of the project
      --part          [int] Part of the workflow to run. One of [1, 2, 3, 4]


    Part 1 arguments:
      --ref           [file] Reference FASTA
      --indexb        [file] Index in BED format for fast counting
      --design        [file] A CSV file with header and three columns (name,cram,crai)

    Part 2 arguments:
      --index         [file] Index in tab format (index_tab.txt)
      --bindir        [path] Path to the directory of all bin files
      --binlist       [file] A txt file with paths to all bin files (one per row)

    Part 3 arguments:
      --index         [file] Index in tab format (index_tab.txt)
      --bindir        [path] Path to the directory of all bin files
      --gender        [file] Gender file from part 2
      --batch_size         [int]  Batch size for references
      --samples       [file] Samples to process

    Part 4 arguments:
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

if (params.wgs) {
  mem_factor = params.wgs
} else {
  mem_factor = 1
}
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


// Directories

// ! bindir and binlist are mutually exclusive

// Path to a folder of bin/rbin/cor files
if (params.bindir) ch_bin = Channel.value(file(params.bindir))
if (params.rbindir) ch_rbin = Channel.value(file(params.rbindir))
if (params.cordir) ch_cor = Channel.value(file(params.cordir))


// A txt file with one bin file per row
if (params.binlist) {
  Channel
    .fromPath(params.binlist)
    .splitText() { it.trim() }
    .collect()
    .set {ch_bins}
}


// Sample names
// If sample names are specified as a file, using it instead of all sample names in bin_dir or rbin_dir
if (params.samples) {
  Channel
    .fromPath(params.samples)
    .splitText() { it.trim() }
    .set { ch_sample_names }
}

if (!params.samples && params.bindir) {
  all_bins = file(params.bindir).list()
  ch_sample_names = Channel.from(all_bins)
}

if (!params.samples && params.rbindir) {
  all_rbins = file(params.rbindir).list()
  ch_sample_names = Channel.from(all_rbins)
}

if (!params.samples && params.binlist) {
  Channel
    .fromPath(params.binlist)
    .splitText() { file(it).baseName.trim() }
    .set { ch_sample_names }
}


// Helper files
if (params.index) ch_index = Channel.value(file(params.index))
if (params.indexb) ch_index_bed = Channel.value(file(params.indexb))
if (params.gender) ch_gender = Channel.value(file(params.gender))
if (params.cov) ch_cov = Channel.value(file(params.cov))

// Test mode
if (params.test == true) {subset_samples = 5} else {subset_samples = params.test}
if (params.test && params.design) ch_files_sets = ch_files_sets.take(subset_samples)
if (params.test && (params.bindir || params.binlist || params.rbindir || params.samples)) ch_sample_names = ch_sample_names.take(subset_samples)

/*
================================================================================
                                File staging
================================================================================
*/
if (params.binlist) {
  process stage_bins {
      echo true

      input:
        path bins from ch_bins

      output:
        file ("bin") into ch_bin

      shell:
      '''
      ls ./ > all_files
      mkdir -p bin
      cat ./all_files | while read f
      do
          mv $f bin/
      done
      rm bin/all_files
      '''
  }
}

/*
================================================================================
                                Main parts
================================================================================
*/
if (params.part == 0) {
  ch_bedgz = Channel.value(file("$baseDir/data/hg38.1kb.baits.bed.gz"))

  process step0 {
    tag "${params.project}"
    echo true

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
  process step1 {
    tag "${params.project}"
    publishDir "results/", mode: params.mode, pattern: "${params.project}/index*"
    echo true

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

if (params.part == 1) {
  process step2 {
    tag "id:${name}-file:${file_path}-index:${index_path}"
    publishDir "results/", mode: params.mode
    echo true

    input:
    set val(name), file(file_path), file(index_path) from ch_files_sets
    file("genome.fa") from ch_ref
    path "${params.project}/index.bed" from ch_index_bed

    output:
    path "${params.project}/bin/$name"

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
}

if (params.part == 2) {

  process gender_qc {
    echo true
    publishDir "results/", mode: params.mode
    time '10h'

    input:
    path bin_dir from ch_bin
    path index from ch_index

    output:
    path "gender_qc.txt"
    path "gender_classification.txt"
    path "mean_coverage.txt"

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

if (params.bindir) ch_bin_files = Channel.fromPath("${params.bindir}/*")

println "Total number of samples in bin directory - "
number_of_bin_files = ch_bin_files.count().view().val
number_of_batches = (int) Math.ceil(number_of_bin_files/params.batch_size)
println "Number of batches to run - "
println number_of_batches

if(params.run_first_n_batch_for_test){
  Channel
  .of( 0..params.run_first_n_batch_for_test-1)
  .map {
    (it * params.batch_size) + 1
  }
  .set { ch_start_pos }
}else{
  Channel
  .of( 0..number_of_batches-1)
  .map {
    (it * params.batch_size) + 1
  }
  .set { ch_start_pos }
  // this above channel will produce 1, 11, 21, 31 ... for starting position
}



if (params.part == 3) {
  process logR_ratio {
    tag "start_pos_${start_pos}"
    echo true
    publishDir "results/", mode: params.mode
//    memory { 1.GB * params.batch * mem_factor / 100 }
//    time { 40.m * params.batch * mem_factor / 100  }

    input:
    path bin_dir from ch_bin
    path index from ch_index
    val(start_pos) from ch_start_pos
    //path gender from ch_gender
    //val sample_name from ch_sample_names

    output:
    path "${params.project}/cor/*" into ch_cor_files

    script:
    """
      mkdir -p ${params.project}/cor/
      cnest_dev.py step4 \
        --bindir $bin_dir \
        --indextab $index \
        --batch ${params.batch_size} \
        --tlen ${params.target_size} \
        --spos ${start_pos} \
        --cordir ${params.project}/cor/
    """
  }
}

if (params.part == 4){
  process hmm_call {
    tag "${sample_name}"
    echo true
    publishDir "results/", mode: params.mode
    // memory { 5.GB * params.batch * mem_factor / 100 }
    // time { 40.m * params.batch * mem_factor / 100  }

    input:
    path rbin_dir from ch_rbin
    path cor_dir from ch_cor
    path index from ch_index
    path gender_file from ch_gender
    path cov_file from ch_cov
    val sample_name from ch_sample_names

    output:
    path "${params.project}/cnv/${sample_name}/${sample_name}_mixed_calls.txt"
    path "${params.project}/cnv/${sample_name}/${sample_name}_mixed_states.txt"

    script:
    """
      echo "Processing sample $sample_name"
      mkdir -p ${params.project}/cnv/
      cnest.py step5 \
        --indextab $index \
        --rbindir $rbin_dir \
        --cordir $cor_dir \
        --cnvdir ${params.project}/cnv/ \
        --cov    $cov_file \
        --sample $sample_name \
        --gender $gender_file \
        --batch $params.batch_size
    """
  }
}
