#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/index_BC
========================================================================================
mpozuelo/index_BC Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/mpozuelo/index_BC
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/index_BC --input '*.txt' -profile docker

    Mandatory arguments
      --input [file]                Samplesheet with run and lane information to get the barcodes
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.
      --cluster_path                Cluster path to store data and to search for input data (Default: /datos/ngs/dato-activo)

    Optional:
      --single_index                In case the sequencing has been done single index although there were two indexes (i5 and i7). Only run with i7 and avoid the second index2 removal

    Demultiplexing parameters:
      --save_untrimmed              Saves untrimmed reads when demultiplexing (Default: FALSE)

    QC:
      --skipQC                      Skip all QC steps apart from MultiQC
      --skipFastQC                  Skip FastQC

    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
   custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }



// Validate inputs

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

if (!params.outdir) {
  params.outdir = params.run
}

cluster_path = params.cluster_path


// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Input'] = params.input
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-index_BC-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/index_BC Workflow Summary'
    section_href: 'https://github.com/mpozuelo/index_BC'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}




 process modify_samplesheet {
   publishDir "${cluster_path}/05_QC/${project}/samplesheet/", mode: params.publish_dir_mode

   input:
   path samplesheet from ch_input

   output:
   path "samplesheet_validated.csv" into ch_samplesheet

   script:
   out = "samplesheet_validated.csv"

   """
   modify_samplesheet.py $samplesheet $out
   """
 }



 def validate_input(LinkedHashMap line) {
   def run = line.run
   def lane = line.lane
   def bcsize = line.bc
   def indexsize = line.index
   def index2size = line.index2
   def libsize = line.size
   def fastq1 = line.fastq1
   def fastq2 = line.fastq2



   def array = []
   array = [ run, lane, [file(fastq1, checkIfExists: true), file(fastq2, checkIfExists: true)], bcsize, indexsize, index2size, libsize ]

   return array
 }

 /*
 * Create channels for input fastq files
 */
 ch_samplesheet
 .splitCsv(header:true, sep:',')
 .map { validate_input(it) }
 .into { ch_extra
         ch_index
         ch_index2
         ch_bc }


/*
 * STEP 1 - Change header and cell ranger
 */
//Detect index in the end of read2


process extraseq {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/02_pfastq/${platform}/${run}/${lane}/barcodes", mode: 'copy',

  input:
  set val(run), val(lane), path(reads), val(bcsize), val(indexsize), val(index2size), val(libsize) from ch_extra

  output:
  file("*.txt")

  script:

  seqsiz = libsize + 1

  """
  zcat ${reads[1]} | awk 'NR % 4 == 2' - | cut -c $seqsize- | sort | uniq -c | sort -k1,nr1 | head -1000 > "${run}_${lane}.read2.index.rank.txt"
  """

}


process index {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/02_pfastq/${platform}/${run}/${lane}/barcodes", mode: 'copy',

  input:
  set val(run), val(lane), path(reads), val(bcsize), val(indexsize), val(index2size), val(libsize) from ch_index

  output:
  file("*.txt")

  script:

  """
  zcat ${reads[1]} | awk 'NR % 4 == 2' - | rev | cut -c -$indexsize | rev | sort | uniq -c | sort -k1,nr1 | head -1000 > "${run}_${lane}.read2.index.rank.txt"
  """

}


process index2 {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/02_pfastq/${platform}/${run}/${lane}/barcodes", mode: 'copy',

  input:
  set val(run), val(lane), path(reads), val(bcsize), val(indexsize), val(index2size), val(libsize) from ch_index2

  output:
  file("*.txt")

  script:

  """
  zcat ${reads[1]} | awk 'NR % 4 == 2' - | rev | cut -c $indexsize- | cut -c -$index2size | rev | sort | uniq -c | sort -k1,nr1 | head -1000 > "${run}_${lane}.read2.index.rank.txt"
  """

}


process bc {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/02_pfastq/${platform}/${run}/${lane}/barcodes", mode: 'copy',

  input:
  set val(run), val(lane), path(reads), val(bcsize), val(indexsize), val(index2size), val(libsize) from ch_bc

  output:
  file("*.txt")

  script:

  """
  zcat ${reads[0]} | awk 'NR % 4 == 2' - | cut -c 1-$bcsize | sort | uniq -c | sort -k1,nr1 | head -1000 > "${run}_${lane}.read1.BC.${bcsize}bp.rank.txt"
  """

}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/index_BC] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/index_BC] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/index_BC]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/index_BC]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/index_BC v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
