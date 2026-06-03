nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
<<<<<<< HEAD
include {FASTQC} from './modules/local/fastqc'
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
=======
>>>>>>> parent of c2f6680 (feat: Added first pass of fastQC to pipeline)



workflow version {
    println "RENEE_NF ${workflow.manifest.version}"
}

workflow LOG {
    log.info """\
            RENEE_NF $workflow.manifest.version
            =============
            cmd line     : $workflow.commandLine
            start time   : $workflow.start
            launchDir    : $workflow.launchDir
            input        : ${params.input}
            genome       : ${params.genome}
            """
            .stripIndent()
    log.info paramsSummaryLog(workflow)
}


process yeet {
    container "${params.containers.base}"

    output:
    stdout

    script:
    """
    echo ${params.input}
    """
}

workflow {
    LOG()
    validateParameters()

    workflow.onComplete = {
        if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
            def message = Utils.spooker(workflow)
            if (message) {
                println message
            }
        }
    }


    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id       : "${row.sample}_${row.replicate}",
                sample   : row.sample,
                replicate: row.replicate
            ]

            tuple(meta, [file(row.fastq_1), file(row.fastq_2)])
        }

    ch_reads | view


<<<<<<< HEAD
    FASTQC(ch_reads)
    FASTQC.out.html

    // Perform cutadapt on the reads
    CUTADAPT(ch_reads)


=======
    yeet | view
>>>>>>> parent of c2f6680 (feat: Added first pass of fastQC to pipeline)
}
