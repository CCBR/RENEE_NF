nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include {FASTQC as FASTQC_RAW} from './modules/local/fastqc'
include {FASTQC as FASTQC_TRIMMED} from './modules/local/fastqc'
include {CUTADAPT} from './modules/CCBR/cutadapt'




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




    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id       : "${row.sample}_${row.replicate}",
                sample   : row.sample,
                replicate: row.replicate,
                qc_stage : 'raw' // can be used to track the stage of QC (e.g. raw, trimmed, etc.)
            ]

            tuple(meta, [file(row.fastq_1), file(row.fastq_2)])
        }

    // QC and trimming steps
    FASTQC_RAW(ch_reads)
    CUTADAPT(ch_reads)

    // update metadata for trimmed reads to reflect the new QC stage
    trimmed_reads = CUTADAPT.out.reads.map { meta, reads ->
        def qc_meta = meta + [qc_stage: 'trimmed']
        tuple(qc_meta, reads)
    }

    FASTQC_TRIMMED(trimmed_reads)


    workflow.onComplete = {
        if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
            def message = Utils.spooker(workflow)
            if (message) {
                println message
            }
        }
    }
}
