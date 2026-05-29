nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'



workflow version {
    println "RENEE_NF_ALEC ${workflow.manifest.version}"
}

workflow LOG {
    log.info """\
            RENEE_NF_ALEC $workflow.manifest.version
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


    yeet | view
}
