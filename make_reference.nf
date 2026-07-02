nextflow.enable.dsl = 2
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { prepare_genome_workflow as PREPARE_GENOME }     from './subworkflows/local/prepare_genome/main.nf'

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
            """
            .stripIndent()
    log.info paramsSummaryLog(workflow)
}

workflow MAKE_REFERENCE {
    PREPARE_GENOME()
}

workflow {
    main:
        LOG()
        // validateParameters()

        PREPARE_GENOME()

        workflow.onComplete = {
            if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
                def message = Utils.spooker(workflow)
                if (message) {
                    println message
                }
            }
        }

    publish:
        prepare_genome_conf = PREPARE_GENOME.out.conf
}

output {
    prepare_genome_conf { path { file -> "genome/" } }
}
