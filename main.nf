nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include {FASTQC as FASTQC_RAW} from './modules/local/fastqc'
include {FASTQC as FASTQC_TRIMMED} from './modules/local/fastqc'
include {CUTADAPT} from './modules/CCBR/cutadapt'
include {validate_fastqs as VALIDATE_FASTQS} from './subworkflows/local/validate_fastqs/main'




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
    main:
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
                    fastq_1  : new File(row.fastq_1.toString()).name,
                    fastq_2  : new File(row.fastq_2.toString()).name
                ]

                tuple(meta, [file(row.fastq_1), file(row.fastq_2)])
            }
        // Sample validation gate
        VALIDATE_FASTQS(ch_reads)
        ch_validated_reads = VALIDATE_FASTQS.out.reads


        // QC and trimming steps
        FASTQC_RAW(ch_validated_reads)

        CUTADAPT(ch_validated_reads)

        FASTQC_TRIMMED(CUTADAPT.out.reads)

        workflow.onComplete = {
            if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
                def message = Utils.spooker(workflow)
                if (message) {
                    println message
                }
            }
        }

    publish:
        fastqc_raw = FASTQC_RAW.out.html.mix(FASTQC_RAW.out.zip)
        fastqvalidator = VALIDATE_FASTQS.out.report.mix(VALIDATE_FASTQS.out.log_r1).mix(VALIDATE_FASTQS.out.log_r2)
        cutadapt_reads = CUTADAPT.out.reads
        cutadapt_log = CUTADAPT.out.log
        fastqc_trimmed = FASTQC_TRIMMED.out.html.mix(FASTQC_TRIMMED.out.zip)
}

output {
    fastqc_raw {
        path { meta, file -> "fastqc/raw/" }
    }

    fastqvalidator {
        path { meta, file -> "fastqvalidator/${meta.id}/" }
    }

    cutadapt_reads {
        path { meta, reads -> "cutadapt/${meta.id}/" }
    }

    cutadapt_log {
        path { meta, log -> "cutadapt/${meta.id}/" }
    }

    fastqc_trimmed {
        path { meta, file -> "fastqc/trimmed/" }
    }
}
