nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include {FASTQC as FASTQC_RAW} from './modules/local/fastqc'
include {FASTQC as FASTQC_TRIMMED} from './modules/local/fastqc'
include {BBTOOLS_BBMERGE} from './modules/local/bbtools'
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
                def has_fastq_2 = row.fastq_2 && row.fastq_2.toString().trim()
                def meta = [
                    id       : "${row.sample}_${row.replicate}",
                    sample   : row.sample,
                    replicate: row.replicate,
                    single_end   : !has_fastq_2
                ]

                def reads = [file(row.fastq_1)]
                if (has_fastq_2) {
                    reads << file(row.fastq_2)
                }

                tuple(meta, reads)
            }
        // Split each sample read list into one fastq per emitted tuple for validation.
        individual_fastq_ch = ch_reads.transpose()


        // Sample validation gate
        VALIDATE_FASTQS(individual_fastq_ch, ch_reads)
        ch_validated_reads = VALIDATE_FASTQS.out.reads

        // ch_validated_reads.view()

        // QC and trimming steps
        FASTQC_RAW(ch_validated_reads)

        CUTADAPT(ch_validated_reads)

        FASTQC_TRIMMED(CUTADAPT.out.reads)

        BBTOOLS_BBMERGE(CUTADAPT.out.reads)

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
        fastqvalidator = VALIDATE_FASTQS.out.logs
        cutadapt_reads = CUTADAPT.out.reads
        cutadapt_log = CUTADAPT.out.log
        fastqc_trimmed = FASTQC_TRIMMED.out.html.mix(FASTQC_TRIMMED.out.zip)
        bbtools_ihist = BBTOOLS_BBMERGE.out.ihist
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

    bbtools_ihist {
        path { meta, ihist -> "bbtools/${meta.id}/" }
    }
}
