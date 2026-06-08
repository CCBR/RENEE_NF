nextflow.enable.dsl = 2

// Plugins
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include {FASTQC as FASTQC_RAW} from './modules/local/fastqc'
include {FASTQC as FASTQC_TRIMMED} from './modules/local/fastqc'
include {CUTADAPT} from './modules/CCBR/cutadapt'
include {STAR_ALIGN as STAR_ALIGN_PASS1} from './modules/nf-core/star/align'
include {STAR_ALIGN as STAR_ALIGN_PASS2} from './modules/nf-core/star/align'
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

process STAR_SJDB_FILTER {
    label 'process_low'

    input:
    path sj_tabs

    output:
    path 'uniq.filtered.SJ.out.tab', emit: sjdb

    script:
    """
    cat ${sj_tabs.join(' ')} | \
        sort | \
        uniq | \
        awk -F "\t" '{if ($5>0 && $6==1) {print}}' | \
        cut -f1-4 | \
        sort | \
        uniq | \
        grep "^chr" | \
        grep -v "^chrM" > uniq.filtered.SJ.out.tab
    """
}

workflow {
    main:
        LOG()
        validateParameters()

        def genome_cfg = (params.genomes instanceof Map && params.genome) ? params.genomes[params.genome] : null
        def star_index_path = params.star_index ?: genome_cfg?.star_index
        def star_gtf_path = params.star_gtf ?: genome_cfg?.genes_gtf

        if (!star_index_path) {
            error "Missing STAR index. Set --star_index or define params.genomes['${params.genome}'].star_index."
        }
        if (!star_gtf_path) {
            error "Missing STAR GTF. Set --star_gtf or define params.genomes['${params.genome}'].genes_gtf."
        }

        ch_star_index = Channel.value(tuple([id: params.genome ?: 'custom'], file(star_index_path, checkIfExists: true)))
        ch_star_gtf = Channel.value(tuple([id: params.genome ?: 'custom'], file(star_gtf_path, checkIfExists: true)))
        ch_sjdb_placeholder = ch_star_gtf.map { meta, gtf -> gtf }

        ch_reads = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def has_fastq_2 = row.fastq_2 && row.fastq_2.toString().trim()
                def meta = [
                    id       : "${row.sample}_${row.replicate}",
                    sample   : row.sample,
                    replicate: row.replicate,
                    layout   : has_fastq_2 ? 'paired' : 'single',
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

        ch_star_reads = CUTADAPT.out.reads
            .map { meta, reads -> tuple(meta + [single_end: meta.layout == 'single'], reads) }

        STAR_ALIGN_PASS1(ch_star_reads, ch_star_index, ch_star_gtf, false, false, ch_sjdb_placeholder)

        STAR_SJDB_FILTER(
            STAR_ALIGN_PASS1.out.spl_junc_tab
                .map { meta, sj -> sj }
                .collect()
        )

        ch_sjdb = STAR_SJDB_FILTER.out.sjdb.first()

        STAR_ALIGN_PASS2(ch_star_reads, ch_star_index, ch_star_gtf, false, true, ch_sjdb)

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
        star_pass1_sj = STAR_ALIGN_PASS1.out.spl_junc_tab
        star_pass1_log = STAR_ALIGN_PASS1.out.log_final
        star_sjdb = STAR_SJDB_FILTER.out.sjdb
        star_pass2_log = STAR_ALIGN_PASS2.out.log_final.mix(STAR_ALIGN_PASS2.out.log_out).mix(STAR_ALIGN_PASS2.out.log_progress)
        star_pass2_sj = STAR_ALIGN_PASS2.out.spl_junc_tab
        star_pass2_reads_per_gene = STAR_ALIGN_PASS2.out.read_per_gene_tab
        star_pass2_bam = STAR_ALIGN_PASS2.out.bam_sorted_aligned
        star_pass2_transcript_bam = STAR_ALIGN_PASS2.out.bam_transcript
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

    star_pass1_sj {
        path { meta, file -> 'STAR_files/' }
    }

    star_pass1_log {
        path { meta, file -> 'STAR_files/' }
    }

    star_sjdb {
        path { file -> 'STAR_files/' }
    }

    star_pass2_log {
        path { meta, file -> 'STAR_files/' }
    }

    star_pass2_sj {
        path { meta, file -> 'STAR_files/' }
    }

    star_pass2_reads_per_gene {
        path { meta, file -> 'STAR_files/' }
    }

    star_pass2_bam {
        path { meta, file -> 'STAR_files/' }
    }

    star_pass2_transcript_bam {
        path { meta, file -> 'bams/' }
    }
}
