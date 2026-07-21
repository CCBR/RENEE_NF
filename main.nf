nextflow.enable.dsl = 2

// Plugins
// Modules
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { SAMTOOLS_FLAGSTAT } from './modules/CCBR/samtools/flagstat/main.nf'
// Subworkflows
include { STAR_ALIGN }         from './subworkflows/local/star_align/main'
include { PREPARE_GENOME }     from './subworkflows/local/prepare_genome/main.nf'
include { INITIAL_QC }         from './subworkflows/local/initial_qc/main'
include { CHECK_INPUT }        from './subworkflows/local/read_samples/main'
include { PICARD_INITIAL_QC }  from './subworkflows/local/picard_initial_qc/main'
include { RSEQC_QC }           from './subworkflows/local/rseqc_qc/main'
include { RSEM }           from './subworkflows/local/rsem/main'
include { arriba as ARRIBA } from './subworkflows/local/arriba/main'

// Modules
include { BAM2STRANDEDBW } from './modules/local/bam2strandedbw/main'








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

workflow {
    main:
        LOG()
        validateParameters()

        // build genome first, if set to build mode stop after this step
        PREPARE_GENOME()

        if (params.build_genome) {
            log.info "Build genome only mode enabled. Stopping workflow after genome preparation."
            prepare_genome_conf = PREPARE_GENOME.out.conf



        } else {
            log.info "Genome preparation complete. Continuing with workflow."
            prepare_genome_conf = Channel.empty() // dont save genome copy if not in build mode

            // Read samplesheet and emit channel of reads
            CHECK_INPUT(params.input)

            // Initial QC, trimming, and FastQ Screen steps
            INITIAL_QC(CHECK_INPUT.out.reads)


            // STAR alignment steps ----------------------------------------------------------

            STAR_ALIGN(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf
            )

            // post-alignment steps ----------------------------------------------------------
            PICARD_INITIAL_QC(STAR_ALIGN.out.pass2_bam)

            // RSeQC QC: strandedness, read distribution, inner distance, TIN
            ch_bam_bai = PICARD_INITIAL_QC.out.bam
                .join(PICARD_INITIAL_QC.out.bai)


            RSEQC_QC(
                ch_bam_bai,
                PREPARE_GENOME.out.bed_ref,
                PREPARE_GENOME.out.tin_ref
            )

            // RSEM quantification -----------------------------------------------
            // strand_info is emitted as empty channel until RSeQC is integrated;
            // strandedness defaults to unstranded (--forward-prob 0.5)
            RSEM(
                STAR_ALIGN.out.pass2_transcript_bam,
                RSEQC_QC.out.infer_experiment,
                PREPARE_GENOME.out.rsem_ref,
                PREPARE_GENOME.out.annotate
            )
            // BAM to stranded BigWig files
            BAM2STRANDEDBW(
                PICARD_INITIAL_QC.out.bam
                    .join(PICARD_INITIAL_QC.out.bai)
                    .join(RSEQC_QC.out.infer_experiment)
            )
            // SAMTOOLS_FLAGSTAT expects a tuple: [ meta, bam, bai ]
            picard_bam_bai_ch = PICARD_INITIAL_QC.out.bam.join(PICARD_INITIAL_QC.out.bai)

            SAMTOOLS_FLAGSTAT(picard_bam_bai_ch)


            // Arriba gene-fusion calling (only when genome supplies a blacklist) ----------

            ARRIBA(
                INITIAL_QC.out.trimmed_reads,
                PREPARE_GENOME.out.star_index,
                PREPARE_GENOME.out.genes_gtf,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.fusion_blacklist,
                PREPARE_GENOME.out.fusion_known_fusions,
                PREPARE_GENOME.out.fusion_cytoband,
                PREPARE_GENOME.out.fusion_protdomain
            )

        }
        workflow.onComplete = {
            if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
                def message = Utils.spooker(workflow)
                if (message) {
                    println message
                }
            }
        }

    publish:
        // In build genome mode, only publish the genome conf file, otherwise publish all outputs
        prepare_genome_conf = params.build_genome ? prepare_genome_conf : Channel.empty()

        fastqc_raw     = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqc_raw
        fastqvalidator = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqvalidator
        cutadapt_reads = params.build_genome ? Channel.empty() : INITIAL_QC.out.cutadapt_reads
        cutadapt_log   = params.build_genome ? Channel.empty() : INITIAL_QC.out.cutadapt_log
        fastqc_trimmed = params.build_genome ? Channel.empty() : INITIAL_QC.out.fastqc_trimmed
        bbtools_ihist  = params.build_genome ? Channel.empty() : INITIAL_QC.out.bbtools_ihist
        fqscreen_1_txt = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_1_txt
        fqscreen_1_png = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_1_png
        fqscreen_2_txt = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_2_txt
        fqscreen_2_png = params.build_genome ? Channel.empty() : INITIAL_QC.out.fqscreen_2_png

        star_pass1_sj             = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass1_sj
        star_pass1_log            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass1_log
        star_sjdb                 = params.build_genome ? Channel.empty() : STAR_ALIGN.out.sjdb
        star_pass2_log            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_log
        star_pass2_sj             = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_sj
        star_pass2_reads_per_gene = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_reads_per_gene
        star_pass2_bam            = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_bam
        star_pass2_transcript_bam = params.build_genome ? Channel.empty() : STAR_ALIGN.out.pass2_transcript_bam

        rsem_genes_results    = params.build_genome ? Channel.empty() : RSEM.out.genes_results
        rsem_isoforms_results = params.build_genome ? Channel.empty() : RSEM.out.isoforms_results
        rsem_gene_counts      = params.build_genome ? Channel.empty() : RSEM.out.gene_counts
        rsem_gene_fpkm        = params.build_genome ? Channel.empty() : RSEM.out.gene_fpkm
        rsem_gene_tpm         = params.build_genome ? Channel.empty() : RSEM.out.gene_tpm
        rsem_isoform_counts   = params.build_genome ? Channel.empty() : RSEM.out.isoform_counts
        rsem_isoform_fpkm     = params.build_genome ? Channel.empty() : RSEM.out.isoform_fpkm
        rsem_isoform_tpm      = params.build_genome ? Channel.empty() : RSEM.out.isoform_tpm
        rsem_reformatted      = params.build_genome ? Channel.empty() : RSEM.out.reformatted
        rsem_gene_matrix      = params.build_genome ? Channel.empty() : RSEM.out.gene_matrix
        rsem_isoform_matrix   = params.build_genome ? Channel.empty() : RSEM.out.isoform_matrix

        picard_bam                = params.build_genome ? Channel.empty() : PICARD_INITIAL_QC.out.bam
        picard_bai                = params.build_genome ? Channel.empty() : PICARD_INITIAL_QC.out.bai

        rseqc_infer_experiment       = params.build_genome ? Channel.empty() : RSEQC_QC.out.infer_experiment
        rseqc_read_distribution      = params.build_genome ? Channel.empty() : RSEQC_QC.out.read_distribution
        rseqc_inner_distance_freq    = params.build_genome ? Channel.empty() : RSEQC_QC.out.inner_distance_freq
        rseqc_inner_distance_dist    = params.build_genome ? Channel.empty() : RSEQC_QC.out.inner_distance_dist
        rseqc_inner_distance_rscript = params.build_genome ? Channel.empty() : RSEQC_QC.out.inner_distance_rscript
        rseqc_tin_txt                = params.build_genome ? Channel.empty() : RSEQC_QC.out.tin_txt
        rseqc_tin_xls                = params.build_genome ? Channel.empty() : RSEQC_QC.out.tin_xls

        bam2bw_fwd = params.build_genome ? Channel.empty() : BAM2STRANDEDBW.out.fwd_bw
        bam2bw_rev = params.build_genome ? Channel.empty() : BAM2STRANDEDBW.out.rev_bw

        flagstat                  = params.build_genome ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.flagstat
        flagstat_versions         = params.build_genome ? Channel.empty() : SAMTOOLS_FLAGSTAT.out.versions

        arriba_fusions      = params.build_genome ? Channel.empty() : ARRIBA.out.fusions
        arriba_fusions_fail = params.build_genome ? Channel.empty() : ARRIBA.out.fusions_fail
        arriba_bam          = params.build_genome ? Channel.empty() : ARRIBA.out.bam
        arriba_pdf          = params.build_genome ? Channel.empty() : ARRIBA.out.pdf
        arriba_star_log     = params.build_genome ? Channel.empty() : ARRIBA.out.star_log
}

output {
    prepare_genome_conf { path { file -> "genome/" } }
    fastqc_raw { path { meta, file -> "fastqc/raw/" } }
    fastqvalidator { path { meta, file -> "fastqvalidator/${meta.id}/" } }
    cutadapt_reads { path { meta, reads -> "cutadapt/${meta.id}/" } }
    cutadapt_log { path { meta, log -> "cutadapt/${meta.id}/" } }
    fastqc_trimmed { path { meta, file -> "fastqc/trimmed/" } }
    bbtools_ihist { path { meta, ihist -> "bbtools/${meta.id}/" }}
    fqscreen_1_txt  { path { meta, file  -> "FQscreen/" } }
    fqscreen_1_png  { path { meta, file  -> "FQscreen/" } }
    fqscreen_2_txt  { path { meta, file  -> "FQscreen2/" } }
    fqscreen_2_png  { path { meta, file  -> "FQscreen2/" } }

    star_pass1_sj { path { meta, file -> 'STAR_files/pass1/' } }
    star_pass1_log { path { meta, file -> 'STAR_files/pass1/' } }
    star_sjdb { path { file -> 'STAR_files/pass1/' } }
    star_pass2_log { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_sj { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_reads_per_gene { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_bam { path { meta, file -> 'STAR_files/pass2/' } }
    star_pass2_transcript_bam { path { meta, file -> 'bams/' } }

    rsem_genes_results    { path { meta, file -> 'DEG_ALL/' } }
    rsem_isoforms_results { path { meta, file -> 'DEG_ALL/' } }
    rsem_gene_counts      { path { file -> 'DEG_ALL/' } }
    rsem_gene_fpkm        { path { file -> 'DEG_ALL/' } }
    rsem_gene_tpm         { path { file -> 'DEG_ALL/' } }
    rsem_isoform_counts   { path { file -> 'DEG_ALL/' } }
    rsem_isoform_fpkm     { path { file -> 'DEG_ALL/' } }
    rsem_isoform_tpm      { path { file -> 'DEG_ALL/' } }
    rsem_reformatted      { path { file -> 'DEG_ALL/' } }
    rsem_gene_matrix      { path { file -> 'DEG_ALL/' } }
    rsem_isoform_matrix   { path { file -> 'DEG_ALL/' } }

    picard_bam { path { meta, file -> 'bams/' } }
    picard_bai { path { meta, file -> 'bams/' } }

    rseqc_infer_experiment       { path { meta, file -> 'RSeQC/' } }
    rseqc_read_distribution      { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_freq    { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_dist    { path { meta, file -> 'RSeQC/' } }
    rseqc_inner_distance_rscript { path { meta, file -> 'RSeQC/' } }
    rseqc_tin_txt                { path { meta, file -> 'RSeQC/' } }
    rseqc_tin_xls                { path { meta, file -> 'RSeQC/' } }

    bam2bw_fwd { path { meta, file -> 'bigwigs/' } }
    bam2bw_rev { path { meta, file -> 'bigwigs/' } }

    flagstat { path { meta, file -> 'log_files/' } }
    flagstat_versions { path { file -> "log_files/versions/${file.getParent().getFileName()}/" } }

    arriba_fusions      { path { meta, file -> 'fusions/' } }
    arriba_fusions_fail { path { meta, file -> 'fusions/' } }
    arriba_bam          { path { meta, bam, bai -> 'fusions/' } }
    arriba_pdf          { path { meta, file -> 'fusions/' } }
    arriba_star_log     { path { meta, file -> 'STAR_files/arriba/' } }
}
