process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "nciccbr/ccbr_kraken_v2.1.1:v0.0.1"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt')                  , emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                           , emit: report
    tuple val(meta), path('*krona.html')                           , emit: krona_html
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions, emit: versions_kraken2
    tuple val("${task.process}"), val('KronaTools'), eval('ktImportTaxonomy 2>&1 | sed -nE "s/.*KronaTools ([0-9]+(\\.[0-9]+)*).*/\\1/p"'), topic: versions, emit: versions_kronatools
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | head -1 | sed "s/gzip //g"'), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired       = meta.single_end ? "" : "--paired"
    def classified   = meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""

    // the original zip method was "pigz -p $task.cpus *.fastq"
    // This is a parallel compression method that is faster and more efficient than gzip
    // TODO : Include pigz into our container for faster use. For now, we will use gzip which is slower and less efficient
    def compress_reads_command = save_output_fastqs ? "gzip *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --output ${prefix}.kraken2.classifiedreads.txt \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $paired \\
        $args \\
        $reads

    # extract Read ID and taxonomy ID from the Kraken2 output file
    cut -f2,3 ${prefix}.kraken2.classifiedreads.txt | \
        ktImportTaxonomy - -o ${prefix}.kraken2.krona.html

    # Compress the output fastq files if they were generated
    $compress_reads_command
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def classified   = meta.single_end ? "${prefix}.classified.fastq.gz"   : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"

    """
    touch ${prefix}.kraken2.report.txt
    touch ${prefix}.kraken2.out.txt
    touch ${prefix}.kraken2.krona.html
    touch ${prefix}.kraken2.classifiedreads.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    """

}
