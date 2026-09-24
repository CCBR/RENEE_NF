process BUILD_ANNOTATE {
    label 'process_single'
    container "${params.containers.build_rnaseq}"

    input:
        path(gtf)

    output:
        path("annotate.genes.txt"),            emit: genes
        path("annotate.isoforms.txt"),         emit: isoforms
        path("refFlat.txt"),                   emit: refflat
        path("genes.ref.bed"),                 emit: bed_ref
        path("geneinfo.bed"),                  emit: gene_info
        path("karyobeds"),                     emit: karyobeds
        path("karyoplot_gene_coordinates.txt"), emit: karyoploter

    when:
        task.ext.when == null || task.ext.when

    script:
    // Ports RENEE (classic) workflow/rules/build.smk rules `annotate`,
    // `karyo_beds`, and `karyo_coord` -- combined into one process since all
    // three take just the GTF as input. make_refFlat.py and make_geneinfo.py
    // take no arguments -- they read annotate.isoforms.txt / annotate.genes.txt
    // and genes.ref.bed from the current directory, so the first five outputs
    // have to be produced together, in this order.
    """
    get_gene_annotate.py ${gtf} > annotate.genes.txt
    get_isoform_annotate.py ${gtf} > annotate.isoforms.txt
    gtfToGenePred -ignoreGroupsWithoutExons ${gtf} genes.genepred
    genePredToBed genes.genepred genes.bed12
    sort -k1,1 -k2,2n genes.bed12 > genes.ref.bed
    make_refFlat.py > refFlat.txt
    make_geneinfo.py ${gtf} > geneinfo.bed

    mkdir -p karyobeds
    cd karyobeds
    get_karyoplot_beds.py ../${gtf}
    cd ..

    get_karyoplot_gene_coordinates.py ${gtf} > karyoplot_gene_coordinates.txt
    """

    stub:
    """
    touch annotate.genes.txt annotate.isoforms.txt refFlat.txt genes.ref.bed geneinfo.bed
    mkdir -p karyobeds
    touch karyobeds/karyobed.bed
    touch karyoplot_gene_coordinates.txt
    """
}
