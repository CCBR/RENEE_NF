process MULTIQCPARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "nciccbr/ccbr_multiqc_1.15:v2"

    input:
    tuple val(meta), path(multiqc_data, stageAs: "multiqc_data"), path(inner_distance_freq, stageAs: "inner_distance/*"), path(tin_summary, stageAs: "tin/*"), path(fastq_info, stageAs: "fastq_info/*")

    output:
    tuple val(meta), path("multiqc_matrix.tsv")        , emit: matrix
    tuple val(meta), path("rseqc_inner_distances.txt"), emit: inner_distances
    tuple val(meta), path("rseqc_median_tin.txt")     , emit: median_tin
    tuple val(meta), path("fastq_flowcell_lanes.txt") , emit: flowcell_lanes
    tuple val("${task.process}"), val('python'), eval('python3 --version 2>&1 | cut -d " " -f 2'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p inner_distance tin fastq_info

    printf 'Sample\\tInner_Dist_Maxima\\n' > rseqc_inner_distances.txt
    while IFS= read -r -d '' f; do
        sample=\$(basename "\${f}")
        inner_dist_maxima=\$(sort -k3,3nr "\${f}" | awk -F '\\t' 'NR == 1 { print \$1 }')
        printf '%s\\t%s\\n' "\${sample}" "\${inner_dist_maxima}"
    done < <(find -L inner_distance -maxdepth 1 -type f -print0 | sort -z) \\
        >> rseqc_inner_distances.txt

    printf 'Sample\\tmedian_tin\\n' > rseqc_median_tin.txt
    while IFS= read -r -d '' f; do
        awk -F '\\t' '\$1 != "Bam_file" { printf "%s\\t%.3f\\n", \$1, \$3 }' "\${f}"
    done < <(find -L tin -maxdepth 1 -type f -print0 | sort -z) \\
        | LC_ALL=C sort -t \$'\\t' -k1,1 \\
        >> rseqc_median_tin.txt

    printf 'Sample\\tflowcell_lanes\\n' > fastq_flowcell_lanes.txt
    while IFS= read -r -d '' f; do
        awk -F '\\t' -v OFS='\\t' 'NR == 2 { print \$1, \$5 }' "\${f}"
    done < <(find -L fastq_info -maxdepth 1 -type f -print0 | sort -z) \\
        | LC_ALL=C sort -t \$'\\t' -k1,1 \\
        >> fastq_flowcell_lanes.txt

    python3 ${moduleDir}/resources/usr/bin/pyparser.py \\
        multiqc_data/*.txt \\
        rseqc_inner_distances.txt \\
        rseqc_median_tin.txt \\
        fastq_flowcell_lanes.txt \\
        .
    """

    stub:
    """
    touch multiqc_matrix.tsv
    printf 'Sample\\tInner_Dist_Maxima\\n' > rseqc_inner_distances.txt
    printf 'Sample\\tmedian_tin\\n' > rseqc_median_tin.txt
    printf 'Sample\\tflowcell_lanes\\n' > fastq_flowcell_lanes.txt
    """
}
