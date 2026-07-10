process WRITE_GENOME_CONFIG {
    label 'process_single'
    container "${params.containers.base}"

    input:
        tuple val(meta), path(fasta)
        tuple val(meta2), path(genes_gtf)
        path(gene_info)
        tuple val(meta4), path(star_index)
        val(organism)
        path(annotate)
        path(annotate_isoforms)
        path(refflat)
        path(bed_ref)
        path(qualimap_info)
        path(karyobeds)
        path(karyoploter)
        val(rsem_ref)
        path(rrna_list)
        path(tin_ref)
        path(fusion_blacklist)
        path(fusion_cytoband)
        path(fusion_protdomain)
        path(fusion_known_fusions)

    output:
        path("*.config"), emit: conf
        path("custom_genome/"), emit: files

    script:
    def genome_name = 'custom_genome'
    """
    #!/usr/bin/env python
    import os
    import shutil

    genome_name = "${genome_name}"
    os.makedirs(genome_name, exist_ok=True)

    def stage_path(src, destdir):
        if not src or not os.path.exists(src):
            return
        dst = os.path.join(destdir, os.path.basename(src))
        if os.path.isdir(src):
            shutil.copytree(src, dst)
        else:
            shutil.copy(src, dst)

    # Stage required files
    for f in ("${fasta}", "${genes_gtf}", "${gene_info}"):
        stage_path(f, genome_name)

    # Stage STAR index into a fixed subdirectory name
    if os.path.isdir("${star_index}"):
        shutil.copytree("${star_index}", os.path.join(genome_name, "star_index"))

    # Stage optional files / directories
    for f in ("${annotate}", "${annotate_isoforms}", "${refflat}",
              "${bed_ref}", "${qualimap_info}", "${karyobeds}",
              "${karyoploter}", "${rrna_list}", "${tin_ref}",
              "${fusion_blacklist}", "${fusion_cytoband}", "${fusion_protdomain}",
              "${fusion_known_fusions}"):
        stage_path(f, genome_name)

    idx = "\${params.index_dir}"
    genome = {
        "fasta":      f'"{idx}/{genome_name}/${fasta}"',
        "genes_gtf":  f'"{idx}/{genome_name}/${genes_gtf}"',
        "star_index": f'"{idx}/{genome_name}/star_index"',
        "gene_info":  f'"{idx}/{genome_name}/${gene_info}"',
        "organism":   '"${organism}"',
    }

    opt_paths = [
        ("${annotate}",          "annotate"),
        ("${annotate_isoforms}", "annotate_isoforms"),
        ("${refflat}",           "refflat"),
        ("${bed_ref}",           "bed_ref"),
        ("${qualimap_info}",     "qualimap_info"),
        ("${karyobeds}",         "karyobeds"),
        ("${karyoploter}",       "karyoploter"),
        ("${rrna_list}",         "rrna_list"),
        ("${tin_ref}",           "tin_ref"),
        ("${fusion_blacklist}",     "fusion_blacklist"),
        ("${fusion_cytoband}",      "fusion_cytoband"),
        ("${fusion_protdomain}",    "fusion_protdomain"),
        ("${fusion_known_fusions}", "fusion_known_fusions"),
    ]
    for src, key in opt_paths:
        if src and os.path.exists(src):
            suffix = "/" if os.path.isdir(src) else ""
            genome[key] = f'"{idx}/{genome_name}/{os.path.basename(src)}{suffix}"'

    if "${rsem_ref}":
        genome["rsem_ref"] = '"${rsem_ref}"'

    with open(f"{genome_name}.config", "w") as out:
        out.write("params {\\n")
        out.write('\\tindex_dir = "\${outputDir}/genome/"\\n')
        out.write("\\tgenomes {\\n")
        out.write(f"\\t\\t'{genome_name}' {{\\n")
        for k, v in genome.items():
            out.write(f"\\t\\t\\t{k:<20} = {v}\\n")
        out.write("\\t\\t}\\n")
        out.write("\\t}\\n")
        out.write("}\\n")
    """

    stub:
    """
    mkdir custom_genome/
    touch custom_genome.config custom_genome/genome.fa
    """
}
