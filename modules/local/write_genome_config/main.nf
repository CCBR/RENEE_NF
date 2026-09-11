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
        path(rsem_ref_dir)
        path(rrna_list)
        path(tin_ref)
        path(fusion_blacklist)
        path(fusion_cytoband)
        path(fusion_protdomain)
        path(fusion_known_fusions)

    output:
        path("*.config"), emit: conf
        path("${genome_name}/"), emit: files

    script:
    genome_name = params.genome ?: 'custom_genome'
    // Absolute, so the generated config keeps working regardless of what
    // --outputDir a later run uses -- workflow.outputDir is the fully
    // resolved native output dir (see main.nf's sharedResourcesDir()),
    // matching exactly where this process's own `files` output gets
    // published (prepare_genome_conf -> "genome/" in main.nf's output{} block).
    def index_dir = "${workflow.outputDir}/genome/"
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

    # rsem_ref is a --reference *prefix*, not a single file/dir, so it's handled
    # separately from opt_paths above. An explicit rsem_ref string always wins
    # (matches conf/genomes/*.config, which point at a pre-existing external
    # RSEM reference and never stage it); otherwise, if BUILD_RSEM_REF produced
    # one, stage that whole directory and point at the prefix inside it.
    if "${rsem_ref}":
        genome["rsem_ref"] = '"${rsem_ref}"'
    elif os.path.isdir("${rsem_ref_dir}"):
        shutil.copytree("${rsem_ref_dir}", os.path.join(genome_name, "rsemref"))
        genome["rsem_ref"] = f'"{idx}/{genome_name}/rsemref/{genome_name}"'

    with open(f"{genome_name}.config", "w") as out:
        out.write("params {\\n")
        out.write('\\tindex_dir = "${index_dir}"\\n')
        out.write("\\tgenomes {\\n")
        out.write(f"\\t\\t'{genome_name}' {{\\n")
        for k, v in genome.items():
            out.write(f"\\t\\t\\t{k:<20} = {v}\\n")
        out.write("\\t\\t}\\n")
        out.write("\\t}\\n")
        out.write("}\\n")
    """

    stub:
    genome_name = params.genome ?: 'custom_genome'
    """
    mkdir ${genome_name}/
    touch ${genome_name}.config ${genome_name}/genome.fa
    """
}
