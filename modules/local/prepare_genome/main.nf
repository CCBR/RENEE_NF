
process GTF2BED {
    tag { gtf }
    label 'process_single'

    container "${params.containers.base}"

    input:
        path(gtf)

    output:
        path("${gtf.baseName}.bed"), emit: bed

    script:
    """
    # remove comment lines and convert GTF to BED format.
    # looks for transcript_id for test case
    grep -v "^#" ${gtf} | grep 'transcript_id' | gtf2bed > ${gtf.baseName}.bed

    """
    stub:
    """
    touch ${gtf.baseName}.bed
    """
}
process SPLIT_REF_CHROMS {
    tag { fasta }
    label 'process_single'
    container "${params.containers.base}"

    input:
        path(fasta)

    output:
        path("${fasta.baseName}.chrom.sizes"), emit: chrom_sizes
        path("chroms/")                      , emit: chrom_dir

    script:
    """
    splitRef.py ${fasta} ${fasta.baseName}.chrom.sizes chroms
    """
    stub:
    """
    touch ${fasta.baseName}.chrom.sizes
    mkdir -p chroms/
    touch chroms/chr1.fa
    """
}

process RENAME_FASTA_CONTIGS {
    """
    Convert ensembl to UCSC contig names in a fasta file
    """
    tag { fasta }

    container "${params.containers.base}"

    input:
        path(fasta)
        path(map)

    output:
        path("*_renamed*"), emit: fasta

    script:
    def renamed_fasta= "${fasta.getSimpleName()}_renamed.${fasta.getExtension()}"
    """
    #!/usr/bin/env python

    with open("${map}", "r") as mapfile:
        contig_map = {line.strip().split()[0]: line.strip().split()[1] for line in mapfile}

    with open("${fasta}", "r") as infile:
        with open("${renamed_fasta}", "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    old_contig = line.strip(">").strip()
                    contig = contig_map[old_contig] if old_contig in contig_map.keys() else old_contig
                    outfile.write(f">{contig}\\n")
                else:
                    outfile.write(line)

    """

    stub:
    def renamed_fasta= "${fasta.getSimpleName()}_renamed.${fasta.getExtension()}"
    """
    touch ${renamed_fasta}
    """
}

process RENAME_DELIM_CONTIGS {
    """
    Convert ensembl to UCSC contig names in a delimited file (e.g. GTF, BED)
    using cvbio https://github.com/clintval/cvbio#updatecontignames
    """
    tag { delim }

    container "nciccbr/ccbr_cvbio_3.0.0:v1.0.1"

    input:
        path(delim)
        path(map)

    output:
        path("*_renamed*"), emit: delim

    script:
    def renamed_delim = "${delim.getSimpleName()}_renamed.${delim.getExtension()}"
    """
    cvbio UpdateContigNames \\
        -i ${delim} \\
        -o ${renamed_delim} \\
        -m ${map} \\
        --comment-chars '#' \\
        --columns 0 \\
        --skip-missing false
    """

    stub:
    def renamed_delim = "${delim.getSimpleName()}_renamed.${delim.getExtension()}"
    """
    touch ${renamed_delim}
    """
}

process WRITE_GENOME_CONFIG {
    label 'process_single'
    container "${params.containers.base}"

    input:
        path(fasta)
        path(genes_gtf)
        path(gene_info)
        path(star_index)
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
              "${fusion_blacklist}", "${fusion_cytoband}", "${fusion_protdomain}"):
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
        ("${fusion_blacklist}",  "fusion_blacklist"),
        ("${fusion_cytoband}",   "fusion_cytoband"),
        ("${fusion_protdomain}", "fusion_protdomain"),
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
