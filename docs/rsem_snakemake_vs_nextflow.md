# RSEM — Snakemake vs Nextflow Comparison

## Overview

| Step                      | Snakemake rule            | Nextflow equivalent                           | Status                   |
| ------------------------- | ------------------------- | --------------------------------------------- | ------------------------ |
| Per-sample quantification | `rsem` (PE) / `rsem` (SE) | `RSEM_CALCULATEEXPRESSION` (nf-core, patched) | ✅ Implemented           |
| Annotated count matrices  | `rsem_merge`              | `RSEM_MERGE` (local)                          | ✅ Implemented           |
| Raw count matrices        | `rsem_data_matrix`        | `RSEM_GENERATE_DATA_MATRIX` (local)           | ✅ Implemented           |
| Sample table              | `samplecondition`         | —                                             | ⏭ Not a pipeline target |
| Filtered count table      | `rsemcounts`              | —                                             | ⏭ Not a pipeline target |
| PCA report                | `pca`                     | —                                             | ⏭ Not a pipeline target |

---

## 1. Per-sample quantification

### Snakemake — `paired-end.smk`

```python
rule rsem:
    input:
        file1=join(workpath, bams_dir, "{name}.p2.Aligned.toTranscriptome.out.bam"),
        file2=join(workpath, rseqc_dir, "{name}.strand.info")
    output:
        out1=join(workpath, degall_dir, "{name}.RSEM.genes.results"),
        out2=join(workpath, degall_dir, "{name}.RSEM.isoforms.results"),
    params:
        prefix=join(workpath, degall_dir, "{name}.RSEM"),
        rsemref=config['references'][pfamily]['RSEMREF'],
        tmpdir=tmpdir,
    threads: int(allocated("threads", "rsem", cluster)),  # 32 threads
    container: config['images']['rsem']                   # nciccbr/ccbr_rsem_1.3.3:v1.0
    shell: """
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT

    fp=$(tail -n1 {input.file2} | awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}')

    rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads} {input.file1} {params.rsemref} {params.prefix} --time \
        --temporary-folder ${{tmp}} --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
    """
```

**Single-end** is identical except `--paired-end` is omitted.

### Nextflow — `modules/nf-core/rsem/calculateexpression/main.nf` (patched)

```groovy
process RSEM_CALCULATEEXPRESSION {
    label 'process_high'                              // 32 CPUs via base.config
    container "nciccbr/ccbr_rsem_1.3.3:v1.0"        // same container as Snakemake

    input:
    tuple val(meta), path(reads)   // toTranscriptome BAM; meta.strandedness set by subworkflow
    path  index                    // parent dir of rsem_ref prefix (staged)

    output:
    tuple val(meta), path("*.genes.results"),    emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    ...

    script:
    def strandedness = meta.strandedness == 'forward' ? '--strandedness forward'
                     : meta.strandedness == 'reverse' ? '--strandedness reverse'
                     : ''                             // equivalent to --forward-prob 0.5

    def alignment_mode = '--alignments'              // always BAM input
    def paired_end     = meta.single_end ? '' : '--paired-end'

    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`
    rsem-calculate-expression \
        --num-threads $task.cpus \
        --temporary-folder ./tmp/ \
        $alignment_mode \
        $strandedness \
        $paired_end \
        $args \          # --no-bam-output --calc-ci --seed 12345 --time --keep-intermediate-files --estimate-rspd
        $reads \
        $INDEX \
        $prefix
    """
}
```

**Strandedness mapping** (`strand.info` last-line last-field → `meta.strandedness`):

| `infer_experiment.py` fraction | Snakemake `--forward-prob` | Nextflow `--strandedness` |
| ------------------------------ | -------------------------- | ------------------------- |
| > 0.75                         | 0.0                        | `reverse`                 |
| < 0.25                         | 1.0                        | `forward`                 |
| 0.25–0.75                      | 0.5                        | _(omitted = unstranded)_  |

**ext.args in `conf/modules.config`:**

```groovy
withName: RSEM_CALCULATEEXPRESSION {
    ext.args   = '--no-bam-output --calc-ci --seed 12345 --time --keep-intermediate-files --estimate-rspd'
    ext.prefix = { "${meta.id}.RSEM" }
}
```

---

## 2. Annotated count matrices

### Snakemake — `common.smk`

```python
rule rsem_merge:
    input:
        files =expand(..., "{name}.RSEM.genes.results",    name=samples),
        files2=expand(..., "{name}.RSEM.isoforms.results", name=samples),
    output:
        gene_counts_matrix = "RSEM.genes.expected_count.all_samples.txt",
        gene_fpkm_matrix   = "RSEM.genes.FPKM.all_samples.txt",
        isoform_fpkm_matrix= "RSEM.isoforms.FPKM.all_samples.txt",
        reformatted        = "RSEM.genes.expected_counts.all_samples.reformatted.tsv",
    container: config['images']['base']
    shell: """
    python merge_rsem_results.py {params.annotate} {params.inputdir} {params.inputdir}
    sed 's/\t/|/1' RSEM.genes.expected_count.all_samples.txt | \
        sed '1 s/^gene_id|GeneName/symbol/' > {output.reformatted}
    """
```

### Nextflow — `modules/local/rsem_merge/main.nf`

```groovy
process RSEM_MERGE {
    container "${params.containers.base}"   // same base container
    input:
    path gene_results      // collected *.RSEM.genes.results
    path isoform_results   // collected *.RSEM.isoforms.results
    path annotate          // annotate.genes.txt
    output:
    path "RSEM.genes.expected_count.all_samples.txt",          emit: gene_counts
    path "RSEM.genes.FPKM.all_samples.txt",                    emit: gene_fpkm
    path "RSEM.genes.TPM.all_samples.txt",                     emit: gene_tpm
    path "RSEM.isoforms.expected_count.all_samples.txt",       emit: isoform_counts
    path "RSEM.isoforms.FPKM.all_samples.txt",                 emit: isoform_fpkm
    path "RSEM.isoforms.TPM.all_samples.txt",                  emit: isoform_tpm
    path "RSEM.genes.expected_counts.all_samples.reformatted.tsv", emit: reformatted
    script:
    """
    merge_rsem_results.py ${annotate} . .
    sed 's/\t/|/1' RSEM.genes.expected_count.all_samples.txt | \
        sed '1 s/^gene_id|GeneName/symbol/' > RSEM.genes.expected_counts.all_samples.reformatted.tsv
    """
}
```

> Note: Nextflow also outputs TPM matrices — `merge_rsem_results.py` produces them but Snakemake's rule didn't declare them explicitly.

---

## 3. Raw count matrices

### Snakemake — `common.smk`

```python
rule rsem_data_matrix:
    input:
        genes   =expand(..., "{name}.RSEM.genes.results",    name=samples),
        isoforms=expand(..., "{name}.RSEM.isoforms.results", name=samples),
    output:
        genes   ="RSEM.genes.expected_counts.all_samples.matrix",
        isoforms="RSEM.isoforms.expected_counts.all_samples.matrix",
    container: config['images']['rsem']
    shell: """
    rsem-generate-data-matrix {input.genes}    > {output.genes}
    rsem-generate-data-matrix {input.isoforms} > {output.isoforms}
    """
```

### Nextflow — `modules/local/rsem_generate_data_matrix/main.nf`

```groovy
process RSEM_GENERATE_DATA_MATRIX {
    container "nciccbr/ccbr_rsem_1.3.3:v1.0"   // same RSEM container
    input:
    path gene_results
    path isoform_results
    output:
    path "RSEM.genes.expected_counts.all_samples.matrix",    emit: gene_matrix
    path "RSEM.isoforms.expected_counts.all_samples.matrix", emit: isoform_matrix
    script:
    """
    rsem-generate-data-matrix ${gene_results.join(' ')}    > RSEM.genes.expected_counts.all_samples.matrix
    rsem-generate-data-matrix ${isoform_results.join(' ')} > RSEM.isoforms.expected_counts.all_samples.matrix
    """
}
```

---

## 4. Downstream DEG rules (not pipeline targets — omitted from Nextflow)

These rules are in `group-info.smk` but **absent from `rule all`** in the Snakefile. They require group/condition metadata and a DESeq2/edgeR R environment not yet in RENEE_NF.

### `samplecondition` — builds DESeq2 sample table

```python
# Writes sampletable.txt: sampleName / fileName / condition / label
# Requires: config['project']['groups']['rgroups/rlabels/rsamps']
```

### `rsemcounts` — CPM filtering + normalization

```python
# Runs rsemcounts.R with DESeq2/edgeR
# Outputs:
#   RawCountFile_RSEM_genes.txt
#   RawCountFile_RSEM_genes_filtered.txt   ← CPM-filtered
#   RSEM_CPM_counts.txt.gz
#   RSEM_CPM_TMM_counts.txt.gz
#   RSEM_CPM_TMM_unfiltered_counts.txt.gz
#   RSEM_rlog_counts.txt.gz
```

### `pca` — PCA HTML report

```python
# Runs pcacall.R → PcaReport.Rmd
# Requires: sampletable.txt + RawCountFile_{dtype}_filtered.txt
```
