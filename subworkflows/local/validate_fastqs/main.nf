include { FASTQVALIDATOR } from "../../../modules/local/fastqvalidator"

workflow validate_fastqs {
    take:
        ch_reads

    main:
        FASTQVALIDATOR(ch_reads)

        ch_validator_check = FASTQVALIDATOR.out.status
            .map { meta, status_file ->
                [
                    id    : meta.id,
                    fastq1: meta.fastq_1,
                    fastq2: meta.fastq_2,
                    status: status_file.text.trim()
                ]
            }
            .collect()
            .map { statuses ->
                def failed = statuses.findAll { entry ->
                    !entry.status.startsWith('PASS')
                }

                if (failed) {
                    def failed_summary = failed.collect { entry ->
                        "${entry.id}: ${entry.fastq1}, ${entry.fastq2} (${entry.status})"
                    }.join(',\n')
                    error "FASTQVALIDATOR failed for FASTQ file(s): ${failed_summary}"
                }

                true
            }

        // Using collect will pause until validation is done and check results before proceeding
        ch_validated_reads = ch_reads
            .collect(flat: false)
            .flatMap { reads -> reads }

    emit:
        reads = ch_validated_reads
        report = FASTQVALIDATOR.out.report
        log_r1 = FASTQVALIDATOR.out.log_r1
        log_r2 = FASTQVALIDATOR.out.log_r2
}
