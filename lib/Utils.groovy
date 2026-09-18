import nextflow.Channel
import nextflow.Nextflow

class Utils {
    // Create a value channel for a configured path, or an empty channel when unset.
    public static optionalPathParam(params, String paramName) {
        return params.containsKey(paramName) && params[paramName]
            ? Channel.value(Nextflow.file(params[paramName], checkIfExists: true))
            : Channel.empty()
    }

    // Auto-detect an Arriba reference file for a genome build by name-matching,
    // mirroring RENEE (classic)'s workflow/rules/build.smk jsonmaker rule:
    // it infers the assembly from substrings in the genome name (hg19/hg38/
    // mm10/mm39 and their aliases) and looks up a preset path per assembly.
    // Here the "preset" is whatever file matching that assembly actually
    // exists in arribaDbDir, so it isn't pinned to one Arriba database
    // version. Returns null when arribaDbDir or genomeName is unset, the
    // name doesn't match a known assembly, or no matching file is found.
    public static String arribaReferenceFile(arribaDbDir, genomeName, String filePrefix, String fileSuffix) {
        if (!arribaDbDir || !genomeName) {
            return null
        }
        def name = genomeName.toString().toLowerCase()
        def build
        if (name.contains('hg19') || name.contains('hs37d') || name.contains('grch37')) {
            build = 'hg19_hs37d5_GRCh37'
        } else if (name.contains('hg38') || name.contains('hs38d') || name.contains('grch38')) {
            build = 'hg38_GRCh38'
        } else if (name.contains('mm10') || name.contains('grcm38')) {
            build = 'mm10_GRCm38'
        } else if (name.contains('mm39') || name.contains('grcm39')) {
            build = 'mm39_GRCm39'
        } else {
            return null
        }
        def dir = new File(arribaDbDir.toString())
        if (!dir.exists() || !dir.isDirectory()) {
            return null
        }
        def match = dir.listFiles()?.find { f ->
            f.name.startsWith(filePrefix) && f.name.contains(build) && f.name.endsWith(fileSuffix)
        }
        return match ? match.path : null
    }

    // run spooker for the workflow
    public static String spooker(workflow) {
        def pipeline_name = "${workflow.manifest.name.tokenize('/')[-1]}"
        def out = new StringBuilder()
        def err = new StringBuilder()
        def spooker_in_path = check_command_in_path("spooker")
        if (spooker_in_path) {
            try {
                println "Running spooker"
                def spooker_command = "spooker --outdir ${workflow.launchDir} --name ${pipeline_name} --version ${workflow.manifest.version} --path ${workflow.projectDir}"
                def command = spooker_command.execute()
                command.consumeProcessOutput(out, err)
                command.waitFor()
            } catch(IOException e) {
                err = e
            }
            new FileWriter("${workflow.launchDir}/log/spooker.log").with {
                write("${out}\n${err}")
                flush()
            }
        } else {
            err = "spooker not found, skipping"
        }
        return err
    }
    // check whether a command is in the path
    public static Boolean check_command_in_path(cmd) {
        def command_string = "command -V ${cmd}"
        def out = new StringBuilder()
        def err = new StringBuilder()
        try {
            def command = command_string.execute()
            command.consumeProcessOutput(out, err)
            command.waitFor()
        } catch(IOException e) {
            err = e
        }
        return err.length()==0

    }
}
