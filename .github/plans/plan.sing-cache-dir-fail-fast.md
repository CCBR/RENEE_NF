## Plan: Fail fast on missing Singularity cache dir

Nextflow only resolves and creates `singularity.cacheDir` lazily at the *first container pull*, so a missing `/data/$USER/singularity` surfaces late as a cryptic `Failed to create Singularity cache directory`. Worse, hardcoding `cacheDir` in the profile **silently overrides** `NXF_SINGULARITY_CACHEDIR` (config wins in Nextflow's resolution order), so the `/data/CCBR_Pipeliner/SIFS` value in [assets/slurm_header_biowulf.sh](assets/slurm_header_biowulf.sh#L13) is ignored — and that line is missing `export`, so it never reaches the nextflow process anyway.

Fix: make the profiles env-var-first, and add a startup check in `Utils` that mirrors Nextflow's resolution order (`singularity.cacheDir` → `NXF_SINGULARITY_CACHEDIR` → `${workDir}/singularity`) and errors with an actionable message.

**Steps**

1. Add `Utils.checkSingularityCacheDir(workflow)` to [lib/Utils.groovy](lib/Utils.groovy) — no-op unless `workflow.containerEngine in ['singularity','apptainer']` or `workflow.stubRun`. Resolve the effective dir, then `Nextflow.error` if: path is blank or has an empty `$USER` segment (`/data//singularity`), parent dir doesn't exist, path exists as a file, or path exists but isn't writable. Message includes the resolved path plus `mkdir -p <path>` and `export NXF_SINGULARITY_CACHEDIR=<path>` hints. Follows the existing static-method style of `spooker()`/`check_command_in_path()`.
2. Call it in [main.nf](main.nf#L52) inside `workflow { main: }`, right after `validateParameters()` and before the mutually-exclusive `--build` check — *depends on step 1*.
3. [conf/biowulf.config](conf/biowulf.config#L27): `cacheDir = System.getenv('NXF_SINGULARITY_CACHEDIR') ?: "/data/${System.getenv('USER')}/.singularity"` (note: also normalizes `singularity` → `.singularity` to match slurmint) — *parallel with step 1*.
4. Same env-var-first pattern in [conf/frce.config](conf/frce.config#L17) (fallback `/mnt/projects/CCBR-Pipelines/SIFs`) and [conf/slurmint.config](conf/slurmint.config#L13) — *parallel with step 3*.
5. Add `export` to `NXF_SINGULARITY_CACHEDIR` in [assets/slurm_header_biowulf.sh](assets/slurm_header_biowulf.sh#L13) and [assets/slurm_header_frce.sh](assets/slurm_header_frce.sh#L12) — *parallel*.
6. Add an nf-test asserting the helpful error when the cache dir's parent is missing; add a [CHANGELOG.md](CHANGELOG.md) entry — *depends on steps 1-2*.

**Relevant files**

- `lib/Utils.groovy` — new static validation method alongside `spooker`/`optionalPathParam`
- `main.nf` — invoke check at workflow entry
- `conf/biowulf.config`, `conf/frce.config`, `conf/slurmint.config` — `singularity.cacheDir` blocks
- `assets/slurm_header_biowulf.sh`, `assets/slurm_header_frce.sh` — missing `export`
- `CHANGELOG.md` — bug fix entry

**Verification**

1. `NXF_SINGULARITY_CACHEDIR=/nonexistent/parent/cache nextflow run main.nf -profile biowulf,test -preview` → clear, actionable error before any process launches.
2. Unset env var + valid `/data/$USER` → run proceeds unchanged.
3. `export NXF_SINGULARITY_CACHEDIR=/data/CCBR_Pipeliner/SIFS; nextflow config -profile biowulf | grep cacheDir` → shows the shared SIFS path (proves the override bug is fixed).
4. `nf-test test` and `pytest tests/test_cli.py` still pass.

**Decisions**

- Env-var-first with per-user fallback; no auto-`mkdir` — the pipeline reports, the user decides where their cache lives.
- Check is Groovy-side only, so it applies to bare `nextflow run` as well as `renee_nf run`. The Python CLI (`src/__main__.py`) and `ccbr_tools.pipeline.cache` helpers are deliberately left untouched.
- Read-only shared caches (e.g. `/data/CCBR_Pipeliner/SIFS`) must still pass the check — require readable, not writable, when the dir already exists and is not user-owned.

**Further Considerations**

1. Should the check also verify `singularity.libraryDir` / `NXF_SINGULARITY_LIBRARYDIR` when set? Option A: validate both (more thorough). Option B: cache dir only (recommended — libraryDir is optional and read-only).
2. Should the resolved cache dir be printed in the `LOG` workflow banner alongside `launchDir`? Recommended yes — cheap and makes support requests easier to triage.
