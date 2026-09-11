# RFMix-reader v0.7.0: benchmark revision and application analysis plan

Target: revised `rfmix_reader-benchmarking` repository + manuscript submitted within 8 weeks.

Target journal: either **HGG Advances** or **Bioinformatics Advances**. The core analysis is shared; final framing will depend on whether the strongest manuscript contribution is the human-genetics application (HGG Advances) or the validated software/data-layer advance (Bioinformatics Advances).

Prepared 2026-09-11.

---

## 0. Starting state and release strategy

### 0.1 Software

`heart-gen/rfmix_reader` is currently at `0.6.0` on `main` (`pyproject.toml`), GPL-3.0-or-later, after a four-stage rebuild recorded in `CHANGELOG.md`:

| Release | Change |
|---|---|
| 0.3.2 | fixed silent data corruption in the `.fb.tsv` path; ancestry axis follows the tool's header order |
| 0.4.0 | `open_rfmix` / `open_flare` / `open_simu` / `open_local_ancestry` / `convert` returning a lazy `xarray.Dataset`; streaming parsers; per-chromosome Zarr cache |
| 0.5.0 | Dataset operations (`to_bed`, `at_positions`, `to_parquet`, `interpolate`, `to_tagore`); gnomix-style phasing on haplotype codes |
| 0.6.0 | legacy triple API and `.bin` cache removed; `logging` throughout; package re-layout |
| unreleased | `counts_from_hap_codes` works in input dtype (~4x less temporary memory per dask block) |

The manuscript benchmark will validate and freeze a **v0.7.0 manuscript release**. The benchmark is **not** a gate to v1.0.0.

**Release policy for this revision:**

- `0.7.0`: benchmarked manuscript release after the new correctness suite passes.
- `0.7.x`: bug-fix releases during manuscript preparation if needed.
- `0.8.0` (only if needed): substantive API changes requested during peer review.
- **`1.0.0`: only after peer review is complete / the manuscript is accepted**, after rerunning all correctness and benchmark analyses affected by reviewer-driven changes.

The benchmark and correctness suite should become the long-term regression suite used before the eventual v1.0.0 release, but the paper itself does not claim that v1.0.0 has already been reached.

### 0.2 Benchmark repository

The present benchmarking repository predates the Dataset rewrite. Its analysis scripts target removed APIs such as `read_rfmix`, `binary_dir=`, `generate_binary=`, and `generate_tagore_bed`. They therefore need to be rewritten rather than patched.

The recorded benchmark results are currently `metrics/_h/extracted_metrics.csv`: 19 rows, one observation per condition, two generic parsers, two datasets, and no estimate of dispersion. Real-data inputs are symlinks into `/projects/b1213/resources/processed-data/local-ancestry/` on Quest.

The original preprint was posted in 2024 and evaluated an earlier benchmark implementation. This revision should be treated as a new validation and benchmarking study of the current Dataset architecture.

### 0.3 Scope of the revision

The revised paper should answer three questions:

1. **Correctness:** does the current reader faithfully recover local-ancestry data across supported representations and downstream operations?
2. **Performance:** how do parsing, persistent storage, reopening, querying, and common transformations scale relative to task-appropriate alternatives?
3. **Utility:** what analyses become simpler or more tractable when local ancestry is available as a lazy, haplotype-resolved Dataset rather than repeatedly reparsed flat files?

The paper should not be framed as "a faster parser." Its defensible contribution is a **validated local-ancestry data layer** with persistent storage, random access, haplotype-resolved representation, and a common API across RFMix and FLARE outputs.

---

## 1. Why the present benchmark will not survive review

The existing benchmark has eight major weaknesses.

1. **No replication.** One observation per benchmark cell provides no estimate of run-to-run variability.
   - **Fix:** 10 independent timing replicates per benchmark cell; report median and IQR.

2. **Failures are encoded as measurements.** Existing `cuDF` OOM runs are represented by numeric memory values and missing time.
   - **Fix:** explicit `status` field (`ok`, `oom`, `timeout`, `unsupported`, `error`) and censored display.

3. **Missing modern raw-parser baseline.** A prior Bioinformatics review specifically requested **Polars**.
   - **Fix:** include both **pandas and Polars** as required generic parsing baselines.

4. **Comparator tasks are not always equivalent.** A `convert(..., cache_dir)` call performs parsing, normalization, and persistence, while a raw `read_csv` performs only parsing.
   - **Fix:** split parse, cache construction, reopen, and downstream operations into distinct tasks and compare only equivalent work.

5. **Too few scale points.** Two fixed datasets cannot support a scaling claim.
   - **Fix:** a predeclared scale ladder over samples, markers/segments, and ancestry count.

6. **No independent correctness arm.** Performance claims are insufficient without demonstrating that supported readers and Dataset operations recover known expected values.
   - **Fix:** current-version correctness tests built from independent synthetic fixtures and naive reference implementations. No legacy software version is required.

7. **One workload is benchmarked while the package ships many operations.**
   - **Fix:** benchmark parsing, cache construction, reopen, counts, positional query, interpolation, export, and haplotype-resolved access separately.

8. **Measurement and provenance are uncontrolled.** Page-cache state, thread counts, hardware, and hand-transcribed metrics can all distort the result.
   - **Fix:** one child process per run, pinned threads, homogeneous node class, structured JSONL output, checksums, and automatic aggregation.

---

## 2. Repository restructure — house style

The reference implementation for repository organization is `~/Projects/isograph-brain-aging-benchmarking`.

### 2.1 Target layout

```text
.here
AGENTS.md
ANALYSIS_MAP.md
TODO.md
README.md

configs/
  scale_ladder.yaml
  tools.yaml
  tasks.yaml
  hardware.yaml
  application_loci.yaml
  panels.yaml

env/
  environment.lock.txt
  requirements-analysis.txt

inputs/
  README.md
  raw/
  processed/

rfmix_reader_benchmark/
  __init__.py
  paths.py
  config.py
  harness.py
  manifest.py

  adapters/
    rfmix_reader_07.py
    admix_kit.py
    tractor.py
    pandas_naive.py
    polars_naive.py
    cyvcf2.py

  correctness/
    fixtures.py
    parse_exactness.py
    posterior_checks.py
    cross_format.py
    zarr_roundtrip.py
    ancestry_axis.py
    missing_calls.py
    operations_reference.py
    phase.py

  application/
    ancestry_of_allele.py
    locus_query.py
    diagnosis_summary.py
    design_table.py

  stats/
    summarize.py
    scaling_fit.py

  figures/
    benchmark.R
    application.R

01_simulation_design/
  00_design/{_h,_m}
  01_simulate/{_h,_m}

02_correctness/
  {_h,_m,_m/logs}

03_benchmark/
  01_cpu/{_h,_m}
  02_lowmem/{_h,_m}
  03_realdata/{_h,_m}
  04_metrics/{_h,_m,figures}

04_application/
  01_brainseq/{_h,_m}
  02_public_panels/{_h,_m}
  03_ancestry_of_allele/{_h,_m}
  04_locus_query/{_h,_m}
  05_figures/{_h,_m,figures}

manuscript/
  README.md
  MANUSCRIPT_PLAN.md
  FIGURE_ORDERING.md
  _m/

reports/
scripts/
tests/
zenodo/

archive/v1_preprint/
develop/
```

### 2.2 What is removed from the prior plan

The revised repository does **not** contain:

- a legacy 0.3.1 environment;
- a legacy reader adapter;
- any correctness test against an old software release;
- a GPU benchmark stage.

The historical software versions are not needed to validate the current package and were not part of the original preprint benchmark.

### 2.3 `paths.py`

All analyses write through one registry:

```python
OUTPUT_DIRS: dict[str, tuple[str, ...]] = {
    "sim.design":       ("01_simulation_design", "00_design", "_m"),
    "sim.data":         ("01_simulation_design", "01_simulate", "_m"),
    "correctness":      ("02_correctness", "_m"),
    "bench.cpu":        ("03_benchmark", "01_cpu", "_m"),
    "bench.lowmem":     ("03_benchmark", "02_lowmem", "_m"),
    "bench.real":       ("03_benchmark", "03_realdata", "_m"),
    "bench.metrics":    ("03_benchmark", "04_metrics", "_m"),
    "app.brainseq":     ("04_application", "01_brainseq", "_m"),
    "app.public":       ("04_application", "02_public_panels", "_m"),
    "app.allele":       ("04_application", "03_ancestry_of_allele", "_m"),
    "app.query":        ("04_application", "04_locus_query", "_m"),
    "manuscript":       ("manuscript", "_m"),
    "tmp":              ("03_benchmark", "04_metrics", "_m", "tmp"),
}
```

### 2.4 Benchmark manifest

Do not hard-code a fixed Slurm array size.

`python -m rfmix_reader_benchmark.manifest` should expand:

- scale ladder;
- tool support matrix;
- tasks;
- hardware class;
- 10 replicates;

into a deterministic `run_manifest.parquet`.

Each row is one benchmark run. The Slurm submission helper reads the manifest and submits `--array=0-(N-1)`.

This makes the run count auditable and prevents array indices from silently becoming inconsistent with YAML configuration.

### 2.5 Quest wrapper skeleton

CPU is the primary benchmark environment.

```bash
#!/usr/bin/env bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=rfmixr-bench-cpu
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --constraint=quest10
#SBATCH --output=03_benchmark/01_cpu/_m/logs/bench-cpu-%A_%a.log

set -euo pipefail

PROJECT_ROOT="${RFMIX_READER_BENCHMARK_ROOT:-${PWD}}"
cd "${PROJECT_ROOT}"

[[ -f .here && -d rfmix_reader_benchmark ]] || {
    echo "ERROR: submit from repository root."
    exit 1
}

export RFMIX_READER_BENCHMARK_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

module purge
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/rfmix-bench

python -m rfmix_reader_benchmark.harness \
    --manifest 03_benchmark/01_cpu/_m/run_manifest.parquet \
    --row "${SLURM_ARRAY_TASK_ID}"
```

### 2.6 Repository actions

- Move original preprint scripts and figures to `archive/v1_preprint/`.
- Create `.here` at repository root.
- Rewrite `.gitattributes` to match large artifacts by extension rather than literal path.
- Drop obsolete `.bin` cache patterns from `.gitignore`.
- Add `la_cache/`, `*.zarr/`, `develop/`, `inputs/raw/`, and per-run output directories to `.gitignore`.
- Keep only compact result tables and manuscript-facing figures under version control.
- Stage heavy regenerable data and run-level artifacts for Zenodo.
- Add `tests/` for:
  - JSONL schema;
  - benchmark-manifest determinism;
  - status assignment;
  - timeout/OOM detection;
  - path registry;
  - adapter output schema.
- Reconcile the benchmark repository license with the package's GPL-3.0-or-later license before submission.

---

## 3. Benchmark design

## 3.1 Core principle: compare equivalent work

The benchmark should not combine operations that perform different amounts of work.

Each task is therefore defined by its **semantic output**, not merely by a function name.

A comparator is timed only if it can produce the same or a clearly equivalent result.

Unsupported functionality is recorded in the feature matrix, but an unsupported task is **not treated as a performance loss**.

---

## 3.2 Factors

| Factor | Levels | Rationale |
|---|---|---|
| Source format | RFMix `.msp.tsv`, RFMix `.fb.tsv`, FLARE `.anc.vcf.gz`, haptools `simgenotype` where applicable | covers supported input families |
| Samples | 100, 500, 2,000, 10,000 | pilot through biobank-scale arm |
| Marker/segment scale | chr21, chr1, genome-wide | small, large chromosome, full analysis |
| Ancestries K | 2, 3, 5 | posterior width and hard-call complexity |
| Posteriors | kept / discarded | isolates posterior-memory cost |
| Hardware | standard CPU, low-memory CPU | substantiates both throughput and memory-light claims |
| Replicates | **10** | stable median/IQR and paired ratio estimates |

Do not run a full factorial.

Use a **screening + focused design**:

1. full sample/marker ladder on `.msp.tsv`, K=3, CPU;
2. full sample/marker ladder on `.fb.tsv`, K=3, CPU;
3. K sensitivity at selected medium and large scales;
4. posterior kept/discarded sensitivity;
5. low-memory replay of the most informative large cells;
6. real-data replay on BrainSeq.

This should provide enough points for scaling without spending most compute on redundant combinations.

---

## 3.3 Core comparators

### Required timed comparators

| Tool | Role | Timed where |
|---|---|---|
| `rfmix-reader` 0.7.0 | subject | all supported tasks |
| `admix-kit` | domain-specific local-ancestry data layer | overlapping local-ancestry load/query operations |
| **Polars** | modern high-performance generic parser requested in prior review | equivalent flat-file parse / reshape / aggregation tasks |
| pandas | conventional generic parser baseline | same generic flat-file tasks as Polars |
| `cyvcf2` | optimized VCF parser | FLARE VCF ingestion and selected extraction tasks |
| Tractor | downstream ancestry-aware workflow comparator | only where it performs semantically comparable haplotype/allele extraction |

### Not a core timed comparator

`LAIT` can be discussed as prior art or included in a feature table if useful, but it does not need to appear in every timed benchmark.

### GPU

The GPU backend is optional and is **not part of the manuscript benchmark**.

No dedicated GPU benchmark is planned.

This keeps the main claims aligned with the default CPU execution path that all users can access.

---

## 3.4 Benchmark tasks

### B1. Cold parse to an equivalent in-memory representation

**Question:** how quickly and with how much peak memory can each tool parse the source representation into a usable in-memory representation containing the same essential information?

Examples:

- RFMix MSP → hard ancestry calls/segments;
- RFMix FB → posterior matrix plus hard calls when requested;
- FLARE ANC VCF → ancestry calls.

For pandas and Polars, the baseline includes the minimal deterministic transformation needed to reach the declared canonical table/array representation. Timing a bare `read_csv` against a fully normalized Dataset is not allowed.

### B2. Parse + persistent query-store construction

**Question:** what is the one-time cost to build rfmix-reader's persistent representation?

For rfmix-reader:

```python
convert(path, fmt, cache_dir)
```

This task is reported as a separate operation from B1.

Comparators are included only if they expose a persistent local-ancestry store with sufficiently similar semantics.

### B3. Reopen persistent representation

**Question:** after preprocessing once, what is the cost of reopening data for a new analysis session?

For rfmix-reader:

```python
open_local_ancestry(cache_dir)
```

This is a core architectural claim.

### B4. Genome-wide ancestry counts

**Question:** what is the cost of computing genome-wide/sample-level ancestry counts from the loaded representation?

For rfmix-reader:

```python
ds.la.counts.compute()
```

This directly evaluates the lower-memory integer implementation.

### B5. Random positional query

**Question:** what is the cost of retrieving local ancestry at an arbitrary list of genomic positions?

For rfmix-reader:

```python
ds.la.at_positions(loci_df)
```

Use several query sizes, e.g.:

- 10 loci;
- 100 loci;
- 1,000 loci;
- 10,000 loci.

This should be one of the headline benchmarks because repeated random access is a major reason to persist the Dataset.

### B6. Interpolation to a variant grid

For rfmix-reader:

```python
ds.la.interpolate(variants_df, outdir)
```

Compare only against tools/workflows that produce an equivalent ancestry assignment on the target variant grid.

### B7. Export / downstream handoff

Benchmark equivalent output generation:

- Parquet;
- BED where semantically appropriate.

Do not combine unlike output formats into one timing estimate.

### B8. Haplotype-resolved ancestry / ancestry-of-allele extraction

This task measures the workflow needed to join phased genotype alleles to local-ancestry haplotypes.

It is included because it directly supports the real-data application.

For comparisons such as Tractor, define one fixed output schema and time the full workflow required to produce that output.

---

## 3.5 Correctness arm

The correctness suite validates the current software directly. It does **not** compare against any old rfmix-reader release.

### C1. Hand-constructed reader fixtures

Create small, human-auditable fixtures for each supported source format from a canonical ancestry tensor with:

- known variant positions;
- known sample order;
- known haplotype order;
- known ancestry codes;
- known ancestry-name ordering;
- known posterior values where applicable.

**Pass criterion:** exact recovery of discrete ancestry calls and metadata; floating-point tolerance for posterior probabilities.

### C2. Cross-format semantic equivalence

Encode the **same canonical local-ancestry tensor** into equivalent supported test representations.

Parse every representation and transform back to the canonical schema.

**Pass criterion:** identical hard calls and locus/sample/haplotype order after normalization.

This tests readers rather than differences between RFMix and FLARE inference algorithms.

### C3. Posterior checks

For synthetic `.fb.tsv` fixtures:

- ancestry-axis order;
- posterior normalization;
- argmax ancestry;
- posterior-margin calculation.

**Pass criterion:** exact expected ancestry assignment and bounded numeric tolerance on posterior values.

### C4. Zarr round-trip

```text
source -> convert -> reopen -> canonical representation
```

**Pass criterion:** bitwise identity for integer ancestry calls and agreed float tolerance for posterior arrays.

### C5. Missing-call handling

Inject:

- all-zero posterior rows;
- missing segments;
- explicit missing ancestry codes;
- boundary positions.

**Pass criterion:** missingness is preserved and excluded from counts exactly as documented.

### C6. Downstream operations against a naive reference

For small fixtures, independently implement simple reference versions of:

- ancestry counts;
- `at_positions`;
- interpolation;
- export content checks.

**Pass criterion:** rfmix-reader output exactly matches the naive reference.

### C7. Phase operation

Use simulated haplotypes with known switch errors.

Measure:

- switch-error rate before;
- switch-error rate after;
- whether any correctly phased segments are degraded.

### C8. End-to-end integration smoke test

Run a small haptools → RFMix / FLARE workflow and confirm that parsed outputs are internally consistent and usable in all downstream operations.

This is an integration test, **not** a claim that disagreement between RFMix and FLARE is reader error.

### Release rule

All deterministic correctness checks required for the manuscript must pass before freezing `0.7.0`.

A correctness failure blocks the manuscript benchmark until fixed.

It does **not** trigger or imply a v1.0.0 release.

---

## 3.6 Measurement protocol

### Process isolation

Each timed run is launched in a fresh subprocess.

No two tools share one interpreter process.

### Pairing / blocking

To make "paired on the same node" real rather than nominal:

- one Slurm allocation handles one benchmark block;
- comparator order is randomized within the allocation;
- each comparator is still launched as an independent child process;
- input, node, thread limits, and resource ceilings are held constant.

This avoids confounding tool performance with between-node differences.

### Page-cache state

Cold and warm measurements are separate conditions.

For cold reads, either:

- use a distinct physical copy per replicate; or
- use an explicit cache-eviction strategy that is validated on Quest.

Do not average cold and warm runs together.

### Threading

Record and pin:

```text
OMP_NUM_THREADS
OPENBLAS_NUM_THREADS
MKL_NUM_THREADS
POLARS_MAX_THREADS
Dask scheduler/worker configuration
```

A default benchmark should use one declared CPU concurrency policy for all tools unless a task is explicitly a scaling-by-threads experiment.

### Timing

Use:

- `time.perf_counter()` inside the child;
- `/usr/bin/time -v` outside the child.

Remove line-level `memory_profiler` instrumentation from every timed path.

### Memory

Headline CPU metric:

- `/usr/bin/time -v` maximum resident set size.

Secondary checks:

- `resource.getrusage`;
- persistent cache size on disk as a separate metric.

Never conflate peak RAM with on-disk cache size.

### Failures

Per-cell wall and memory ceilings are defined in configuration.

Each run emits:

```text
ok
oom
timeout
unsupported
error
```

OOM and timeout are treated as censored outcomes, not numeric completion measurements.

### Provenance

Every run records:

- benchmark repository SHA;
- rfmix-reader SHA;
- tool version;
- environment hash;
- input checksum;
- CPU model;
- RAM;
- host;
- thread settings;
- cache state.

---

## 3.7 Statistical reporting

### Per benchmark cell

Report:

- median;
- IQR;
- min/max for diagnostic tables;
- `n=10`;
- completion/failure counts.

Do not report a mean that includes or substitutes for censored runs.

### Paired performance ratios

When two tools complete the same benchmark block, calculate paired ratios within replicate/block.

Headline comparisons should report the distribution of:

```text
competitor wall time / rfmix-reader wall time
competitor peak RSS / rfmix-reader peak RSS
```

### Uncertainty

With 10 replicates, use:

- median + IQR as primary descriptive statistics;
- bootstrap confidence intervals for headline paired median ratios where useful.

Avoid presenting a bootstrap CI as more precise than the underlying 10-run design permits.

### Scaling models

Do **not** force every format into one `samples × loci` predictor.

Use format-appropriate complexity measures.

For dense posterior data such as `.fb.tsv`:

```text
log(time) ~ log(n_samples × n_markers × K)
```

For segment-oriented hard-call representations:

```text
log(time) ~ log(n_samples × n_segments)
```

or an empirically justified equivalent.

Report:

- slope;
- confidence interval;
- plotted raw data.

Use robust regression / robust standard errors if residual behavior warrants it.

The manuscript should describe these as **empirical scaling relationships**, not theoretical complexity proofs.

---

## 3.8 Metrics schema

One JSON object per run:

```text
run_id
timestamp
benchmark_repo_sha
rfmix_reader_sha
tool
tool_version
env_hash

task
source_format
n_samples
n_markers
n_segments
n_ancestries
keep_posteriors

hardware_class
cpu_model
n_threads
host_ram_gb
replicate
block_id
within_block_order
cache_state

status
wall_s
cpu_user_s
cpu_sys_s
peak_rss_mb
cache_bytes_on_disk

input_checksum
stderr_tail
```

`03_benchmark/04_metrics/_h/01.aggregate.py` reads JSONL and writes:

```text
03_benchmark/04_metrics/_m/metrics.parquet
```

Every manuscript benchmark figure reads only the aggregated Parquet or a documented derivative.

No manually transcribed benchmark table remains in the workflow.

---

## 4. Real-data application

## 4.1 Primary application: BrainSeq

Use BrainSeq as the primary real-data application and real-data benchmark because it provides:

- a substantially larger admixed cohort than the relevant 1000 Genomes populations;
- phased/genotyped donors suitable for local-ancestry analysis;
- schizophrenia and control labels;
- a biologically relevant cohort already used by the group.

The manuscript does **not** need to preserve the earlier AD application concept if the available data make a schizophrenia application stronger.

The application should remain a demonstration of what the Dataset architecture enables, not an attempt to reproduce a full schizophrenia GWAS in a modest cohort.

### Core biological/workflow question

> Can phased genotype alleles at schizophrenia-relevant loci be assigned efficiently and reproducibly to their local-ancestry backgrounds, and can the same persistent local-ancestry Dataset support repeated locus-level queries without reparsing the original RFMix output?

This directly exercises the feature that is difficult to represent with counts-only or one-off flat-file workflows.

---

## 4.2 BrainSeq locus set

Create `configs/application_loci.yaml`.

The locus set should be fixed before looking at ancestry-of-allele results.

Potential sources include:

- genome-wide significant schizophrenia lead variants;
- credible-set variants from a current schizophrenia GWAS;
- loci already central to the group's schizophrenia work;
- a compact set of biologically interpretable loci used for figure display.

The analysis should distinguish:

1. **predeclared full locus set** used for all summary statistics;
2. **small display set** used for manuscript figures.

Do not select figure loci based on the most interesting ancestry pattern observed after analysis.

---

## 4.3 BrainSeq analyses

### A1. Cohort-level local-ancestry summary

Using the persistent Dataset:

- local-ancestry proportion per donor;
- missing-call rate;
- ancestry-call coverage;
- number of ancestry transitions / segments as a QC summary.

This establishes that the real-data Dataset is usable before locus-level analysis.

### A2. Ancestry of risk / alternate alleles

At each predeclared schizophrenia locus:

1. retrieve local ancestry with `at_positions`;
2. join local-ancestry haplotypes to phased genotype alleles;
3. count reference and alternate/risk alleles by local-ancestry background.

Output:

```text
locus
variant
allele
local_ancestry
n_allele_copies
n_donors
frequency
confidence_interval
```

This is the primary haplotype-resolved software demonstration.

### A3. Diagnosis-aware descriptive summary

Because BrainSeq has schizophrenia and control status, report the ancestry-of-allele counts stratified by diagnosis.

This is initially descriptive.

Do not imply that a difference is a replicated ancestry-specific schizophrenia association.

### A4. Conditional exploratory association

Only fit a genotype × local-ancestry or ancestry-specific allelic model at a locus if a predeclared minimum-information criterion is met.

For example, require adequate counts of the tested allele across the relevant ancestry backgrounds and diagnosis groups.

If no loci meet the criterion, the manuscript loses nothing: A2/A3 remain the application.

If loci do meet the criterion, association results are explicitly labeled **exploratory** and are not required for the software paper.

### A5. Real-data query workload

Use the same BrainSeq Dataset to benchmark repeated queries:

- 10 loci;
- 100 loci;
- 1,000 loci;
- full predeclared locus set or a large variant grid.

Compare:

1. repeated reparsing workflow;
2. reopen persistent Dataset + `at_positions`.

This closes the loop between the software benchmark and the biological application.

---

## 4.4 Public 1000 Genomes companion analysis

1000 Genomes remains useful because reviewers can reproduce it without restricted BrainSeq access.

However, it is no longer the primary biological dataset.

Use it as a **public reproducibility companion**, for example:

- a small ancestry-of-allele demonstration;
- a public end-to-end tutorial;
- a validation that the same workflow runs outside BrainSeq.

Do not build the manuscript's central biological claim around rare allele counts in small admixed 1000 Genomes populations.

If three-way ancestry inference in Latino populations requires a reference panel beyond unadmixed 1000 Genomes populations, document and freeze the external reference-panel definition in `configs/panels.yaml`.

If that adds disproportionate complexity, narrow the public companion analysis rather than allowing reference-panel construction to become a second paper.

---

## 4.5 Application scope guard

The application should show that the software enables a scientifically useful operation.

It should not become a separate underpowered disease-association manuscript.

### For Bioinformatics Advances

Lead with the operation:

```text
persistent Dataset -> position query -> phased allele/local-ancestry join
```

BrainSeq is the motivating real dataset.

The biological result demonstrates software utility.

### For HGG Advances

Lead with the human-genetics question:

```text
on which local-ancestry backgrounds are schizophrenia-relevant alleles observed in an admixed brain cohort?
```

Then show that the Dataset architecture makes the analysis reproducible and scalable.

The same core analyses therefore support both target journals.

---

## 5. Manuscript display items

## 5.1 Main figures

### Figure 1 — Software architecture and correctness

Possible panels:

- source formats → canonical Dataset → persistent cache → downstream operations;
- C1/C2 exact concordance;
- posterior / ancestry-axis checks;
- Zarr round-trip;
- phase simulation.

This figure establishes trust before speed.

### Figure 2 — Scaling and memory

CPU wall time and peak RSS across the sample/data-size ladder.

Include:

- rfmix-reader;
- admix-kit where equivalent;
- Polars;
- pandas;
- cyvcf2 where relevant.

Plot raw points + median/IQR.

Show censored OOM/timeout outcomes explicitly.

### Figure 3 — Parse-once, query-many

Focus on the architectural advantage:

- cache construction;
- reopen;
- repeated positional queries at increasing query-set size;
- crossover point relative to repeated flat-file parsing.

This should be one of the central manuscript figures.

### Figure 4 — Task/feature matrix

Rows:

- B1–B8.

Columns:

- rfmix-reader;
- admix-kit;
- Polars;
- pandas;
- cyvcf2;
- Tractor where relevant.

Separate **feature support** from **performance**.

Do not encode unsupported functionality as an arbitrarily poor performance value.

### Figure 5 — BrainSeq ancestry-of-allele application

Possible panels:

- cohort ancestry summary;
- ancestry background of alleles at the predeclared locus set;
- diagnosis-stratified descriptive counts for selected display loci;
- real-data query performance.

If HGG Advances becomes the target, this figure can move earlier in the manuscript.

---

## 5.2 Tables

### Table 1 — Feature/comparator matrix

Formats × operations × tools.

### Table 2 — Benchmark design and completion status

Scale cells, replicates, failures, and censoring.

### Table 3 — Correctness suite

Check, fixture, expected output, pass/fail.

### Table 4 — BrainSeq locus summary

For every predeclared locus:

- variant;
- allele;
- local-ancestry counts;
- diagnosis-stratified counts;
- whether exploratory association criterion was met.

---

## 5.3 Supplement

Include:

- full `metrics.parquet`;
- run manifest;
- correctness report;
- environment lockfile;
- benchmark configuration YAML;
- all predeclared application loci;
- public-panel workflow if retained;
- Zenodo archive.

Restricted BrainSeq genotype/local-ancestry data are described through the appropriate controlled-access mechanism rather than redistributed.

---

## 6. Eight-week execution timeline

| Week | Work | Gate |
|---|---|---|
| 1 | restructure repository; implement manifest/harness; adapters for rfmix-reader, pandas, **Polars**; define canonical outputs | toy benchmark produces valid JSONL and deterministic manifest |
| 2 | implement C1–C8 correctness suite; add admix-kit/cyvcf2/Tractor adapters where needed | **all required correctness checks pass → freeze/tag v0.7.0 manuscript release** |
| 3 | generate scale-ladder simulation inputs; stage BrainSeq real-data inputs; freeze application locus list | checksums and configs frozen |
| 4 | CPU benchmark blocks, 10 repeats per cell | primary timing/memory dataset complete |
| 5 | low-memory replay + BrainSeq real-data benchmark + application A1–A5 | all runs complete or explicitly censored |
| 6 | aggregate statistics; scaling models; figures/tables | manuscript result set frozen |
| 7 | manuscript draft + reproducibility documentation + Zenodo staging | full draft circulated |
| 8 | revision, journal-specific framing, cover letter, submission | **manuscript submitted** |

### Post-submission / peer review

- reviewer-requested bug fixes: `0.7.x`;
- substantive API changes if required: `0.8.0`;
- rerun every affected correctness and benchmark cell;
- **release v1.0.0 only after peer review is complete / manuscript acceptance**.

---

## 7. Journal strategy

The analysis plan should remain common until results are complete.

## 7.1 Bioinformatics Advances framing

Best fit if the strongest result is:

> rfmix-reader provides a validated, format-independent, persistent local-ancestry Dataset that substantially reduces the cost of repeated downstream queries and provides a unified API for analyses that otherwise require bespoke parsing and transformation.

Suggested Results order:

1. correctness;
2. architecture and supported operations;
3. scaling/memory;
4. parse-once/query-many benchmark;
5. BrainSeq application.

Keep the application explicitly tied to software operations.

Polars should be visible in the main benchmark because its absence was raised in prior Bioinformatics review.

## 7.2 HGG Advances framing

Best fit if BrainSeq produces a stronger genetics/resource result.

Suggested Results order:

1. local-ancestry analysis problem in an admixed human brain cohort;
2. BrainSeq ancestry-of-allele result;
3. correctness of the data layer;
4. performance/scaling;
5. implications for ancestry-aware human-genetic study design.

The software remains central, but the manuscript opens with the genetics problem rather than the engineering abstraction.

## 7.3 Decision rule

Choose the journal after Figure 5 is available.

**Prefer Bioinformatics Advances if:**

- correctness + query architecture are the strongest results;
- BrainSeq mainly demonstrates the API;
- biological findings are descriptive.

**Prefer HGG Advances if:**

- the BrainSeq ancestry-of-allele analysis yields a coherent, interpretable human-genetics result;
- the case/control or ancestry-background summaries provide information beyond software demonstration;
- the genetics story can lead without overstating power.

Society membership does not need to influence the scientific choice because the user is a member of the relevant society for either target.

---

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Correctness fixtures expose a current-reader bug | This is the suite working. Fix before benchmark freeze and rerun affected checks. |
| Polars substantially outperforms rfmix-reader on cold parse | Expected and acceptable if rfmix-reader's value is persistent reuse, queryability, memory behavior, or functionality. Report the crossover rather than hiding it. |
| rfmix-reader loses B1 at small scale | Expected for a richer representation. The main claim is not universal cold-parse dominance. |
| Domain comparator cannot perform a task with equivalent semantics | Mark unsupported in the feature table; do not invent an unfair timing comparison. |
| BrainSeq has insufficient allele counts for diagnosis-aware association | Keep A2/A3 descriptive. Exploratory association is conditional and nonessential. |
| BrainSeq restricted access weakens reproducibility | Retain a small public 1000 Genomes companion workflow and publish all code/configuration. |
| 1000 Genomes reference-panel construction becomes complicated | Narrow the public companion instead of expanding it into a separate reference-panel project. |
| Quest node heterogeneity affects timings | Run paired comparator blocks within the same allocation and record hardware per run. |
| Ten replicates increase compute burden | Use screening + focused design; avoid a full factorial. |
| Eight-week schedule slips | Correctness + v0.7.0 freeze + primary CPU benchmark have standalone value. Drop nonessential public-panel expansion or exploratory association before compromising the core benchmark. |
| Reviewer asks why not use Polars directly | Main figures explicitly show Polars for equivalent flat-file tasks and distinguish raw parsing from persistent local-ancestry operations. |
| Reviewer asks why not use admix-kit | Benchmark overlapping operations directly and define the contribution in terms of supported representations, persistence, random access, and end-to-end workflow rather than claiming no alternative local-ancestry data layer exists. |

---

## 9. Manuscript north star

The repository `AGENTS.md` should carry this statement:

> **RFMix-reader is a validated data layer for local-ancestry output, not a local-ancestry inference method and not a GWAS framework.**
>
> The benchmark must distinguish raw parsing from the richer work of constructing a persistent local-ancestry representation. It is acceptable for pandas or Polars to win a small cold-read benchmark. The manuscript's central test is whether rfmix-reader faithfully represents local ancestry and reduces the computational and implementation burden of repeated ancestry-aware analyses through persistent, random-access, haplotype-resolved data.
>
> All benchmark comparisons must perform equivalent work. Unsupported functionality is reported as unsupported rather than converted into a performance claim.
>
> The manuscript release is v0.7.0. **v1.0.0 is reserved for the post-peer-review stable release.**
