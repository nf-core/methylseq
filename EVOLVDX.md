# EvolvDx additions to nf-core/methylseq

This fork tracks `nf-core/methylseq` (`upstream` remote). All EvolvDx changes live on the
`evolvdx-main` branch; `master` mirrors nf-core. Keep upstream files untouched where possible
so that new releases merge cleanly.

## Running

```bash
nextflow run EVolvDx/methylseq -r evolvdx-main -profile evolvdx,docker --input samplesheet.csv --outdir results
```

| Profile | Purpose |
|---|---|
| `evolvdx` | Production defaults: `bwameth`, `grch38_core_bs_controls`, `clip_r1 20`, `clip_r2 15`, EvolvDx QC steps on |
| `evolvdx_test` | Small NSQCAM529 subsample from S3, EvolvDx steps and Picard HS report on |
| `evolvdx_aws`, `evolvdx_aws_500gb` | AWS Batch queues (`conf/evolvdx/`) |

Plain nf-core behaviour is unchanged unless one of the flags below is set.

## What was added

| Flag | Component | Notes |
|---|---|---|
| `--run_methylqc` | `subworkflows/local/methylqc` | Cohort-level methylKit/HDF5 merge, control stats, read stats, tiled stats, Twist and repeat annotation, per-chromosome histograms, all summarised in the MultiQC report ("EvolvDx methylation QC") with the raw tables and plots in `methqc/`. The samplesheet's optional `group` column sets the methylKit treatment groups. Needs `--aligner bwameth`; runs a second MethylDackel extraction for the methylKit files (do not also pass `--methyl_kit`, which would suppress the bedGraph output). |
| `--run_cpg_coverage` | `subworkflows/local/cpg_cov` | CpG island intersect and per-base coverage. Needs `--cpg_island_bed`. |
| `--run_methsnsv` | `modules/local/exo/methsnsv.nf` | Methylation SNV/SV statistics from MethylDackel bedGraphs. |
| `--run_picardhs_report` | `modules/local/exo/custom_multiqc` | EvolvDx Picard HS tables in MultiQC. Needs `--run_targeted_sequencing --collecthsmetrics`. |

Also: `conf/evolvdx/custom_genomes.config` (S3 genomes under `--exogenomes_base`), `bin/` analysis
scripts, `assets/multiqc/` headers, `env/` Dockerfiles, `analysis/notebooks/`, `scripts/`, and
`tests/samplesheets`, `tests/scripts`.

## Deliberately not carried over from the old fork

- `PICARD_HS` subworkflow and local Picard modules: replaced by upstream `TARGETED_SEQUENCING`
  (`--run_targeted_sequencing --target_regions_file ... --collecthsmetrics`). Upstream uses one
  regions file for both bait and target; the old fork had separate `probe_bed` and `target_bed`.
- `SUBREAD_FEATURECOUNTS` (gene counts on the methylation BAMs).
- `cegx` and `epignome` library presets (use explicit `--clip_r1/--clip_r2/...`).
- `with_hybcap`, `probe_bed`, `target_bed`, `fai` params, and the `max_cpus/max_memory/max_time` params.

## Known upstream issues (present in nf-core/methylseq 4.2.0 itself)

`nextflow lint` reports errors in `modules/nf-core/rastair/mbiasparser/main.nf`, the `aws_batch`
profile include (`conf/aws/batch/nextflow.config` is missing) and `nf-test.config`.
