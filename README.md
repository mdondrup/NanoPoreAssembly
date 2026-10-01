# L. salmonis Snakemake scaffold

Raw sequencing data lives in `Lsalmonis_raw_data/`. Runs have different nesting depths; future rules should discover inputs beneath this root rather than assume a fixed run path. Each run contains `fastq_pass/` and `fastq_fail/` directories, with POD5 and sequencing reports alongside them.

`Snakefile` loads `config/config.yaml`. The first rule basecalls the `FK2-AdFe-j13-LSK*` sample only, creating a separate BAM for `pod5/` and `pod5_skip/` in each run under `results/basecalling/`. Standard error goes to `logs/basecalling/`; raw data is never modified. Keeping runs and signal sources separate leaves merging and later analysis for subsequent steps.

Set `dorado_model` in the config to a Dorado model name or local model path before running, or override it on the command line:

```sh
snakemake -n --config dorado_model=YOUR_MODEL
snakemake --cores 1 --config dorado_model=YOUR_MODEL
```

Dorado must be on `PATH` and CUDA GPUs must be available. Use a single concurrent job for now because each call requests `cuda:all`. Add later rules in `workflow/rules/`, scripts in `workflow/scripts/`, and environment definitions in `workflow/envs/`.

Currently the ten matching `pod5/` and `pod5_skip/` directories in this workspace are empty. The dry-run will fail with a missing POD5 input until those files are available in the configured raw-data directory.