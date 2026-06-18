# Influenza ONT LDSP Pipeline

This repository contains a Bash pipeline for processing Oxford Nanopore influenza sequencing data from POD5 input through basecalling, demultiplexing, alignment, variant calling, consensus generation, and optional Nextclade analysis.

## Requirements

- Linux environment
- Conda / Mamba
- Miniconda or Anaconda installed
- A reference FASTA file for your target influenza genome(s)
- POD5 files for the sequencing run

## Environment setup

Run the installer once to create the required environments:

```bash
bash install_pipeline.sh
```

The installer creates and validates:
- `pod5_env` for core bioinformatics tools
- `nanoplot_clean` for NanoPlot QC reports
- `nextclade_env` for Nextclade analysis

## Running the pipeline

Use the pipeline script with your run directory and reference FASTA:

```bash
bash run_pipeline.sh -i 22052026_FLU -r reference.fa
```

### Common arguments

- `-i RUN_DIR`  : sequencing run directory (for example `22052026_FLU`)
- `-r REF_FASTA`: reference FASTA file
- `-o OUTPUT_DIR`: optional output directory (defaults to `ont_pipeline_results`)
- `-m MODEL`: optional Dorado model override
- `-k KIT`: optional kit name override
- `-t THREADS`: optional thread count override
- `-c COVERAGE`: optional coverage threshold for reference selection
- `-d DORADO_BIN`: optional Dorado executable path
- `-f`: force rerun all steps

## Expected input layout

The pipeline expects a run directory containing a `no_sample_id` folder with sequencing outputs, for example:

```text
22052026_FLU/
  no_sample_id/
    .../pod5/
    .../fastq/
```

## Outputs

By default, outputs are placed under:
- `ont_pipeline_results/` for main pipeline outputs
- `nanoplot_results/` for QC reports
- `nextclade_results/` for Nextclade outputs

## Notes

- The pipeline uses checkpoint files under `ont_pipeline_results/.checkpoints/` so reruns can resume from completed steps.
- If you need to restart from scratch, use `-f`.
- Large generated data files are ignored by the repository `.gitignore` file.

## Troubleshooting

- If a tool is missing, rerun `bash install_pipeline.sh`.
- If basecalling or demultiplexing fails, confirm that the POD5 files are present and readable.
- If you see an environment mismatch, verify that the correct Conda environment is being activated for the step you are running.
