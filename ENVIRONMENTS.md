# Influenza ONT Pipeline - Environment Setup

This pipeline requires three separate conda environments to avoid dependency conflicts.

## Environments

### 1. `pod5_env` - Core bioinformatics tools
Contains: pod5, minimap2, samtools, bcftools, bedtools

Note: `dorado` is not currently available from the accessible conda channels and must be installed separately if your pipeline uses ONT basecalling.

If `dorado` is installed outside of Conda, pass its path with `-d /path/to/dorado` when running `run_pipeline.sh`, or ensure it is on `PATH`.

**Used in:**
- POD5 merge
- STEP 1: Basecalling
- STEP 2: Demultiplexing
- STEP 4: Mapping
- STEP 5: Personalized references
- STEP 6: Remapping
- STEP 7: Variant calling
- STEP 8: Zero-coverage detection
- STEP 9: Consensus generation
- STEP 10: Multi-fasta aggregation

### 2. `nanoplot_env` - Quality control visualization
Contains: nanoplot

**Used in:**
- STEP 3: NanoPlot QC reports

### 3. `nextclade_env` - Phylogenetic analysis
Contains: nextclade

**Used in:**
- STEP 11: Nextclade processing

## Setup Instructions

### Create all three environments at once:

```bash
cd /mnt/d/influenza_LDSP

# Create pod5 environment
mamba env create -f pod5_env.yml --channel-priority flexible

# Create nanoplot environment
conda env create -f nanoplot_env.yml

# Create nextclade environment
conda env create -f nextclade_env.yml
```

### Or create individually:

```bash
# Core bioinformatics tools
mamba env create -f pod5_env.yml -n pod5_env --channel-priority flexible

# Quality control
conda env create -f nanoplot_env.yml -n nanoplot_env

# Phylogenetics
conda env create -f nextclade_env.yml -n nextclade_env
```

## Pipeline Execution

The pipeline automatically handles environment switching:

```bash
./run_pipeline.sh -i 23042026_INFLU -r consensus_irma.fasta
```

The script will:
1. Activate `pod5_env` at startup
2. Switch to `nanoplot_env` for STEP 3
3. Return to `pod5_env` for STEPS 4-10
4. Switch to `nextclade_env` for STEP 11
5. Deactivate all environments at completion

## Troubleshooting

### If a tool is not found in an environment:

Check which environment is active:
```bash
conda info --envs
```

Verify tool installation:
```bash
conda activate pod5_env
which minimap2

conda activate nanoplot_env
which NanoPlot

conda activate nextclade_env
which nextclade
```

### To update an environment:

```bash
conda env update -f pod5_env.yml --prune
```

### To remove an environment:

```bash
conda remove --name pod5_env --all
conda remove --name nanoplot_env --all
conda remove --name nextclade_env --all
```

## Version Information

- **pod5**: Latest from nanoporetech channel
- **dorado**: Not available from the accessible conda channels; install manually using ONT's official distribution or another supported installer
- **minimap2**: Latest from bioconda
- **samtools/bcftools/bedtools**: Latest from bioconda
- **nanoplot**: Latest from bioconda
- **nextclade**: Latest from bioconda
