#!/usr/bin/env bash

set -euo pipefail

echo "========================================="
echo "Installing Influenza ONT Pipeline"
echo "========================================="

# ----------------------------------------
# Load conda
# ----------------------------------------

if [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
else
    echo "ERROR: conda.sh not found."
    echo "Check your Miniconda installation."
    exit 1
fi

# ----------------------------------------
# Check mamba
# ----------------------------------------

if ! command -v mamba &> /dev/null; then
    echo "Mamba not found. Installing into base environment..."

    conda activate base
    conda install -y -n base -c conda-forge mamba
fi

# ----------------------------------------
# Create required environments
# ----------------------------------------

if conda env list | grep -q "^pod5_env"; then
    echo "pod5_env already exists. Skipping."
else
    echo "Creating pod5_env..."
    mamba env create -f pod5_env.yml
fi

if conda env list | grep -q "^nanoplot_clean"; then
    echo "nanoplot_clean already exists. Skipping."
else
    echo "Creating nanoplot_clean..."
    mamba env create -f nanoplot_env.yml
fi

if conda env list | grep -q "^nextclade_env"; then
    echo "nextclade_env already exists. Skipping."
else
    echo "Creating nextclade_env..."
    mamba env create -f nextclade_env.yml
fi

# ----------------------------------------
# Validate pod5_env
# ----------------------------------------

echo ""
echo "Validating pod5_env..."

conda activate pod5_env

for tool in minimap2 samtools bcftools bedtools seqkit pod5; do
    if ! command -v $tool &> /dev/null; then
        echo "ERROR: $tool missing in pod5_env"
        exit 1
    fi
done

echo "pod5_env OK"

# ----------------------------------------
# Validate nanoplot_clean
# ----------------------------------------

echo ""
echo "Validating nanoplot_clean..."

conda activate nanoplot_clean

if ! command -v NanoPlot &> /dev/null; then
    echo "ERROR: NanoPlot missing in nanoplot_clean"
    exit 1
fi

echo "nanoplot_clean OK"

# ----------------------------------------
# Validate nextclade_env
# ----------------------------------------

echo ""
echo "Validating nextclade_env..."

conda activate nextclade_env

if ! command -v nextclade &> /dev/null; then
    echo "ERROR: nextclade missing in nextclade_env"
    exit 1
fi

echo "nextclade_env OK"

# ----------------------------------------
# Finish
# ----------------------------------------

conda deactivate

echo ""
echo "========================================="
echo "All environments installed successfully."
echo "Pipeline ready to use."
echo "========================================="