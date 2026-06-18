#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

########################################
# DEFAULT CONFIGURATION
########################################

BASECALL_MODEL="dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
KIT_NAME="SQK-NBD114-96"
THREADS=8
COVERAGE_THRESHOLD=50
DORADO_BIN="dorado"

########################################
# HELP
########################################

usage() {
    cat << EOF

ONT Influenza Pipeline

Usage:
    $(basename "$0") -i RUN_DIRECTORY -r REFERENCE_FASTA [options]

Required arguments:
    -i    Sequencing run directory
          Example: 22052026_FLU

    -r    Reference FASTA

Optional arguments:
    -o    Output directory
          Default: ont_pipeline_results

    -m    Dorado model

    -k    Kit name

    -t    Threads
          Default: 8

    -c    Coverage threshold for reference selection
          Default: 50

    -d    Dorado binary or executable path
          Default: dorado (must be installed separately if not in Conda)

    -f    Force restart (skip checkpoints, rerun all steps)

    -h    Show help message

Example:
    ./run_pipeline.sh \
        -i 22052026_FLU \
        -r reference.fa

    # Resume pipeline after it broke at step 7
    ./run_pipeline.sh \
        -i 22052026_FLU \
        -r reference.fa

    # Force complete restart (ignore checkpoints)
    ./run_pipeline.sh \
        -i 22052026_FLU \
        -r reference.fa \
        -f

EOF
}

########################################
# ARGUMENTS
########################################

OUTPUT_DIR="ont_pipeline_results"
FORCE_RESTART=0

while getopts "i:r:o:m:k:t:c:d:fh" opt; do
    case $opt in
        i) RUN_DIR="$OPTARG" ;;
        r) REFERENCE_FASTA="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        m) BASECALL_MODEL="$OPTARG" ;;
        k) KIT_NAME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        c) COVERAGE_THRESHOLD="$OPTARG" ;;
        d) DORADO_BIN="$OPTARG" ;;
        f)
            FORCE_RESTART=1
            echo "Force restart enabled - will rerun all steps"
            ;;
        h)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

########################################
# VALIDATION
########################################

if [[ -z "${RUN_DIR:-}" ]]; then
    echo "ERROR: Sequencing run directory required."
    usage
    exit 1
fi

if [[ -z "${REFERENCE_FASTA:-}" ]]; then
    echo "ERROR: Reference FASTA required."
    usage
    exit 1
fi

########################################
# ACTIVATE CONDA ENVIRONMENT (EARLY)
########################################

echo "Loading conda shims"
source ~/miniconda3/etc/profile.d/conda.sh

check_command() {
    local cmd="$1"
    if [[ "$cmd" == */* ]]; then
        [[ -x "$cmd" ]] || {
            echo "ERROR: $cmd is not executable or not found."
            exit 1
        }
    else
        command -v "$cmd" >/dev/null 2>&1 || {
            echo "ERROR: $cmd is not installed or not in PATH."
            exit 1
        }
    fi
}

echo "Checking dependencies..."
echo "Note: pod5 CLI not required for file collection (we'll move files)."
check_command "$DORADO_BIN"
check_command minimap2
check_command samtools
check_command bcftools
# bedtools is installed in pod5_env; validate later when pod5_env is active

########################################
# DEFINE PATHS
########################################

RAW_DIR="${RUN_DIR}/no_sample_id"

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: no_sample_id directory not found:"
    echo "$RAW_DIR"
    exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found:"
    echo "$REFERENCE_FASTA"
    exit 1
fi

########################################
# OUTPUT STRUCTURE
########################################

BASECALL_OUTPUT="${OUTPUT_DIR}/reads.bam"
DEMUX_DIR="${OUTPUT_DIR}/fastq"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
CONSENSUS_DIR="${OUTPUT_DIR}/consensus"
POD5_DIR="${OUTPUT_DIR}/pod5"
MERGED_POD5=""
LOG_FILE="${OUTPUT_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"

mkdir -p \
    "$OUTPUT_DIR" \
    "$DEMUX_DIR" \
    "$ALIGN_DIR" \
    "$CONSENSUS_DIR"

# Redirect output to log file
exec &> >(tee -a "$LOG_FILE")

########################################
# CHECKPOINTING
########################################

CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "$CHECKPOINT_DIR"

check_step() {
    local step_name="$1"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_name}.done"
    
    if [ $FORCE_RESTART -eq 1 ]; then
        return 1  # Always run if force restart is enabled
    fi
    
    if [ -f "$checkpoint_file" ]; then
        return 0  # Step already completed
    fi
    return 1  # Step not completed yet
}

mark_step_complete() {
    local step_name="$1"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_name}.done"
    touch "$checkpoint_file"
}

skip_step() {
    local step_name="$1"
    echo "[CHECKPOINT] $step_name already completed. Skipping..."
}

echo "Checkpoint directory: $CHECKPOINT_DIR"

########################################
# FIND POD5 FILES
########################################

echo "Searching for POD5 files in: $RAW_DIR"
echo "(searching recursively in all subdirectories)"

mapfile -t POD5_FILES < <(
    find "$RAW_DIR" -type f -name "*.pod5" 2>/dev/null | sort
)

if [[ ${#POD5_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No POD5 files found in:"
    echo "$RAW_DIR"
    echo ""
    echo "Please ensure POD5 files exist in subdirectories like:"
    echo "$RAW_DIR/*/pod5/"
    exit 1
fi

echo "Found ${#POD5_FILES[@]} POD5 files from different directories:"
for i in "${!POD5_FILES[@]}"; do
    if [ $((i)) -lt 5 ]; then
        echo "  - ${POD5_FILES[$i]}"
    elif [ $((i)) -eq 5 ]; then
        echo "  ... and $((${#POD5_FILES[@]} - 5)) more"
        break
    fi
done
echo "Log file: $LOG_FILE"

########################################
# MERGE POD5
########################################

echo "================================="
echo "COLLECTING POD5 FILES INTO COMMON DIRECTORY"
echo "================================="

if check_step "collect_pod5"; then
    skip_step "Collect POD5"
else
    # Check whether all pod5 files already live in a single directory
    mapfile -t PARENT_DIRS < <(
        printf "%s
" "${POD5_FILES[@]}" | xargs -n1 dirname | sort -u
    )

    if [[ ${#PARENT_DIRS[@]} -eq 1 ]]; then
        POD5_DIR="${PARENT_DIRS[0]}"
        echo "All pod5 files already in single directory: $POD5_DIR"
        echo "Skipping move; using existing directory as common pod5 directory."
    else
        mkdir -p "$POD5_DIR"
        echo "Moving ${#POD5_FILES[@]} pod5 files to $POD5_DIR"
        for p in "${POD5_FILES[@]}"; do
            dest="$POD5_DIR/$(basename "$p")"
            if [ "$p" != "$dest" ]; then
                mv "$p" "$dest"
            else
                echo "Already in destination: $p"
            fi
        done
    fi

    mark_step_complete "collect_pod5"
fi

########################################
# BASECALLING
########################################

echo "================================="
echo "STEP 1: BASECALLING"
echo "================================="

if check_step "step_1_basecalling" || [ -s "$BASECALL_OUTPUT" ]; then
    if [ -s "$BASECALL_OUTPUT" ] && [ ! -f "${CHECKPOINT_DIR}/step_1_basecalling.done" ]; then
        mark_step_complete "step_1_basecalling"
    fi
    skip_step "Step 1: Basecalling"
else
    # Ensure basecalling runs in conda base (not in a custom env)
    if command -v conda >/dev/null 2>&1; then
        conda deactivate || true
        conda activate base || true
    fi

    # Ensure the common POD5 directory is valid for Dorado
    if [[ ! -d "$POD5_DIR" ]]; then
        echo "ERROR: POD5 directory not found: $POD5_DIR"
        exit 1
    fi

    # Run dorado basecaller using the POD5 directory path
    echo "Running dorado basecaller with verbose logging..."
    "$DORADO_BIN" basecaller \
        --verbose \
        "$BASECALL_MODEL" \
        "$POD5_DIR" \
        --kit-name "$KIT_NAME" \
        --barcode-both-ends \
        --recursive \
        > "$BASECALL_OUTPUT"

    if [ ! -s "$BASECALL_OUTPUT" ]; then
        echo "ERROR: Dorado basecalling did not produce output: $BASECALL_OUTPUT"
        exit 1
    fi

    mark_step_complete "step_1_basecalling"
fi

########################################
# DEMULTIPLEXING
########################################

echo "================================="
echo "STEP 2: DEMULTIPLEXING"
echo "================================="

if check_step "step_2_demultiplexing"; then
    skip_step "Step 2: Demultiplexing"
else
    "$DORADO_BIN" demux \
        -o "$DEMUX_DIR" \
        --kit-name "$KIT_NAME" \
        --emit-fastq \
        "$BASECALL_OUTPUT"
    mark_step_complete "step_2_demultiplexing"
fi

echo "================================="
echo "STEP 2.5: CONSOLIDATE DEMULTIPLEXED FASTQ"
echo "================================="

if check_step "merge_demux_fastq"; then
    skip_step "Step 2.5: Consolidate demultiplexed FASTQ"
else
    FASTQ_PASS_DIRS=$(find "$DEMUX_DIR" -type d -name "fastq_pass" 2>/dev/null)

    if [ -z "$FASTQ_PASS_DIRS" ]; then
        echo "Warning: No fastq_pass directories found in $DEMUX_DIR. Skipping FASTQ consolidation."
    else
        # Create consolidated directory for FASTQ outputs
        CONSOLIDATED_FASTQ_DIR="${OUTPUT_DIR}/fastq_consolidated"
        mkdir -p "$CONSOLIDATED_FASTQ_DIR"

        # Build a barcode -> FASTQ file mapping from every demux output
        declare -A BARCODE_FASTQ_FILES=()
        while IFS= read -r fastq_file; do
            barcode_dir=$(dirname "$fastq_file")
            barcode=$(basename "$barcode_dir")

            if [[ "$barcode" =~ ^barcode[0-9]+$ ]]; then
                BARCODE_FASTQ_FILES["$barcode"]+="$fastq_file"$'\n'
            fi
        done < <(find "$DEMUX_DIR" -type f -name "*.fastq" | sort)

        for barcode in $(printf '%s\n' "${!BARCODE_FASTQ_FILES[@]}" | sort); do
            mapfile -t ALL_FASTQ_FILES < <(printf '%b' "${BARCODE_FASTQ_FILES[$barcode]}" | sort -u)

            if [[ ${#ALL_FASTQ_FILES[@]} -eq 0 ]]; then
                continue
            fi

            consolidated_fastq="${CONSOLIDATED_FASTQ_DIR}/${barcode}_consolidated.fastq"
            echo "Consolidating ${#ALL_FASTQ_FILES[@]} FASTQ files for $barcode"
            cat "${ALL_FASTQ_FILES[@]}" > "$consolidated_fastq"

            if [[ -s "$consolidated_fastq" ]]; then
                echo "  ✓ Created: $consolidated_fastq"
            else
                echo "  ✗ ERROR: consolidated FASTQ is empty for $barcode"
                exit 1
            fi
        done

        echo "All FASTQ files consolidated in: $CONSOLIDATED_FASTQ_DIR"
    fi

    mark_step_complete "merge_demux_fastq"
fi

echo "================================="
echo "STEP 3: NanoPlot QC reports"
echo "================================="

if check_step "step_3_nanoplot"; then
    skip_step "Step 3: NanoPlot QC reports"
else
    # Switch to nanoplot environment
    conda deactivate
    conda activate nanoplot_clean
    check_command NanoPlot
    
    NANOPLOT_DIR="nanoplot_results"
    mkdir -p "$NANOPLOT_DIR"

    # Find all fastq_pass directories within the demux output (handles nested structures)
    FASTQ_PASS_DIRS=$(find "$DEMUX_DIR" -type d -name "fastq_pass" 2>/dev/null)

    if [ -z "$FASTQ_PASS_DIRS" ]; then
        echo "Warning: No fastq_pass directories found in $DEMUX_DIR. Skipping NanoPlot."
    else
        # Process each fastq_pass directory found
        while IFS= read -r FASTQ_PASS; do
            echo "Processing barcodes in: $FASTQ_PASS"
            
            # Loop through each barcode directory
            for barcode_dir in "${FASTQ_PASS}"/barcode*/; do
                [ -d "$barcode_dir" ] || continue
                barcode=$(basename "$barcode_dir")
                
                echo "Processing NanoPlot for $barcode"
                
                # Concatenate all fastq files in this barcode directory
                combined_fastq="${barcode}_combined.fastq"
                cat "${barcode_dir}"/*.fastq > "$combined_fastq" 2>/dev/null || true
                
                if [ ! -s "$combined_fastq" ]; then
                    echo "✗ No fastq data found for $barcode"
                    rm -f "$combined_fastq"
                    continue
                fi
                
                # Create output directory for this barcode
                output_dir="$NANOPLOT_DIR/${barcode}"
                mkdir -p "$output_dir"
                
                # Run NanoPlot on the combined fastq
                NanoPlot --fastq "$combined_fastq" --outdir "$output_dir" --plots kde hex dot
                
                if [ $? -eq 0 ]; then
                    echo "✓ Completed NanoPlot for $barcode"
                else
                    echo "✗ Error running NanoPlot for $barcode"
                fi
                
                # Clean up temporary combined fastq
                rm -f "$combined_fastq"
            done
        done <<< "$FASTQ_PASS_DIRS"
        
        echo "All NanoPlot runs completed. Results in $NANOPLOT_DIR/"
    fi

    mark_step_complete "step_3_nanoplot"
fi

# Switch back to pod5 environment for remaining steps (STEPS 4-10)
conda deactivate
conda activate pod5_env
check_command bedtools

echo "================================="
echo "STEP 4: Mapping FASTQ files"
echo "================================="
find "$DEMUX_DIR" -name "*.fastq" | sort | while read -r fq; do
    base=$(basename "${fq%.fastq}")
    checkpoint_file="${CHECKPOINT_DIR}/step_4_map_${base}.done"
    
    if [ $FORCE_RESTART -eq 0 ] && [ -f "$checkpoint_file" ]; then
        echo "[CHECKPOINT] Already mapped $base. Skipping..."
        continue
    fi
    
    echo "Processing $base"
    minimap2 -t "$THREADS" -ax map-ont "$REFERENCE_FASTA" "$fq" |
        samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${base}.sorted.bam"
    samtools index "$ALIGN_DIR/${base}.sorted.bam"
    touch "$checkpoint_file"
done

mark_step_complete "step_4_mapping"
echo "================================="
echo "STEP 5: Create personalized references per barcode"
echo "================================="

if check_step "step_5_personalized_refs"; then
    skip_step "Step 5: Personalized references"
else
    REFS_DIR="refs_by_barcode"
    mkdir -p "$REFS_DIR"

echo "Using coverage threshold: $COVERAGE_THRESHOLD reads"

# Load reference FASTA into associative arrays (contig -> sequence)
declare -A ref_contigs
declare -a contig_order
current_contig=""
current_seq=""

while IFS= read -r line; do
    if [[ "$line" =~ ^">"(.*)$ ]]; then
        # Save previous contig if it exists
        if [ -n "$current_contig" ]; then
            ref_contigs["$current_contig"]="$current_seq"
            contig_order+=("$current_contig")
        fi
        current_contig="${BASH_REMATCH[1]}"
        current_seq=""
    else
        current_seq+="$line"
    fi
done < "$REFERENCE_FASTA"

# Don't forget the last contig
if [ -n "$current_contig" ]; then
    ref_contigs["$current_contig"]="$current_seq"
    contig_order+=("$current_contig")
fi

echo "Loaded ${#ref_contigs[@]} contigs from reference FASTA"

# Process each BAM file to create per-barcode references
for bam in "$ALIGN_DIR"/*.sorted.bam; do
    [ -e "$bam" ] || continue
    base=$(basename "$bam" .sorted.bam)
    
    echo "Processing reference selection for $base"
    
    # Extract barcode from basename (e.g., "barcode01" from "barcode01_75bccf3a")
    barcode=$(echo "$base" | grep -o 'barcode[0-9]\+')
    
    if [ -z "$barcode" ]; then
        echo "  Warning: Could not extract barcode from $base"
        continue
    fi
    
    # Get idxstats and filter contigs by coverage threshold
    selected_contigs=()
    while IFS=$'\t' read -r contig rlen mapped unmapped; do
        if [ "$contig" != "*" ] && [ "$mapped" -ge "$COVERAGE_THRESHOLD" ]; then
            selected_contigs+=("$contig")
        fi
    done < <(samtools idxstats "$bam")
    
    if [ ${#selected_contigs[@]} -eq 0 ]; then
        echo "  Warning: No contigs with coverage >= $COVERAGE_THRESHOLD for $barcode"
        continue
    fi
    
    # Create per-barcode reference FASTA
    ref_file="$REFS_DIR/${barcode}.fasta"
    > "$ref_file"  # Clear the file
    
    for contig in "${selected_contigs[@]}"; do
        if [ -n "${ref_contigs[$contig]:-}" ]; then
            echo ">${contig}" >> "$ref_file"
            echo "${ref_contigs[$contig]}" >> "$ref_file"
        fi
    done
    
    echo "All personalized references created in $REFS_DIR/"
  done

  mark_step_complete "step_5_personalized_refs"
fi

echo "================================="
echo "STEP 6: Remapping FASTQ to personalized references"
echo "================================="

if check_step "step_6_remapping"; then
    skip_step "Step 6: Remapping"
else
    REMAP_DIR="remap_results"
    CONSOLIDATED_FASTQ_DIR="${OUTPUT_DIR}/fastq_consolidated"
    mkdir -p "$REMAP_DIR"

    if [ ! -d "$CONSOLIDATED_FASTQ_DIR" ]; then
        echo "Warning: Consolidated FASTQ directory not found: $CONSOLIDATED_FASTQ_DIR"
        echo "Skipping STEP 6."
    else
    # For each personalized reference, find and map the corresponding fastq files
    for ref_file in "$REFS_DIR"/*.fasta; do
        [ -e "$ref_file" ] || continue
        
        # Extract barcode from filename (e.g., "barcode01" from "barcode01.fasta")
        barcode=$(basename "$ref_file" .fasta)
        barcode=$(echo "$barcode" | grep -o 'barcode[0-9]\+')
        
        if [ -z "$barcode" ]; then
            echo "  Warning: Could not extract barcode from $(basename "$ref_file")"
            continue
        fi
        
        echo "Processing remapping for $barcode"
        
        # Find consolidated FASTQ for this barcode
        consolidated_fastq="${CONSOLIDATED_FASTQ_DIR}/${barcode}_consolidated.fastq"
        
        if [ ! -f "$consolidated_fastq" ]; then
            echo "  Warning: No consolidated FASTQ found for $barcode at $consolidated_fastq"
            continue
        fi
        
        if [ ! -s "$consolidated_fastq" ]; then
            echo "  Warning: Consolidated FASTQ is empty for $barcode"
            continue
        fi
        
        # Map to personalized reference
        echo "  Mapping $barcode to personalized reference"
        minimap2 -t "$THREADS" -ax map-ont "$ref_file" "$consolidated_fastq" |
            samtools sort -@ "$THREADS" -o "$REMAP_DIR/${barcode}.bam"
        
        # Index the BAM file
        samtools index "$REMAP_DIR/${barcode}.bam"
        
        echo "  Completed remapping for $barcode"
    done
    
    echo "All remapping completed. Results in $REMAP_DIR/"
    fi
fi

echo "================================="
echo "STEP 7: Variant calling and consensus"
echo "================================="

if check_step "step_7_variant_calling"; then
    skip_step "Step 7: Variant calling"
else
    for bam in "$REMAP_DIR"/*.bam; do
        [ -e "$bam" ] || continue
        base=$(basename "$bam" .bam)
        checkpoint_file="${CHECKPOINT_DIR}/step_7_variants_${base}.done"
        
        if [ $FORCE_RESTART -eq 0 ] && [ -f "$checkpoint_file" ]; then
            continue
        fi
    
    # Extract barcode
    barcode=$(echo "$base" | grep -o 'barcode[0-9]\+')
    if [ -z "$barcode" ]; then
        echo "Warning: Could not extract barcode from $base"
        continue
    fi
    
    # Find the personalized reference for this barcode
    ref_file="$REFS_DIR/${barcode}.fasta"
    if [ ! -f "$ref_file" ]; then
        echo "Warning: Reference file not found for $barcode: $ref_file"
        continue
    fi
    
    echo "Calling variants for $barcode"
    
    # Variant calling with filtering
    bcftools mpileup -f "$ref_file" -a AD,DP -d 10000 -Q 10 "$bam" |
        bcftools call -mv --ploidy 1 -Ou |
        bcftools filter -i 'QUAL>20 && DP>20' -Oz -o "$REMAP_DIR/${barcode}.filtered.vcf.gz"
    
    bcftools index -f "$REMAP_DIR/${barcode}.filtered.vcf.gz"
    
    echo "  Completed variant calling for $barcode"
        touch "$checkpoint_file"
    done
    mark_step_complete "step_7_variant_calling"
fi
echo "================================="
for bam in "$REMAP_DIR"/*.bam; do
    [ -e "$bam" ] || continue
    base=$(basename "$bam" .bam)
    checkpoint_file="${CHECKPOINT_DIR}/step_8_zerocov_${base}.done"
    
    if [ $FORCE_RESTART -eq 0 ] && [ -f "$checkpoint_file" ]; then
        continue
    fi
    
    echo "Generating zero-coverage BED for $base"
    bedtools genomecov -ibam "$bam" -bga |
        awk '$4==0' > "$REMAP_DIR/${base}.zero_cov.bed"
    touch "$checkpoint_file"
done

echo "================================="
echo "STEP 9: Generate consensus sequences"
echo "================================="
for filtered_vcf in "$REMAP_DIR"/*.filtered.vcf.gz; do
    [ -e "$filtered_vcf" ] || continue
    base=$(basename "$filtered_vcf" .filtered.vcf.gz)
    checkpoint_file="${CHECKPOINT_DIR}/step_9_consensus_${base}.done"
    
    if [ $FORCE_RESTART -eq 0 ] && [ -f "$checkpoint_file" ]; then
        continue
    fi
    
    # Extract barcode
    barcode=$(echo "$base" | grep -o 'barcode[0-9]\+')
    if [ -z "$barcode" ]; then
        echo "Warning: Could not extract barcode from $base"
        continue
    fi
    
    # Find the personalized reference for this barcode
    ref_file="$REFS_DIR/${barcode}.fasta"
    zero_cov_bed="$REMAP_DIR/${barcode}.zero_cov.bed"
    
    if [ ! -f "$ref_file" ]; then
        echo "Warning: Reference file not found for $barcode"
        continue
    fi
    
    echo "Generating consensus for $barcode"
    
    # Generate consensus with masking for zero-coverage regions
    if [ -f "$zero_cov_bed" ]; then
        bcftools consensus -f "$ref_file" -m "$zero_cov_bed" "$filtered_vcf" > "$REMAP_DIR/${barcode}.consensus.fasta"
    else
        bcftools consensus -f "$ref_file" "$filtered_vcf" > "$REMAP_DIR/${barcode}.consensus.fasta"
    fi
    
    echo "  Completed consensus for $barcode"
    touch "$checkpoint_file"
done

echo "Consensus generation complete. Results in $REMAP_DIR/"

echo "================================="
echo "STEP 10: Create multi-fasta files per segment"
echo "================================="

if check_step "step_10_multifasta"; then
    skip_step "Step 10: Multi-fasta aggregation"
else
    MULTIFASTA_DIR="multi_fasta_segments"
    mkdir -p "$MULTIFASTA_DIR"

# Define segment list (matching the notebook)
declare -a SEGMENT_LIST=(
    'A_MP' 'A_NP' 'A_NS' 'A_PA' 'A_PB1' 'A_PB2'
    'A_HA_H1' 'A_HA_H10' 'A_HA_H11' 'A_HA_H12' 'A_HA_H13' 'A_HA_H14' 'A_HA_H15' 'A_HA_H16'
    'A_HA_H2' 'A_HA_H3' 'A_HA_H4' 'A_HA_H5' 'A_HA_H6' 'A_HA_H7' 'A_HA_H8' 'A_HA_H9'
    'A_NA_N1' 'A_NA_N2' 'A_NA_N3' 'A_NA_N4' 'A_NA_N5' 'A_NA_N6' 'A_NA_N7' 'A_NA_N8' 'A_NA_N9'
    'B_HA' 'B_MP' 'B_NA' 'B_NP' 'B_NS' 'B_PA' 'B_PB1' 'B_PB2'
)

# Initialize multi-fasta files (clear them if they exist)
for segment in "${SEGMENT_LIST[@]}"; do
    > "$MULTIFASTA_DIR/${segment}_cv.fasta"
done

# Process each consensus file
for consensus_file in "$REMAP_DIR"/*.consensus.fasta; do
    [ -e "$consensus_file" ] || continue
    
    # Extract barcode from filename (e.g., "barcode01" from "barcode01.consensus.fasta")
    barcode=$(basename "$consensus_file" .consensus.fasta)
    barcode=$(echo "$barcode" | grep -o 'barcode[0-9]\+')
    
    if [ -z "$barcode" ]; then
        echo "  Warning: Could not extract barcode from $(basename "$consensus_file")"
        continue
    fi
    
    # Parse consensus file and extract sequences
    current_segment=""
    current_seq=""
    
    while IFS= read -r line; do
        # Check if this is a header line
        if [[ "$line" =~ ^">"(.*)$ ]]; then
            # If we have a previous sequence, append it to the appropriate segment file
            if [ -n "$current_segment" ] && [ -n "$current_seq" ]; then
                # Find matching segment file
                for segment in "${SEGMENT_LIST[@]}"; do
                    if [[ "$current_segment" == "$segment" ]]; then
                        echo ">${current_segment}_${barcode}_cv" >> "$MULTIFASTA_DIR/${segment}_cv.fasta"
                        echo "$current_seq" >> "$MULTIFASTA_DIR/${segment}_cv.fasta"
                        break
                    fi
                done
            fi
            
            # Start new sequence
            current_segment="${BASH_REMATCH[1]}"
            current_seq=""
        else
            # Accumulate sequence
            current_seq+="$line"
        fi
    done < "$consensus_file"
    
    # Don't forget the last sequence
    if [ -n "$current_segment" ] && [ -n "$current_seq" ]; then
        for segment in "${SEGMENT_LIST[@]}"; do
            if [[ "$current_segment" == "$segment" ]]; then
                echo ">${current_segment}_${barcode}_cv" >> "$MULTIFASTA_DIR/${segment}_cv.fasta"
                echo "$current_seq" >> "$MULTIFASTA_DIR/${segment}_cv.fasta"
                break
            fi
        done
    fi
done

echo "Multi-fasta aggregation complete. Results in $MULTIFASTA_DIR/"
    mark_step_complete "step_10_multifasta"
fi

echo "================================="
echo "STEP 11: Run Nextclade on segment FASTA files"
echo "================================="

# Switch to nextclade environment
conda deactivate
conda activate nextclade_env
check_command nextclade

OUTPUT_DIR="nextclade_results"
DATASET_CACHE="nextclade_dataset"

mkdir -p "$OUTPUT_DIR" "$DATASET_CACHE"

# Map each HA/NA segment to the Nextclade influenza dataset
declare -A DATASET_BY_SEGMENT=(
    ["A_HA_H1"]="nextstrain/flu/h1n1pdm/ha"
    ["A_HA_H2"]="nextstrain/flu/h2n2/ha"
    ["A_HA_H3"]="nextstrain/flu/h3n2/ha"
    ["A_HA_H5"]="nextstrain/flu/h5n1/ha"
    ["A_HA_H7"]="nextstrain/flu/h7n9/ha"
    ["A_NA_N1"]="nextstrain/flu/h1n1pdm/na"
    ["A_NA_N2"]="nextstrain/flu/h3n2/na"
    ["B_HA"]="nextstrain/flu/vic/ha"
    ["B_NA"]="nextstrain/flu/vic/na"
)

# Process each segment multi-fasta file
for segment_file in "$MULTIFASTA_DIR"/*_cv.fasta; do
    [ -e "$segment_file" ] || continue
    
    fname=$(basename "$segment_file" _cv.fasta)
    
    # Skip if file is empty
    if [ ! -s "$segment_file" ]; then
        echo "Skipping $fname: file is empty"
        continue
    fi
    
    # Check if file has sequence content
    seq_data=$(grep -v '^>' "$segment_file" 2>/dev/null | tr -d '[:space:]' || true)
    if [ -z "$seq_data" ]; then
        echo "Skipping $fname: no sequence content found"
        continue
    fi
    
    # Look up dataset for this segment
    dataset_name=${DATASET_BY_SEGMENT[$fname]:-}
    if [ -z "$dataset_name" ]; then
        echo "Warning: No Nextclade dataset mapping for segment '$fname' - skipping"
        continue
    fi
    
    echo "Processing $fname with dataset $dataset_name"
    
    # Download dataset if not cached
    dataset_dir="$DATASET_CACHE/${dataset_name//\//_}"
    mkdir -p "$dataset_dir"
    
    if [ -z "$(ls -A "$dataset_dir" 2>/dev/null)" ]; then
        echo "  Downloading Nextclade dataset: $dataset_name"
        nextclade dataset get --name "$dataset_name" --output-dir "$dataset_dir"
    else
        echo "  Using cached dataset for $fname"
    fi
    
    # Run Nextclade
    segment_output_dir="$OUTPUT_DIR/${fname}_cv"
    mkdir -p "$segment_output_dir"
    
    echo "  Running Nextclade on $fname"
    nextclade run \
        --input-dataset "$dataset_dir" \
        --output-all "$segment_output_dir" \
        --output-basename "${fname}_cv" \
        --output-selection all \
        --output-tsv "$segment_output_dir/${fname}_cv.tsv" \
        --output-json "$segment_output_dir/${fname}_cv.json" \
        --output-ndjson "$segment_output_dir/${fname}_cv.ndjson" \
        "$segment_file"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Completed Nextclade for $fname"
    else
        echo "  ✗ Error running Nextclade for $fname"
    fi
done

echo "Nextclade processing complete. Results in $OUTPUT_DIR/"

# Deactivate nextclade environment

echo "================================="
echo "STEP 12: Generate analysis plots and tables"
echo "================================="

if [ -f "plot_results.py" ] && command -v python3 >/dev/null 2>&1; then
    echo "Running analysis and visualization script..."
    python3 plot_results.py
    if [ $? -eq 0 ]; then
        echo "✓ Analysis plots generated successfully"
    else
        echo "Warning: Analysis script completed with errors"
    fi
else
    if [ ! -f "plot_results.py" ]; then
        echo "Warning: plot_results.py not found in current directory"
    fi
    if ! command -v python3 >/dev/null 2>&1; then
        echo "Warning: python3 not available"
    fi
    echo "Skipping analysis plots. Run plot_results.py manually if needed."
fi

# Deactivate all conda environments
conda deactivate

echo "================================="
echo "PIPELINE COMPLETE"
echo "================================="
echo ""
echo "Log file: $LOG_FILE"
echo "Results directory: $OUTPUT_DIR"
echo ""
