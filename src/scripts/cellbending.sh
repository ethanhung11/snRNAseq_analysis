#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory_path> [output_file]"
    exit 1
fi

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"/../..

# ================== PROCESS ARGS ==================

# default args
cores=10

show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i, --input DIR         Directory > experiments > samples to process (required)"
    echo "  -o, --output FILE       Output file (optional)"
    echo "  --cores CORES           Cores (optional), default is 10"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --input experiment --output output.out"
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            input_dir="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

if [ -v input_dir ]; then
  "No input entered! Please use the -i / --input flag."
fi

if [ -n "$output_file" ]; then
    output=./outs/"$output_file"
    exec > "$output"

    echo "========== Script Output =========="
    echo "Date: $(date)"
    echo "Processing directory: $directory"
    echo "==================================="
fi

# ================== BEGIN SCRIPT ==================

echo
echo "Rye=="`rye version`

echo "Searching $input_dir for experiments..."
for experimentdir in "$input_dir"/*; do
    echo "Searching $experimentdir for samples..."

    if ! [ -d "$experimentdir" ]; then
        echo "Experiment $experimentdir not found!"
        exit 1
    fi

    echo "Samples found in $experimentdir: $(ls $experimentdir)"

    for sampledir in "$experimentdir"/*; do
        
        if [ -f "$sampledir/outs/filtered_feature_bc_matrix.h5" ]; then
            sampleFile="$sampledir/outs/filtered_feature_bc_matrix.h5"
        elif [ -d "$sampledir/outs/filtered_feature_bc_matrix" ]; then
            sampleFile="$sampledir/outs/filtered_feature_bc_matrix"
        else
            echo "Sample file not found within directory!"
            echo "'filtered_feature_bc_matrix' or 'filtered_feature_bc_matrix.h5' must exist at 'dir/exp/samp/outs/'!"
            exit 1
        fi
        
        sampleOutput="./data/cellbender/$(basename $sampledir)"
        mkdir -p "$sampleOutput"

        time rye run cellbender remove-background \
        --input "$sampleFile" \
        --output "$sampleOutput/output.h5" \
        --cpu-threads $cores \
        --checkpoint-mins 10 \
        --low-count-threshold 10

        rye run ptrepack --complevel 5 "./$sampleOutput/output.h5:/matrix" "./$sampleOutput/output_seurat.h5:/matrix"
        rye run ptrepack --complevel 5 "./$sampleOutput/output_filtered.h5:/matrix" "./$sampleOutput/output_filtered_seurat.h5:/matrix"
    done
done
echo "COMPLETE!"
echo