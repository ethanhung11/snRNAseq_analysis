#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"/../..

# ================== PROCESS ARGS ==================

# default args
name="preprocessing"
input_dir="default"
input_type="default"

show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i, --input DIR         Files to process (required)"
    echo "  -n, --name FILE         Processing name (required)"
    echo "  -t, --inputtype FILE    Input filetype (required)"
    echo "  -o, --output FILE       Output file (optional)"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --input ./data/cellbender/testInput --output preprocessing.out"
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
        -n|--name)
            name="$2"
            shift 2
            ;;
        -t|--inputtype)
            input_type="$2"
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

if [ -n "$output_file" ]; then
    output=./outs/"$output_file"
    exec > "$output"

    echo "========== Script Output =========="
    echo "Date: $(date)"
    echo "Processing directory: $directory"
    echo "==================================="
fi

input=""
if [ "$input_dir" == "default" ] || [ "$input_type" == "default" ]; then
    input="$input_dir"
elif [ "$input_type" == "h5" ]; then
    for file in `ls ./$input_dir`; do
        input+="$input_dir/$file/output_filtered_seurat.h5,"
    done
elif [ "$input_type" == "folder" ]; then
    for file in `ls ./$input_dir`; do
        input+="$input_dir/$file,"
    done
fi

echo $input

# ================== BEGIN SCRIPT ==================

Rscript --version
time Rscript ./src/R/preprocess.R -i "$input" -n "$name" -t "$input_type"
echo "COMPLETE!"