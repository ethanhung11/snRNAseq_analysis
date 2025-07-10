#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory_path> [output_file]"
    exit 1
fi

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"/../..

# ================== PROCESS ARGS ==================

# arg defaults
transcriptome="./references/refdata-gex-GRCm39-2024-A"
cores=10
memusage=5d0

show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i, --input DIR         Directory to process (required)"
    echo "  -o, --output FILE       Output file (required)"
    echo "  --transcriptome FILE    Transcriptome (optional), default is ./data/refdata-gex-GRCm39-2024-A"
    echo "  --cores CORES           Cores (optional), default is 10"
    echo "  --mem MBS               Memory Usage (optional), default is 50"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --input experiment_dir --output output.out"
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            input_files="$2"
            shift 2
            ;;
        -o|--output)
            output_file="$2"
            shift 2
            ;;
        --transcriptome)
            transcriptome="$2"
            shift 2
            ;;
        --cores)
            cores="$2"
            shift 2
            ;;
        --mem)
            memusage="$2"
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

if [ -n "$readout" ]; then
    output_file=./outs/"$readout"
    mkdir -p "$output_file"
    exec > "$output_file"

    echo "========== Script Output =========="
    echo "Date: $(date)"
    echo "Processing directory: $directory"
    echo "==================================="
fi

# ================== BEGIN SCRIPT ==================

for experiment in $(echo $input_files | awk -F, '{for (i=1; i<=NF; i++) print $i}' ); do

    inputdir="./data/raw/$experiment"

    if ! [ -d "$inputdir" ]; then
        echo "Experiment $experiment does not exist in './data/raw/'."
        exit 1
    fi
 
    outputdir="./data/cellranger/$experiment"
    mkdir -p "$outputdir"
    cellranger count --id "$experiment" \
            --create-bam false \
            --output-dir "$outputdir" \
            --transcriptome $transcriptome \
            --fastqs "$inputdir" \
            # --sample "$sample" \
            --localcores $cores \
            --localmem $memusage
done

echo "COMPLETE!"