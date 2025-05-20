#!/bin/bash

id="GSM7747185"
cellranger count --id $id \
                --create-bam false \
                --output-dir "./data/cellranger/$id" \
                --transcriptome ~/apps/refdata-gex-GRCm39-2024-A \
                --fastqs "./data/raw/$id" \
                --sample SRR25819158 \
                --localcores 15 \
                --localmem 100

id="GSM7747186"
cellranger multi --id $id \
                --output-dir "./data/cellranger/$id" \
                --transcriptome ~/apps/refdata-gex-GRCm39-2024-A \
                --fastqs "./data/raw/$id" \
                --sample SRR25819156,SRR25819157 \
                --localcores 15 \
                --localmem 100


id="GSM7747187"
cellranger count --id $id \
                --create-bam false \
                --output-dir "./data/cellranger/$id" \
                --transcriptome ~/apps/refdata-gex-GRCm39-2024-A \
                --fastqs "./data/raw/$id" \
                --sample SRR25819158 \
                --localcores 15 \
                --localmem 100
                
id="GSM7747188"
cellranger count --id $id \
                --create-bam false \
                --output-dir "./data/cellranger/$id" \
                --transcriptome ~/apps/refdata-gex-GRCm39-2024-A \
                --fastqs "./data/raw/$id" \
                --sample SRR25819158 \
                --localcores 15 \
                --localmem 100