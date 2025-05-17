#! /usr/bin/env sh

rye version

time nice -n 10 rye run cellbender remove-background \
--input ./data/CellRanger/GSM7747185_Chow-eWAT \
--output ./data/CellBender/GSM7747185/output.h5 \
--checkpoint-mins 10 \
--cpu-threads 10

time nice -n 10 rye run cellbender remove-background \
--input ./data/CellRanger/GSM7747186_Chow-iWAT \
--output ./data/CellBender/GSM7747186/output.h5 \
--checkpoint-mins 10 \
--cpu-threads 10

time nice -n 10 run cellbender remove-background \
--input ./data/CellRanger/GSM7747187_Chow-eWAT \
--output ./data/CellBender/GSM7747187/output.h5 \
--checkpoint-mins 10 \
--cpu-threads 10

time nice -n 10 rye run cellbender remove-background \
--input ./data/CellRanger/GSM7747188_Chow-iWAT \
--output ./data/CellBender/GSM7747188/output.h5 \
--checkpoint-mins 10 \
--cpu-threads 10

echo "COMPLETE!"