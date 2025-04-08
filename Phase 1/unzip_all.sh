#!/bin/bash
# Move to the parent directory where mmCIF is located
cd "$(dirname "$0")/.."

# Optional: Verify that the mmCIF folder exists
if [ ! -d "./mmCIF" ]; then
    echo "Error: mmCIF directory not found in $(pwd)"
    exit 1
fi

# Unzip all .gz files inside the mmCIF folder
find ./mmCIF -name '*.gz' -exec gunzip {} \;