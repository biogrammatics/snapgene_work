#!/bin/bash
#
# Update SnapGene FASTA database for BLAST
# Runs weekly via cron, deletes old file only on success
#

SNAPGENE_DIR="/Users/studio/Library/CloudStorage/Dropbox/SnapGene"
SCRIPT_DIR="/Users/studio/Claude/snapgene-work"
DATE=$(date +%Y-%m-%d)
NEW_FILE="${SNAPGENE_DIR}/${DATE}-snapgene-collection.fa"
LOG_FILE="${SCRIPT_DIR}/snapgene_fasta.log"

# Log start
echo "$(date '+%Y-%m-%d %H:%M:%S') - Starting FASTA update" >> "$LOG_FILE"

# Run the extraction
if python3 "${SCRIPT_DIR}/snapgene_to_fasta.py" --quiet --output "$NEW_FILE"; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Successfully created $NEW_FILE" >> "$LOG_FILE"

    # Delete any older .fa files (keep only today's)
    find "$SNAPGENE_DIR" -maxdepth 1 -name "*-snapgene-collection.fa" ! -name "${DATE}-snapgene-collection.fa" -type f -delete

    # Log cleanup
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Cleaned up old FASTA files" >> "$LOG_FILE"
else
    echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: FASTA extraction failed" >> "$LOG_FILE"
    exit 1
fi
