#!/usr/bin/env bash
set -euo pipefail

# Install Python dependencies into the active environment
python3 -m pip install --upgrade biopython PyYAML

# Download and process the NCBI taxonomy dump
TAXDUMP_DIR="taxdump"
mkdir -p "$TAXDUMP_DIR"

wget -q https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz \
  -O "$TAXDUMP_DIR/new_taxdump.tar.gz"

tar -xzf "$TAXDUMP_DIR/new_taxdump.tar.gz" -C "$TAXDUMP_DIR" fullnamelineage.dmp
rm "$TAXDUMP_DIR/new_taxdump.tar.gz"

sed -e 's/\s\+|\s\+//g' "$TAXDUMP_DIR/fullnamelineage.dmp" \
  | sed 's/|//2g' \
  | awk 'BEGIN{FS="|"}{print $1,$3,$4,$2}' \
  | sed 's/\s/\t/g' \
  | sort -k1,1 > "$TAXDUMP_DIR/fullnamelineage_taxid_sorted.dmp"

rm "$TAXDUMP_DIR/fullnamelineage.dmp"
