#!/bin/bash
# database.sh - Fetch AMRFinder binaries and download the latest database
set -e

REPO_URL="https://github.com/bbeckley-hub/amr.git"
BRANCH="master"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Fetching AMR binaries and helper tools from $REPO_URL"
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Clone only the latest commit (shallow)
git clone --depth 1 --branch "$BRANCH" "$REPO_URL" "$TEMP_DIR"

# Update bin/ and stx/ (overwrite with fork's versions)
for dir in bin stx; do
    if [ -d "$TEMP_DIR/$dir" ]; then
        echo "Updating $dir/"
        rm -rf "$SCRIPT_DIR/$dir"
        cp -r "$TEMP_DIR/$dir" "$SCRIPT_DIR/"
    fi
done

# Ensure the data directory exists
mkdir -p "$SCRIPT_DIR/data/amrfinder_db"

# Download the latest AMR database using amrfinder_update
echo "Downloading the latest AMR database (this may take a few minutes)..."
"$SCRIPT_DIR/bin/amrfinder_update" --database "$SCRIPT_DIR/data/amrfinder_db"

echo "AMR setup complete."
