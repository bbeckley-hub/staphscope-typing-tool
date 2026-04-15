#!/bin/bash
# mlst_database.sh - Fetch MLST database and helper directories from GitHub fork
# Run from the directory containing mlst_module.py

set -e  # exit on error

REPO_URL="https://github.com/bbeckley-hub/mlst.git"
BRANCH="master"          # your fork's default branch
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Fetching MLST data from $REPO_URL (branch $BRANCH)"

# Create a temporary directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Shallow clone the repository (only the latest commit)
git clone --depth 1 --branch "$BRANCH" "$REPO_URL" "$TEMP_DIR"

# Replace the four directories (bin, db, perl5, scripts)
for dir in bin db perl5 scripts; do
    if [ -d "$TEMP_DIR/$dir" ]; then
        echo "Updating $dir/"
        rm -rf "$SCRIPT_DIR/$dir"
        cp -r "$TEMP_DIR/$dir" "$SCRIPT_DIR/"
    else
        echo "Warning: $dir not found in the repository"
    fi
done

echo "Database update complete."
