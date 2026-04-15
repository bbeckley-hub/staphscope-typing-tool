#!/usr/bin/env bash
################################################################################
# sccmec_database.sh - Download/Update SCCmecFinder database
# Author: Brown Beckley <brownbeckley94@gmail.com>
# Description: Clones the latest SCCmec database from Bitbucket into the 
#              'database' folder of the sccmec_module.
# Usage: ./sccmec_database.sh
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATABASE_DIR="${SCRIPT_DIR}/database"
TEMP_CLONE_DIR="/tmp/sccmec_db_$$"

# Repository URL (master branch)
REPO_URL="https://bitbucket.org/genomicepidemiology/sccmecfinder_db.git"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  SCCmec Database Update Tool${NC}"
echo -e "${GREEN}========================================${NC}"

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo -e "${RED}Error: git is not installed. Please install git first.${NC}"
    exit 1
fi

echo -e "${YELLOW}Step 1: Removing old database folder (if exists)...${NC}"
if [ -d "$DATABASE_DIR" ]; then
    rm -rf "$DATABASE_DIR"
    echo -e "${GREEN}  ✓ Removed old database folder.${NC}"
fi

echo -e "${YELLOW}Step 2: Creating temporary clone directory...${NC}"
mkdir -p "$TEMP_CLONE_DIR"
echo -e "${GREEN}  ✓ Temporary directory: $TEMP_CLONE_DIR${NC}"

echo -e "${YELLOW}Step 3: Cloning SCCmec database repository (shallow clone)...${NC}"
git clone --depth 1 "$REPO_URL" "$TEMP_CLONE_DIR"
if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Failed to clone repository. Check your internet connection.${NC}"
    rm -rf "$TEMP_CLONE_DIR"
    exit 1
fi
echo -e "${GREEN}  ✓ Clone successful.${NC}"

echo -e "${YELLOW}Step 4: Copying database files to $DATABASE_DIR...${NC}"
# The repository contains the database files directly in the root (not in a subfolder)
# We copy the three .fasta files and the template_db directory.
mkdir -p "$DATABASE_DIR"
cp "$TEMP_CLONE_DIR"/*.fasta "$DATABASE_DIR/" 2>/dev/null || echo "  (No .fasta files found?)"
if [ -d "$TEMP_CLONE_DIR/template_db" ]; then
    cp -r "$TEMP_CLONE_DIR/template_db" "$DATABASE_DIR/"
    echo -e "${GREEN}  ✓ Copied template_db directory.${NC}"
else
    echo -e "${YELLOW}  ⚠️ template_db directory not found in clone.${NC}"
fi

# Also copy any other needed files (like __init__.py if present)
if [ -f "$TEMP_CLONE_DIR/__init__.py" ]; then
    cp "$TEMP_CLONE_DIR/__init__.py" "$DATABASE_DIR/"
fi

# List the contents to verify
echo -e "${YELLOW}Step 5: Verifying database contents...${NC}"
ls -la "$DATABASE_DIR"

echo -e "${YELLOW}Step 6: Cleaning up temporary files...${NC}"
rm -rf "$TEMP_CLONE_DIR"
echo -e "${GREEN}  ✓ Temporary files removed.${NC}"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  SCCmec database update completed!${NC}"
echo -e "${GREEN}  Database location: $DATABASE_DIR${NC}"
echo -e "${GREEN}========================================${NC}"
