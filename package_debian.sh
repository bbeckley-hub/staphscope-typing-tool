#!/bin/bash
# package_debian.sh — build Debian package for staphscope-typing-tool

set -e

ROOT_DIR=$(pwd)
echo "Building Debian package in $ROOT_DIR"

# Check for debian/ directory
if [ ! -d "debian" ]; then
    echo "ERROR: debian/ directory not found!"
    exit 1
fi

# Only update changelog if it does NOT exist
if [ ! -f "debian/changelog" ]; then
    if command -v dch &> /dev/null; then
        echo "Creating changelog..."
        dch --create -v 0.1.0-1 "Initial Debian release." --package staphscope-typing-tool
    fi
else
    echo "Using existing changelog."
fi

chmod +x debian/rules

# Clean previous build artifacts
echo "Cleaning previous builds..."
debuild clean || true

# Build the source and binary packages
echo "Building Debian package..."
debuild -us -uc -ui

echo "Debian package build completed!"
echo "Packages should be in the parent directory:"
ls -lh ../staphscope-typing-tool_0.1.0-1_all.deb || true
