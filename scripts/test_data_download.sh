#!/bin/bash
set -e  # Exit if any command fails

# Define your bucket path
BUCKET_NAME="pdx_somatic_testing_data"
BUCKET_SUBDIR="data"

# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate one level up to get the project root
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Set destination directory relative to project root
DEST_DIR="$PROJECT_ROOT/tests/data"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

echo "Downloading from gs://$BUCKET_NAME/$BUCKET_SUBDIR to $DEST_DIR..."

# Use gsutil to recursively copy everything from the bucket subdirectory
gsutil -m cp -r "gs://$BUCKET_NAME/$BUCKET_SUBDIR/*" "$DEST_DIR/"

echo "Download complete!"
