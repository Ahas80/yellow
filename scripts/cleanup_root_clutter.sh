#!/usr/bin/env bash
set -euo pipefail

# Define the target directory for clutter
TARGET_DIR="results/wgs/core/legacy_clutter"

# Create the directory if it doesn't exist
echo "Creating directory: $TARGET_DIR"
mkdir -p "$TARGET_DIR"

# Move multi_*.fa.* files
echo "Moving multi_*.fa.* files to $TARGET_DIR..."
count=$(ls multi_*.fa.* 2>/dev/null | wc -l | xargs)

if [ "$count" -gt 0 ]; then
    mv multi_*.fa.* "$TARGET_DIR/"
    echo "Successfully moved $count files."
else
    echo "No multi_*.fa.* files found in the root directory."
fi

# Optional: Check for other potential clutter mentioned in findings (e.g. miniforge installer)
# For now, we only strictly move the multi_ files as per the plan.

echo "Cleanup complete."
ls -ld "$TARGET_DIR"
