#!/bin/bash

# This script reproduces a representative figure.
# Run this script after the project has been built (using the grsi-install.sh script or by following the instructions in README.md).

# Exit on errors.
set -o errexit

echo ""
echo "This program will reproduce the contents of Fig. 8 (a) and (b) in the paper."
echo ""

# Run the program.
cd build
./functional_maps
