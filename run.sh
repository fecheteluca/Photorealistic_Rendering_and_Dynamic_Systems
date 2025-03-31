#!/bin/bash
# Exit immediately if any command fails.
set -e

# Define variables
BUILD_DIR="build"
APP_DIR="app"
EXECUTABLE_NAME="main"  

# Change to the build directory
cd "${BUILD_DIR}"

# Configure the project 
echo "Configuring project with CMake..."
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build the project
echo "Building project..."
cmake --build .

cd "${APP_DIR}"

# Attempt to run the executable
echo "Running the executable..."
if [ -f "${EXECUTABLE_NAME}" ]; then
    ./"${EXECUTABLE_NAME}"
    echo "Executable run success"
else
    echo "Error: Executable '${EXECUTABLE_NAME}' not found in $(pwd)"
fi
