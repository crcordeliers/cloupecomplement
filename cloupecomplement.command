#!/bin/bash

# cLoupeComplement launcher for macOS
# Requirements: Docker Desktop must be installed and running

IMAGE="cloupecomplement"

echo ""
echo "========================================"
echo "  cLoupeComplement - macOS Launcher"
echo "========================================"
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "ERROR: Docker is not running!"
    echo "Please start Docker Desktop and try again."
    echo ""
    read -p "Press Enter to exit..."
    exit 1
fi

echo "Stopping any existing containers..."
docker ps -q --filter ancestor=$IMAGE | xargs -r docker stop > /dev/null 2>&1

echo "Starting cLoupeComplement..."
echo ""

# Start container in background
docker run --rm \
  -p 3838:3838 \
  "$IMAGE" &

CONTAINER_PID=$!

# Wait for app to start
sleep 5

# Open browser (macOS specific)
open http://localhost:3838

echo ""
echo "========================================"
echo "  App Status"
echo "========================================"
echo ""
echo ">> cLoupeComplement is running at: http://localhost:3838"
echo ">> Upload your H5 and CSV files through the web interface"
echo ">> Downloads will save to your browser's default download folder"
echo ""
echo "IMPORTANT: Keep this window open!"
echo "Close this window to stop the app."
echo ""
echo "========================================"
echo ""

# Wait for container to exit
wait $CONTAINER_PID
