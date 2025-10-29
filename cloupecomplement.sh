#!/bin/bash

# cLoupeComplement launcher for Linux
# Requirements: Docker must be installed and running

IMAGE="nicolasalaun/cloupecomplement:latest"
IMAGE_NAME="nicolasalaun/cloupecomplement"

echo ""
echo "========================================"
echo "  cLoupeComplement - Linux Launcher"
echo "========================================"
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "ERROR: Docker is not running!"
    echo "Please start Docker and try again."
    echo ""
    read -p "Press Enter to exit..."
    exit 1
fi

# Check if image exists
if ! docker image inspect $IMAGE > /dev/null 2>&1; then
    echo "Docker image not found locally."
    echo "Downloading from Docker Hub..."
    echo "This may take several minutes (image size: ~2GB compressed)..."
    echo ""

    docker pull $IMAGE

    if [ $? -eq 0 ]; then
        echo ""
        echo "Image downloaded successfully!"
        echo ""
    else
        echo ""
        echo "ERROR: Failed to download image from Docker Hub"
        echo "Please check your internet connection and try again."
        echo ""
        read -p "Press Enter to exit..."
        exit 1
    fi
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

# Open browser (suppress GTK warnings)
xdg-open http://localhost:3838 2>/dev/null

# Show info
echo ""
echo "========================================"
echo "  App Status"
echo "========================================"
echo ""
echo ">> cLoupeComplement is running at: http://localhost:3838"
echo ">> Upload your H5 and CSV files through the web interface"
echo ">> Downloads will save to your browser's default download folder"
echo ""
echo "IMPORTANT: Keep this terminal open!"
echo "Close this terminal to stop the app."
echo ""
echo "========================================"
echo ""

# Wait for container to exit
wait $CONTAINER_PID
