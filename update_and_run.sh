#!/bin/bash

# Script to rebuild and run the cloupecomplement Docker container

IMAGE="cloupecomplement"

echo "🔄 Stopping any running containers..."
docker stop $(docker ps -q --filter ancestor=$IMAGE) 2>/dev/null

echo "🔨 Rebuilding Docker image with latest changes..."
docker build -t $IMAGE .

if [ $? -eq 0 ]; then
    echo "✅ Build successful! Starting the app..."
    ./cloupecomplement.sh
else
    echo "❌ Build failed. Please check the error messages above."
    exit 1
fi
