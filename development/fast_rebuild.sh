#!/bin/bash

# Fast Docker rebuild script with BuildKit and parallel builds
# This script builds the app and saves it as a tar.gz for distribution

IMAGE="cloupecomplement"
IMAGE_FILE="docker/cloupecomplement-docker-image.tar.gz"

echo "🚀 Building with BuildKit (faster, parallel builds)..."

# Enable BuildKit for better caching and parallelization
export DOCKER_BUILDKIT=1

# Build with progress
docker build \
  --progress=plain \
  --build-arg BUILDKIT_INLINE_CACHE=1 \
  -f docker/dockerfile \
  -t $IMAGE \
  .

if [ $? -eq 0 ]; then
    echo "✅ Build successful!"

    echo "💾 Saving Docker image to $IMAGE_FILE..."
    docker save $IMAGE | gzip > $IMAGE_FILE

    if [ $? -eq 0 ]; then
        echo "✅ Image saved successfully!"
        SIZE=$(du -h $IMAGE_FILE | cut -f1)
        echo "📦 File size: $SIZE"
    else
        echo "❌ Failed to save image"
        exit 1
    fi

    echo "🔄 Stopping old containers..."
    docker stop $(docker ps -q --filter ancestor=$IMAGE) 2>/dev/null

    echo "▶️  Starting app..."
    ./cloupecomplement.sh
else
    echo "❌ Build failed"
    exit 1
fi
