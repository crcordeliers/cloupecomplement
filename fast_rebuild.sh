#!/bin/bash

# Fast Docker rebuild script with BuildKit and parallel builds

IMAGE="cloupecomplement"

echo "🚀 Building with BuildKit (faster, parallel builds)..."

# Enable BuildKit for better caching and parallelization
export DOCKER_BUILDKIT=1

# Build with progress
docker build \
  --progress=plain \
  --build-arg BUILDKIT_INLINE_CACHE=1 \
  -t $IMAGE \
  .

if [ $? -eq 0 ]; then
    echo "✅ Build successful!"
    echo "🔄 Stopping old containers..."
    docker stop $(docker ps -q --filter ancestor=$IMAGE) 2>/dev/null

    echo "▶️  Starting app..."
    ./cloupecomplement.sh
else
    echo "❌ Build failed"
    exit 1
fi
