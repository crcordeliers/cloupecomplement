#!/bin/bash

IMAGE="cloupecomplement"

# Start Shiny app with data directories mounted
# Mount specific directories to avoid permission issues
docker run --rm \
  -p 3838:3838 \
  -e DISABLE_AUTH=false \
  -v "/home/$USER/Documents:/data:rw" \
  -v "/home/$USER/Downloads:/downloads:rw" \
  "$IMAGE" &

CONTAINER_PID=$!

# Wait a bit
sleep 3

# Open browser (suppress GTK warnings)
xdg-open http://localhost:3838 2>/dev/null

# Show info
echo ""
echo "🎯 RStudio is running at: http://localhost:3838"
echo "🏠 Your home directory is available inside RStudio."
echo "🔓 Login is disabled — you're automatically logged in."
echo ""
echo "🛑 Close this terminal to stop RStudio."
echo ""

# Wait for container to exit
wait $CONTAINER_PID
