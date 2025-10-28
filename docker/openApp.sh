#!/bin/bash

# Start RStudio inside a Docker container
docker run --rm \
  -p 8787:8787 \
  -e DISABLE_AUTH=true \
  -v "$HOME:/home/rstudio" \
  rocker/rstudio &

CONTAINER_PID=$!

# Wait a bit
sleep 3

# Open browser
xdg-open http://localhost:8787

# Show info
echo ""
echo "🎯 RStudio is running at: http://localhost:8787"
echo "📂 Your Shiny app is inside /home/rstudio/shiny_app."
echo "🔓 Login is disabled — you're automatically logged in."
echo ""
echo "📌 Inside RStudio, run:"
echo "    install.packages('renv')"
echo "    renv::init()"
echo "    renv::snapshot()"
echo ""
echo "🛑 Close this terminal to stop RStudio."
echo ""
