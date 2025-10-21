# cLoupeComplement - Launch Instructions

This document explains how to launch cLoupeComplement on different operating systems.

## Prerequisites

**All platforms require Docker to be installed and running:**

- **Windows**: [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **macOS**: [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/)
- **Linux**: [Docker Engine](https://docs.docker.com/engine/install/)

## Building the Docker Image (First Time Only)

Before launching the app for the first time, you need to build the Docker image:

```bash
docker build -t cloupecomplement .
```

This will take 15-25 minutes the first time. You only need to do this once (or when the code is updated).

## Launching the App

### Windows

1. **Double-click** `cloupecomplement.bat`
2. A command window will open and the app will start
3. Your browser will automatically open to http://localhost:3838
4. **Keep the command window open** while using the app
5. Close the command window to stop the app

### macOS

1. **Right-click** `cloupecomplement.command` and select "Open"
   - The first time, macOS may ask for permission to run the file
   - Go to System Preferences > Security & Privacy and click "Open Anyway"
2. Terminal will open and the app will start
3. Your browser will automatically open to http://localhost:3838
4. **Keep the Terminal window open** while using the app
5. Close the Terminal window to stop the app

**Alternative**: You can also run from Terminal:
```bash
./cloupecomplement.command
```

### Linux

1. Open a terminal in the project directory
2. Run:
   ```bash
   ./cloupecomplement.sh
   ```
3. Your browser will automatically open to http://localhost:3838
4. **Keep the terminal window open** while using the app
5. Press Ctrl+C or close the terminal to stop the app

## Using the App

1. Once the app is running, navigate to http://localhost:3838
2. **Upload your data files:**
   - **H5 file**: The `filtered_feature_bc_matrix.h5` file from your SpaceRanger/CellRanger output
     - For Visium HD: Find it at `binned_outputs/square_008um/filtered_feature_bc_matrix.h5`
     - For standard Visium: Find it at `filtered_feature_bc_matrix.h5`
   - **Cluster CSV**: Your cluster annotation file (max 100MB)
3. **Downloads** will save to your browser's default download folder

## Troubleshooting

### Docker is not running
- Make sure Docker Desktop (Windows/Mac) or Docker Engine (Linux) is running
- On Windows/Mac, you should see the Docker icon in your system tray

### Port 3838 is already in use
- Another application is using port 3838
- Stop that application or change the port in the launcher script

### App won't start
- Check that the Docker image was built successfully
- Try rebuilding: `docker build -t cloupecomplement .`

### Cannot upload files
- Make sure your H5 file is under 100MB
- Check that you're selecting the correct file type (.h5 or .hdf5)

## Fast Rebuild (For Developers)

If you make changes to the code (ui.R, server.R, functions.R), you can quickly rebuild:

```bash
./fast_rebuild.sh
```

This only takes 5-10 seconds since it reuses cached Docker layers.
