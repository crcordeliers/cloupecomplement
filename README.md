
# Shiny App for Spatial Transcriptomic Data Exploration

## Overview

This Shiny app is designed to complement 10x Genomics Loupe Browser for spatial transcriptomic data analysis. It provides an intuitive interface for researchers to perform additional data preprocessing, visualization, and downstream analyses, such as differential expression and pathway analysis, after manual clustering in Loupe Browser.

### Key Features:

- **Data Preprocessing:** Load the `outs` folder from your 10x Genomics pipeline and a CSV file containing cluster information. Filter genes based on a minimum expression threshold.
- **Visualization:** Generate violin plots, beeswarm plots, heatmaps, and dot plots to visualize gene expression across clusters. Pairwise comparisons between clusters can be added interactively for violin and beeswarm plots.
- **Differential Expression Analysis:** Perform differential expression analysis between selected clusters.
- **Pathway Analysis:** Explore enriched pathways using results from the differential expression analysis, with interactive visualization of pathway enrichment.

## Installation and Usage

### Prerequisites:

**Docker** must be installed on your system:
- **Windows**: [Docker Desktop](https://docs.docker.com/desktop/install/windows-install/)
- **macOS**: [Docker Desktop](https://docs.docker.com/desktop/install/mac-install/)
- **Linux**: [Docker Engine](https://docs.docker.com/engine/install/)

All R package dependencies are automatically included in the Docker image (no manual installation needed).

### Quick Start (Using Docker - Recommended):

1. **Clone this repository:**
   ```bash
   git clone https://github.com/crcordeliers/cloupecomplement.git
   cd cloupecomplement
   ```

2. **Launch the app (one of the following methods):**

   **Method A: Using Pre-built Image (Fastest - if provided)**
   - If `cloupecomplement-docker-image.tar.gz` is included:
     - **Windows**: Double-click `cloupecomplement.bat`
     - **macOS**: Double-click `cloupecomplement.command`
     - **Linux**: Run `./cloupecomplement.sh`
   - The launcher will automatically load the pre-built image (first time takes 2-3 minutes)

   **Method B: Build from Source**
   - If you don't have the pre-built image file:
     - **Windows**: Double-click `build.bat` (takes 15-25 minutes)
     - **macOS**: Double-click `build.command` (takes 15-25 minutes)
     - **Linux**: Run `./build.sh` (takes 15-25 minutes)
   - Then use the launcher scripts above

3. The app will automatically open in your browser at http://localhost:3838

### Input Files:

- **H5 File:** Upload the `filtered_feature_bc_matrix.h5` file from your SpaceRanger/CellRanger output
  - For **Visium HD**: Located at `binned_outputs/square_008um/filtered_feature_bc_matrix.h5`
  - For **standard Visium**: Located at `filtered_feature_bc_matrix.h5`
  - Maximum file size: 100MB

- **Cluster CSV:** The CSV file exported from Loupe Browser containing cluster assignments
  - First column should contain cluster identifiers
  - Row names should match cell barcodes
  - Maximum file size: 100MB

### Alternative: Running with RStudio (Not Recommended):

You can also run the app directly in RStudio, but this requires manual installation of all dependencies. The Docker method is strongly recommended for ease of use.

1. Open `app.R` in RStudio
2. Click "Run App"

### App Workflow:

1. **Data Loading:**
   - Upload the H5 file (`filtered_feature_bc_matrix.h5`) containing the gene expression matrix
   - Upload the CSV file with cluster information exported from Loupe Browser
   - Configure preprocessing options:
     - Select species (Human or Mouse)
     - Choose normalization method (LogNormalize or SCTransform)
     - Set minimum gene expression threshold (% of cells)
     - Set minimum genes per spot threshold

2. **Visualization:**
   - **Violin & Beeswarm Plots Tab:** Select genes of interest and visualize their expression across clusters with the option to perform pairwise comparisons. Download the plots in PDF format.
   - **Heatmap & Dotplot Tab:** Visualize the expression of selected genes as heatmaps and dotplots, with clusters defined in the metadata.

3. **Differential Expression:**
   - Perform differential expression analysis between selected clusters.
   - View and download the results in table format.

4. **Pathway Analysis:**
   - Choose a pathway analysis method and run enrichment analysis based on the differential expression results. The pathway analysis will be done on the cluster selected in the Differential Expression tab.
   - View and download the results in table and pdf format.

## Troubleshooting

For detailed troubleshooting, see [docs/LAUNCHER_README.md](docs/LAUNCHER_README.md).

**Common issues:**
- **"Docker image not found"**: Run the build script first (`build.bat`, `build.command`, or `build.sh`)
- **"Docker is not running"**: Start Docker Desktop or Docker Engine
- **Cannot upload files**: Check file size (max 100MB) and format (.h5 for matrix, .csv for clusters)
- **Port 3838 already in use**: Stop other applications using that port

## System Requirements

- Docker Desktop (Windows/Mac) or Docker Engine (Linux)
- 8GB RAM minimum (16GB recommended for large datasets)
- ~7GB disk space for Docker image
- Modern web browser (Chrome, Firefox, Safari, Edge)

## Roadmap

### Planned Features:
- **Customizable Pathway Visualization:** Add an option to control the number of categories shown in the pathway analysis plots, with automatic height scaling to accommodate the chosen number.
- **Method Documentation:** Add a way to give the user the mat & med section of the analysis based on their selection in the app's choices.

## Contributions

Feel free to contribute to this project by submitting pull requests or suggesting new features and improvements.
