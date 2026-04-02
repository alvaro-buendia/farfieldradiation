# Reproducing Figures for the article "Far-field radiation of bulk, edge and corner eigenmodes from a finite 2d su-schrieffer-heeger plasmonic lattice" by the authors, Á. Buendía, J.L. Pura, V. Giannini and J.A. Sánchez-Gil, published in PRB (doi: 10.1103/54vz-h71c).

This guide explains how to set up the Julia environment and run the scripts to reproduce the figures stored in the `results/` directory.

## 📁 Project Structure

*   **`src/`**: Contains the core source code and physics (modules, functions, and calculations).
*   **`scripts/`**: Julia scripts (`.jl`) specifically designed to generate the figures using the source code.
*   **`results/`**: The output directory where generated figures (PNG, PDF, or SVG) are saved.
*   **`Project.toml` / `Manifest.toml`**: Files defining the exact environment and dependencies needed.

## 🚀 How to Run

### 1. Prerequisites
Ensure you have [Julia](https://julialang.org) installed. Open your terminal in the project root folder.

### 2. Instantiate the Environment
Before running the scripts, you must install the required packages defined in the project files:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'


## To generate the figures

You can run figures independently from terminal

julia --project=. scripts/Figure1.jl

You can also generate all the figures, by typing 

julia --project=. scripts/run_all.jl
