# Source Code Documentation (`src/`)

This directory contains the core logic and mathematical framework for the coupled electromagnetic dipole formalism.

## File Descriptions

*   **`2DSSHlattice.jl`**: Calculates the lattice site coordinates for the 2DSSH array, depending on size and beta.
*   ** `GreenMatrix.jl`**: Calculates the Green's matrix for the forementioned lattice.
*   ** `eigensolver.jl`**: Solves the equations to find eigenmodes and eigenfrequencies.
*   **`ClassifyModes.jl`**: Functions to classify eigenmodes as bulk, edge, and corner type.
*   **`FindGammaModes.jl`**: Find the eigenmodes which correspond to the Gamma point for different symmetries/localization.
*   **`Sortinq.jl`**: Sorts eigenmodes in momentum (q) space relying on Fourier transforms of the eigenmodes.
*   **`DispersionBands.jl`**: Calculates dispersion and Q-factor bands using eigensolver.jl and Sortinq.jl.
*   **`RadiationPatterns.jl`**: Calculates radiation patterns for any eigenmode, or concretely Gamma modes.
*   **`Extinction.jl`**: Calculates extinction cross sections of the finite array.

*   **`parameters.jl`**: Contains the general parameters, physical constants and used formulas for figures.
