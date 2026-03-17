# FSSLibrary

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://cmoe-tudelft.github.io/FSSLibrary/)

A comprehensive Python library designed for the **Floating and Submerged Structures** course at TU Delft. This library provides essential tools for marine and ocean engineering calculations, including linear wave theory, structural beam analysis, mooring systems, and wave spectrum processing.

## Features

### 🌊 Wave Analysis
- **Linear Wave Theory**: Dispersion relations, wave length calculations, and 2D linear wave analysis
- **Wave Spectrum Processing**: Comprehensive spectral analysis tools for irregular waves
- **FFT-based Signal Processing**: Single-sided spectrum computation and frequency domain analysis

### 🏗️ Structural Analysis
- **Beam Matrices**: Mass, stiffness, and damping matrices for both 2D and 3D beam elements
- **Finite Element Support**: Essential matrices for structural dynamics of marine structures

### ⚓ Mooring Systems
- **Catenary Analysis**: Mooring line calculations with extensible cable dynamics
- **Static Analysis**: Position and tension calculations for mooring systems

### 🧮 Utilities
- **Physical Constants**: Standard marine engineering constants (gravity, water density, etc.)
- **Optimization Tools**: Integrated scipy-based solvers for engineering problems

## Installation

### Prerequisites
- Python 3.7 or higher
- Required dependencies: numpy, scipy, ipykernel

### Install from source
```bash
git clone https://github.com/CMOE-TUDelft/FSSLibrary.git
cd FSSLibrary
python -m pip install -e .
```

### Install dependencies
```bash
python -m pip install -r requirements.txt
```

## Quick Start

```python
import FSSLibrary as fss
import numpy as np

# Linear wave analysis
wave = fss.LinearWave.LinearWave2D(d=50.0, T=8.0, H=2.0)
print(f"Wave length: {wave.L:.2f} m")

# Wave spectrum processing
fs = 10.0  # Sampling frequency
signal = np.random.randn(1024)  # Example signal
f, amp, psd = fss.FFTBasic.get_single_sided_spectrum(signal, fs)

# Access physical constants
print(f"Gravity: {fss.Constants.g} m/s²")
print(f"Water density: {fss.Constants.rhoW} kg/m³")
```

## Module Overview

| Module | Description |
|--------|-------------|
| `LinearWave` | Linear wave theory, dispersion relations, 2D wave calculations |
| `WaveSpectrum` | Wave spectral analysis and irregular wave generation |
| `BeamMatrices` | 2D/3D beam finite element matrices for structural analysis |
| `MoorLib` | Mooring line catenary analysis and calculations |
| `FFTBasic` | FFT-based signal processing and spectrum analysis |
| `Constants` | Physical constants for marine engineering |

## Documentation

### Generate Documentation with pdoc

To generate and serve the documentation locally using pdoc:

```bash
# Install pdoc if not already installed
python -m pip install pdoc

# Generate static HTML documentation
pdoc --output-dir docs/ src/FSSLibrary

# Or Generate documentation and serve on localhost
pdoc --http localhost:8080 src/FSSLibrary
```

The documentation will be available at `http://localhost:8080` when using the second command.

### Pre-built Documentation

Pre-built HTML documentation is available in the `docs/` directory. Open `docs/index.html` in your browser to view the complete API documentation.

## Examples

### Wave Length Calculation
```python
from FSSLibrary.LinearWave import getWaveLen
from FSSLibrary.Constants import g

# Calculate wave length for given conditions
depth = 50.0  # m
period = 8.0  # s
wave_length = getWaveLen(g, depth, period)
print(f"Wave length: {wave_length:.2f} m")
```

### Beam Matrix Generation
```python
from FSSLibrary.BeamMatrices import Beam2DMatrices

# Define beam properties
mass_per_length = 100.0  # kg/m
axial_stiffness = 1e9    # N
bending_stiffness = 1e6  # N⋅m²
node_coords = ([0, 0], [10, 0])  # Start and end coordinates

# Generate matrices
M, K, Q = Beam2DMatrices(mass_per_length, axial_stiffness, 
                         bending_stiffness, node_coords)
```

### Mooring Analysis
```python
from FSSLibrary.MoorLib import catenary_xz

# Mooring line parameters
weight_per_meter = 50.0   # N/m
elastic_modulus = 1e8     # Pa
vertical_tension = 1000.0 # N
horizontal_tension = 800.0 # N
suspended_length = 100.0   # m

# Calculate catenary position
x, z = catenary_xz(weight_per_meter, elastic_modulus, 
                   vertical_tension, horizontal_tension, suspended_length)
```

## Course Context

This library is specifically designed for the **Floating and Submerged Structures** course and provides:

- Educational implementations of fundamental marine engineering concepts
- Ready-to-use functions for course assignments and projects
- Clear, documented code suitable for learning and teaching
- Integration with Jupyter notebooks for interactive analysis

## Development

### Project Structure
```
FSSLibrary/
├── src/FSSLibrary/     # Main package source
├── docs/               # Pre-built documentation
├── requirements.txt    # Dependencies
├── pyproject.toml     # Package configuration
└── README.md          # This file
```

### Testing
```bash
# Run basic import test
python -c "import FSSLibrary; print('Import successful')"
```

## Authors

- **Shagun Agarwal** - shagun.1994@gmail.com
- **Oriol Colomés** - j.o.colomesgene@tudelft.nl

*Delft University of Technology*

## License

This project is licensed under the MIT License - see the project configuration for details.

## Version

Current version: v0.2.0

---

*For course-related questions and support, please contact the course instructors.*
