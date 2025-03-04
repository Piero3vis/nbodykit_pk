# Power Spectrum Analysis

Tools for computing and analyzing matter power spectra from N-body simulations. The code computes power spectra from simulation data and compares them with theoretical predictions from WMAP9 cosmology using CAMB's non-linear power spectrum.

## Features
- Power spectrum computation from simulation data
- Theoretical power spectrum calculation using WMAP9 cosmology
- Comparison between different resolution simulations (LR, HR, SR)
- Interactive command-line interface
- Automated testing suite

## Project Structure
The code expects simulation data in the following structure:
```
├── data
│   ├── HR
│   │   └── PART_010
│   │       ├── 0
│   │       │   ├── Density
│   │       │   ├── EgyWtDensity
│   │       │   ├── ElectronAbundance
│   │       │   ├── Generation
│   │       │   ├── GroupID
│   │       │   ├── HeIIIIonized
│   │       │   ├── ID
│   │       │   ├── InternalEnergy
│   │       │   ├── Mass
│   │       │   ├── Metallicity
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   ├── SmoothingLength
│   │       │   └── Velocity
│   │       ├── 1
│   │       │   ├── GroupID
│   │       │   ├── ID
│   │       │   ├── Mass
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   └── Velocity
│   │       ├── 2
│   │       │   ├── GroupID
│   │       │   ├── ID
│   │       │   ├── Mass
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   └── Velocity
│   │       ├── 3
│   │       │   ├── GroupID
│   │       │   ├── ID
│   │       │   ├── Mass
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   └── Velocity
│   │       ├── 4
│   │       │   ├── BirthDensity
│   │       │   ├── Generation
│   │       │   ├── GroupID
│   │       │   ├── ID
│   │       │   ├── Mass
│   │       │   ├── Metallicity
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   ├── StarFormationTime
│   │       │   └── Velocity
│   │       ├── 5
│   │       │   ├── BlackholeAccretionRate
│   │       │   ├── BlackholeDensity
│   │       │   ├── BlackholeJumpToMinPot
│   │       │   ├── BlackholeMass
│   │       │   ├── BlackholeMinPotPos
│   │       │   ├── BlackholeMseed
│   │       │   ├── BlackholeMtrack
│   │       │   ├── BlackholeProgenitors
│   │       │   ├── BlackholeSwallowID
│   │       │   ├── BlackholeSwallowTime
│   │       │   ├── Generation
│   │       │   ├── GroupID
│   │       │   ├── ID
│   │       │   ├── Mass
│   │       │   ├── Position
│   │       │   ├── Potential
│   │       │   ├── SmoothingLength
│   │       │   ├── StarFormationTime
│   │       │   ├── Swallowed
│   │       │   └── Velocity
│   │       └── Header
│   ├── LR
│   │   ├── dis_32_ds
│   │   │   ├── Position
│   │   │   └── Velocity
│   │   └── lr_ds_from_sim_64
│   │       ├── Position
│   │       └── Velocity
│   └── SR
│       ├── sr_from_sim_32_x2
│       │   ├── Position
│       │   └── Velocity
│       └── srx2_from_dis_32_ds
│           ├── Position
│           └── Velocity
├── docs
├── outputs
│   ├── plots
│   └── power_spectrum
├── power_spectrum
├── src
│   └── __pycache__
└── tests
    └── __pycache__
```
## Setup

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/power-spectrum-analysis.git
cd power-spectrum-analysis
```

2. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate power_spectrum_env
```

3. Verify installation:
```bash
python -m unittest discover tests
```

## Usage

### Interactive Interface

Run the interactive analysis tool:
```bash
./compute_plot_pk.sh
```

This presents a menu with three options:

1. **Run Environment Check**
   - Verifies that the simulation data is properly organized
   - No additional input needed

2. **Compute Power Spectrum**
   - Interactive selection of simulation directory
   - Parameters:
     - Simulation path (interactive selection)
     - Number of mesh cells (default: 256)
     - Box size in Mpc/h (default: 100.0)
   - Outputs:
     - Power spectrum data in `outputs/power_spectum/`
     - Plot in `outputs/plots/`
   - Compares with WMAP9 non-linear power spectrum

3. **Generate Comparison Plots**
   - Creates comparison plots for LR, HR, and SR simulations
   - Currently uses fixed paths in data/LR/, data/HR/, data/SR/
   - Parameters:
     - Box size in Mpc/h (default: 100.0)
     - Number of mesh cells (default: 256)

### Individual Scripts

You can also run the scripts individually:

```bash
# Environment check
python src/check_sims.py

# Single simulation analysis
python src/main.py --input_path <path> --nmesh <N> --box-size <L>

# Multi-simulation comparison
python src/LR_HR_SR.py --box-size <L> --nmesh <N>
```

## Power Spectrum Analysis

### Theoretical Model
- Uses CAMB for theoretical predictions
- WMAP9 cosmology parameters:
  - h = 0.697
  - Ωm = 0.2814
  - Ωb = 0.0464
  - ns = 0.971
  - σ8 = 0.82
- Both linear and non-linear power spectra are computed
- Current comparison uses the non-linear prediction

### Output Analysis
Each power spectrum plot shows:
1. Top panel: P(k) measurements and theoretical predictions
2. Middle panel: Ratio of measured to theoretical power
3. Bottom panel: Absolute difference between measurement and theory

## Testing

Run the test suite:
```bash
python -m unittest discover tests
```

Tests include:
- Power spectrum calculations
- k-range validations
- Plot generation
- File operations
- Physics sanity checks

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Future Features
- Support for additional theoretical models
- Multiple redshift analysis
- Statistical analysis tools
- Custom cosmology parameters
- Additional plotting options

## Data Structure
The code expects simulation data in the following structure:
```
data/
├── LR/  # Low Resolution simulation
│   └── Position/
├── HR/  # High Resolution simulation
│   └── Position/
└── SR/  # Standard Resolution simulation
    └── Position/
```

## Data Requirements

### Resolution Requirements
- HR and SR simulations must have the same number of particles
- LR simulation should have a number of particles that divides evenly into HR/SR

### Position Requirements
- Positions are in kpc units
- Valid range: [0, 100000] kpc
- **Important**: Positions must be divided by 1000 for Pk calculations (to convert to Mpc)

## Scripts

### check_sims.py
Preliminary check of simulation data:
- Verifies directory structure
- Supports both simple and header formats
- Checks particle counts and cube roots
- Validates position ranges
- Verifies resolution relationships
- Reports simulation metadata (boxsize, redshift, etc. when available)

Usage:
```bash
python check_sims.py
```

### pk.py (coming soon)
Power spectrum calculation and analysis:
- Computes power spectra for each simulation
- Compares with theoretical predictions
- Generates comparison plots

## File Descriptions

### environment.yml
Conda environment specification:
- Python dependencies
- Required packages for power spectrum analysis
- Plotting libraries

Key packages:
- nbodykit: For power spectrum calculations
- bigfile: For reading simulation data
- camb: For theoretical power spectra
- matplotlib: For plotting
- numpy: For numerical operations

## Usage Notes

1. Ensure your simulation data is properly organized in the `data/` directory
2. Run `check_sims.py` first to verify data accessibility
3. Use `pk.py` (when available) for power spectrum analysis

## Troubleshooting

Common issues:
1. Missing `data/` directory
   - Create the directory and ensure proper simulation data structure
2. BigFile read errors
   - Check file permissions
   - Verify data format matches expected structure
3. Resolution mismatches
   - Verify HR and SR have identical particle counts
   - Verify LR particle count divides evenly into HR/SR
4. Position range issues
   - Verify positions are in kpc
   - Check for values outside [0, 100000] range