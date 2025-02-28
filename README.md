# Power Spectrum Analysis

This project provides tools for analyzing power spectra from cosmological simulations.

## Setup

### Requirements
- Python 3.9 or later
- Conda package manager

### Installation
1. Clone this repository
2. Create the conda environment:
```bash
conda env create -f environment.yml
```
3. Activate the environment:
```bash 

```

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

## Future Features
- Multiple simulation comparison
- Theoretical predictions
- Various plotting options
- Statistical analysis tools