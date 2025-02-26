#!/bin/bash

module purge
module load python/3.8.6--gcc--10.2.0
module load anaconda3/2020.07--gcc--8.3.1
source activate nbodykit-env

## Usage: name_program, name of the simulation, cosmology, directory path to the simulations to analyse, snapshot, number of input files of the snapshot, boxsize[Mpc/h], Nmesh (the number of cells per side on the mesh to compute the power spectrum)
python pk_matter.py ../SRS-map2map/sim/PART_010/    
