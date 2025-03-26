import argparse
import numpy as np
from nbodykit.source.catalog import  ArrayCatalog
from nbodykit.lab import FFTPower
import matplotlib.pyplot as plt
import camb
from camb import model, initialpower
import os
from scipy.interpolate import interp1d
import readsnap_mod as rs


def calculate_k_range(box_size, nmesh):
    """Calculate theoretical k_min and k_max given box size and nmesh"""
    k_min = 2 * np.pi / box_size
    k_max = np.pi * nmesh / box_size
    return k_min, k_max


def calculate_power_spectrum(inpath, nmesh=256, box_size=100.0):
    print(f"Reading BigFile from: {inpath}")
    
    # Determine simulation type (LR/HR/SR) from path
    sim_type = None
    for t in ['LR', 'HR', 'SR']:
        if t in inpath:
            sim_type = t
            break
    print(f"Simulation type: {sim_type if sim_type else 'unknown'}")
    
    # Calculate and print theoretical k range
    k_min_theory, k_max_theory = calculate_k_range(box_size, nmesh)
    print(f"\nTheoretical k range formulas:")
    print(f"k_min = 2π/L     (fundamental mode)")
    print(f"k_max = πN/L     (Nyquist frequency)")
    print(f"where L = {box_size} Mpc/h (box size)")
    print(f"      N = {nmesh} (mesh cells)")
    
    print(f"\nCalculated k range for simulation:")
    print(f"k_min = {k_min_theory:.6f} h/Mpc = {k_min_theory:.2e} h/Mpc")
    print(f"     = (2π)/{box_size:.1f}")
    print(f"k_max = {k_max_theory:.6f} h/Mpc = {k_max_theory:.2e} h/Mpc")
    print(f"     = π*{nmesh}/{box_size:.1f}")
    
    try:
        # original simulation has a different format for getting the positions
        if os.path.exists(os.path.join(inpath, 'Position')):
            print("Found simple Position format")
            
        else:
            header = rs.snapshot_header(inpath+'.0')
            boxsize = header.boxsize
            redshift = header.redshift
            #Ng = header.npart[1] ** (1/3)
            #Ng = int(np.rint(Ng))
            Ng = 1024

            cellsize = boxsize / Ng
            # get parameters
            print(f'boxsize, redshift: {boxsize}, {redshift}')
            part_type=1
        #get positions and velocities
            pos_ = rs.read_block(inpath,"POS ",parttype=part_type,verbose=True)
            pid_ = rs.read_block(inpath,"ID ",parttype=part_type,verbose=True)- 1
            print(f'shape of pos_: {np.shape(pos_)}')
            pos = np.empty_like(pos_)
            pos[pid_] = pos_
            pos = pos.reshape(Ng, Ng, Ng, 3)
            positions = pos.reshape(-1, 3)
            del pid_, pos_
            
            
        
        # Convert positions from Kpc/h to Mpc/h
        positions = positions / 1000.0  # Convert from Kpc/h to Mpc/h
        print("Converted positions from Kpc/h to Mpc/h")
        
        data = {'Position': positions}
        arr = ArrayCatalog(data)
        print(f"Created ArrayCatalog with {len(arr)} particles")
        
        print(f"Creating mesh with Nmesh={nmesh}, BoxSize={box_size}")
        mesh = arr.to_mesh(Nmesh=nmesh, BoxSize=box_size, window=None, resampler='tsc', interlaced=True, compensated=True)
        
        print("Calculating power spectrum")
        k_min, k_max = calculate_k_range(box_size, nmesh)
        power_spectrum = FFTPower(mesh, mode='1d', kmin=k_min)
        power_spectrum.shot_noise = True
        k = power_spectrum.power['k']
        print(f'computed k min: {k.min()}, difference between computed and theoretical k min: {k.min() - k_min_theory}')
        pk = power_spectrum.power['power'].real
        shot_noise = power_spectrum.attrs['volume'] / power_spectrum.attrs['N1']
        
        print(f'computed dk: {power_spectrum.attrs["dk"]}')
        # Create directory for power spectrum data if it doesn't exist
        os.makedirs(outpath+'/power_spectrum', exist_ok=True)
        
        # Generate output filename based on parameters
        sim_name = inpath.rstrip('/').split('/')[-1]
        output_file = os.path.join(outpath+'/power_spectrum', 
                                 f"power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.txt")
        
        # Save to a file
        print(f"Saving power spectrum to {output_file}")
        np.savetxt(output_file, np.column_stack([k, pk, np.repeat(shot_noise, len(pk))]), 
                  header="k P(k) shot_noise", fmt='%.8e')
        
        # Print output locations more clearly
        print("\nOutput files:")
        print(f"Power spectrum data: {output_file}")
        theory_file = os.path.join('power_spectrum', f"theory_pk_box{box_size}_nmesh{nmesh}.txt")
        print(f"Theoretical P(k): {theory_file}")
        
        return k, pk, shot_noise, sim_name
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        raise
    
    
def main():
    parser = argparse.ArgumentParser(description='Calculate and plot power spectrum from BigFile')
    parser.add_argument('--inpath', help='Path to the input BigFile', )
    parser.add_argument('--outpath', help='Path to the output directory')
    parser.add_argument('--nmesh', type=int, default=256, help='Number of mesh cells (default: 256)')
    parser.add_argument('--box-size', type=float, default=100.0, 
                       help='Box size in Mpc/h (default: 100.0)')
    parser.add_argument('--no-show', action='store_true', help='Do not show the plot window')
    
    args = parser.parse_args()

    
    k, pk, shot_noise, sim_name = calculate_power_spectrum(
        args.inpath, args.outpath,
        nmesh=args.nmesh, 
        box_size=args.box_size
    )