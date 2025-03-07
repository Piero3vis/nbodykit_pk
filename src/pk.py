import argparse
import numpy as np
from nbodykit.io.bigfile import BigFile
from bigfile import File
from nbodykit.source.catalog import BigFileCatalog, ArrayCatalog
from nbodykit.lab import FFTPower
import matplotlib.pyplot as plt
import logging
import camb
from camb import model, initialpower
import os
from scipy.interpolate import interp1d

def setup_logging():
    # Create a custom formatter that only shows the message
    formatter = logging.Formatter('%(message)s')
    
    # Set up the handler with the custom formatter
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    
    # Configure the root logger
    logging.root.handlers = []
    logging.root.addHandler(handler)
    logging.root.setLevel(logging.INFO)

def calculate_k_range(box_size, nmesh):
    """Calculate theoretical k_min and k_max given box size and nmesh"""
    k_min = 2 * np.pi / box_size
    k_max = np.pi * nmesh / box_size
    return k_min, k_max

def calculate_power_spectrum(input_path, nmesh=256, box_size=100.0):
    logging.info(f"Reading BigFile from: {input_path}")
    
    # Determine simulation type (LR/HR/SR) from path
    sim_type = None
    for t in ['LR', 'HR', 'SR']:
        if t in input_path:
            sim_type = t
            break
    logging.info(f"Simulation type: {sim_type if sim_type else 'unknown'}")
    
    # Calculate and print theoretical k range
    k_min_theory, k_max_theory = calculate_k_range(box_size, nmesh)
    logging.info(f"\nTheoretical k range formulas:")
    logging.info(f"k_min = 2π/L     (fundamental mode)")
    logging.info(f"k_max = πN/L     (Nyquist frequency)")
    logging.info(f"where L = {box_size} Mpc/h (box size)")
    logging.info(f"      N = {nmesh} (mesh cells)")
    
    logging.info(f"\nCalculated k range for simulation:")
    logging.info(f"k_min = {k_min_theory:.6f} h/Mpc = {k_min_theory:.2e} h/Mpc")
    logging.info(f"     = (2π)/{box_size:.1f}")
    logging.info(f"k_max = {k_max_theory:.6f} h/Mpc = {k_max_theory:.2e} h/Mpc")
    logging.info(f"     = π*{nmesh}/{box_size:.1f}")
    
    try:
        # original simulation has a different format for getting the positions
        if os.path.exists(os.path.join(input_path, 'Position')):
            logging.info("Found simple Position format")
            bf = BigFile(input_path)[:]
            positions = bf['Position']
        else:
            logging.info("Found Header format")
            bigf = File(input_path)
            header = bigf.open('Header')
            boxsize = header.attrs['BoxSize'][0]
            redshift = 1./header.attrs['Time'][0] - 1
            
            Ng = header.attrs['TotNumPart'][1] ** (1/3)
            Ng = int(np.rint(Ng))

            cellsize = boxsize / Ng

            pid_ = bigf.open('1/ID')[:] - 1   # so that particle id starts from 0
            pos_ = bigf.open('1/Position')[:]
            pos = np.empty_like(pos_)
            pos[pid_] = pos_
            pos = pos.reshape(Ng, Ng, Ng, 3)
            positions = pos.reshape(-1, 3)
            
            logging.info("Successfully loaded BigFile")
        
        # Convert positions from Kpc/h to Mpc/h
        positions = positions / 1000.0  # Convert from Kpc/h to Mpc/h
        logging.info("Converted positions from Kpc/h to Mpc/h")
        
        data = {'Position': positions}
        arr = ArrayCatalog(data)
        logging.info(f"Created ArrayCatalog with {len(arr)} particles")
        
        logging.info(f"Creating mesh with Nmesh={nmesh}, BoxSize={box_size}")
        mesh = arr.to_mesh(Nmesh=nmesh, BoxSize=box_size, window=None, resampler='tsc', interlaced=True, compensated=True)
        
        logging.info("Calculating power spectrum")
        k_min, k_max = calculate_k_range(box_size, nmesh)
        power_spectrum = FFTPower(mesh, mode='1d', kmin=k_min)
        power_spectrum.shot_noise = True
        k = power_spectrum.power['k']
        print(f'computed k min: {k.min()}, difference between computed and theoretical k min: {k.min() - k_min_theory}')
        pk = power_spectrum.power['power'].real
        shot_noise = power_spectrum.attrs['volume'] / power_spectrum.attrs['N1']
        
        print(f'computed dk: {power_spectrum.attrs["dk"]}')
        # Create directory for power spectrum data if it doesn't exist
        os.makedirs('outputs/power_spectrum', exist_ok=True)
        
        # Generate output filename based on parameters
        sim_name = input_path.rstrip('/').split('/')[-1]
        output_file = os.path.join('outputs/power_spectrum', 
                                 f"power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.txt")
        
        # Save to a file
        logging.info(f"Saving power spectrum to {output_file}")
        np.savetxt(output_file, np.column_stack([k, pk, np.repeat(shot_noise, len(pk))]), 
                  header="k P(k) shot_noise", fmt='%.8e')
        
        # Print output locations more clearly
        logging.info("\nOutput files:")
        logging.info(f"Power spectrum data: {output_file}")
        theory_file = os.path.join('power_spectrum', f"theory_pk_box{box_size}_nmesh{nmesh}.txt")
        logging.info(f"Theoretical P(k): {theory_file}")
        
        return k, pk, shot_noise, sim_name
        
    except Exception as e:
        logging.error(f"Error occurred: {str(e)}")
        raise

def get_camb_linear_pk(box_size, kmin=1e-4, kmax=10.0, npoints=4000):
    """Calculate linear and non-linear matter power spectrum using CAMB with WMAP9 cosmology"""
    logging.info("Calculating linear and non-linear matter power spectrum with CAMB (WMAP9 cosmology)")
    logging.info(f"CAMB k range:")
    logging.info(f"k_min = {kmin:.6f} h/Mpc = {kmin:.2e} h/Mpc")
    logging.info(f"k_max = {kmax:.6f} h/Mpc = {kmax:.2e} h/Mpc")
    
    # Initial amplitude guess
    As = 2e-9
    
    # Set up parameters with WMAP9 cosmology from the paper
    h = 0.697       # Hubble parameter
    Omega_m = 0.2814  # Matter density
    Omega_b = 0.0464  # Baryon density
    
    # Convert to CAMB format
    ombh2 = Omega_b * h**2    # Physical baryon density
    omch2 = (Omega_m - Omega_b) * h**2  # Physical cold dark matter density
    
    # Set initial parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=h*100,     # H0 = 100h km/s/Mpc
        ombh2=ombh2,  # Physical baryon density
        omch2=omch2,  # Physical cold dark matter density
    )
    pars.InitPower.set_params(
        ns=0.971,     # Scalar spectral index from paper
        As=As         # Initial amplitude guess
    )
    
    # First get linear spectrum
    pars.NonLinear = model.NonLinear_none
    
    # Calculate power spectrum
    pars.set_matter_power(redshifts=[0.0], kmax=kmax)
    
    # Get results and calculate fiducial sigma8
    results = camb.get_results(pars)
    sigma8_fid = results.get_sigma8_0()
    
    # Set correct As to match desired sigma8=0.82 from paper
    sigma8_target = 0.82
    pars.InitPower.set_params(
        As=As * (sigma8_target**2 / sigma8_fid**2),
        ns=0.971
    )
    
    # Get linear spectrum first
    results = camb.get_results(pars)
    k_h, z, pk_lin = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=1000)
    
    # Now get non-linear spectrum
    pars.NonLinear = 1
    results = camb.get_results(pars)
    k_h, z, pk_nonlin = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=1000)
    
    return k_h, pk_lin[0], pk_nonlin[0]  # Return k and both P(k) at z=0

def plot_power_spectrum(k, pk, shot_noise, box_size, nmesh, sim_name, show=True):
    """Create a three-panel plot: P(k), ratio, and difference."""
    
    # Calculate theoretical spectra
    k_theory, pk_lin, pk_nonlin = get_camb_linear_pk(box_size)
    
    # Calculate measured P(k) minus shot noise
    pk_measured = pk - shot_noise
    
    # Create interpolation functions for theoretical spectra
    f_lin = interp1d(k_theory, pk_lin, bounds_error=False, fill_value='extrapolate')
    f_nonlin = interp1d(k_theory, pk_nonlin, bounds_error=False, fill_value='extrapolate')
    
    # Calculate ratios and differences
    pk_lin_interp = f_lin(k)
    pk_nonlin_interp = f_nonlin(k)
    
    ratio_lin = pk_measured / pk_lin_interp
    ratio_nonlin = pk_measured / pk_nonlin_interp
    
    diff_lin = pk_measured - pk_lin_interp
    diff_nonlin = pk_measured - pk_nonlin_interp
    
    # Set up the figure
    fig = plt.figure(figsize=(10, 12))
    gs = plt.GridSpec(3, 1, height_ratios=[1.2, 1, 1], hspace=0)
    
    # Calculate k range from theory
    k_min_theory, k_max_theory = calculate_k_range(box_size, nmesh)
    x_min = k_min_theory/1.2
    x_max = k_max_theory*1.1
    
    # Colors
    measured_color = 'navy'
    linear_color = 'orange'
    nonlinear_color = 'brown'
    
    # Top panel: P(k)
    ax1 = plt.subplot(gs[0])
    ax1.loglog(k, pk_measured, '-', color=measured_color, label='Measured', linewidth=2)
    ax1.loglog(k_theory, pk_lin, '-', color=linear_color, label='Linear (WMAP9)', linewidth=2)
    ax1.loglog(k_theory, pk_nonlin, '-', color=nonlinear_color, label='Non-linear (WMAP9)', linewidth=2)
    ax1.set_ylabel('P(k) [(Mpc/h)³]')
    ax1.legend(fontsize=10)
    ax1.set_title(f'Matter Power Spectrum\n(Box={box_size} Mpc/h, Nmesh={nmesh})', pad=10)
    
    # Middle panel: Ratio
    ax2 = plt.subplot(gs[1])
    ax2.semilogx(k, ratio_lin, '-', color=linear_color, label='Measured/Linear', linewidth=2, alpha=0.8)
    ax2.semilogx(k, ratio_nonlin, '-', color=nonlinear_color, label='Measured/Non-linear', linewidth=2, alpha=0.8)
    ax2.axhline(y=1, color='black', linestyle='-', alpha=0.8, linewidth=2)
    ax2.set_ylabel('P(k)$_{measured}$/P(k)$_{theory}$')
    ax2.legend(fontsize=10)
    
    # Bottom panel: Difference
    ax3 = plt.subplot(gs[2])
    ax3.semilogx(k, diff_lin, '-', color=linear_color, label='Measured - Linear', linewidth=2, alpha=0.8)
    ax3.semilogx(k, diff_nonlin, '-', color=nonlinear_color, label='Measured - Non-linear', linewidth=2, alpha=0.8)
    ax3.axhline(y=0, color='black', linestyle='-', alpha=0.8, linewidth=2)
    ax3.set_ylabel('ΔP(k) [(Mpc/h)³]')
    ax3.set_xlabel('k [h/Mpc]')
    ax3.legend(fontsize=10)
    
    # Set consistent x-axis limits and format all panels
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(x_min, x_max)
        ax.grid(False)
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Format x-axis ticks as powers of 10
        ax.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, _: r'$10^{{{:d}}}$'.format(int(np.log10(x))) if x > 0 else '0'))
        
        # Format y-axis ticks
        if ax == ax1:  # Power spectrum (log scale)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(
                lambda y, _: r'$10^{{{:d}}}$'.format(int(np.log10(y))) if y > 0 else '0'))
    
    # Hide x-axis labels for top and middle panels
    ax1.tick_params(labelbottom=False)
    ax2.tick_params(labelbottom=False)
    
    # Save and show plot
    os.makedirs('outputs/plots', exist_ok=True)
    plt.savefig(f'outputs/plots/power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.png', 
                bbox_inches='tight', dpi=300)
    print(f'saved plot to outputs/plots/power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.png')
    if show:
        plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Calculate and plot power spectrum from BigFile')
    parser.add_argument('--input_path', help='Path to the input BigFile')
    parser.add_argument('--nmesh', type=int, default=256, help='Number of mesh cells (default: 256)')
    parser.add_argument('--box-size', type=float, default=100.0, 
                       help='Box size in Mpc/h (default: 100.0)')
    parser.add_argument('--no-show', action='store_true', help='Do not show the plot window')
    
    args = parser.parse_args()
    
    setup_logging()
    
    logging.info(f"Configuration:")
    logging.info(f"Box size: {args.box_size} Mpc/h")
    logging.info(f"Nmesh: {args.nmesh}")
    
    k, pk, shot_noise, sim_name = calculate_power_spectrum(
        args.input_path, 
        nmesh=args.nmesh, 
        box_size=args.box_size
    )
    
    plot_power_spectrum(k, pk, shot_noise, args.box_size, args.nmesh, sim_name, show=not args.no_show)

if __name__ == "__main__":
    main()