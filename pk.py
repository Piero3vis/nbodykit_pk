import argparse
import numpy as np
from nbodykit.io.bigfile import BigFile
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
        bf = BigFile(input_path)[:]
        logging.info("Successfully loaded BigFile")
        
        # Convert positions from Kpc/h to Mpc/h
        positions = bf['Position'] / 1000.0  # Convert from Kpc/h to Mpc/h
        logging.info("Converted positions from Kpc/h to Mpc/h")
        
        data = {'Position': positions}
        arr = ArrayCatalog(data)
        logging.info(f"Created ArrayCatalog with {len(arr)} particles")
        
        logging.info(f"Creating mesh with Nmesh={nmesh}, BoxSize={box_size}")
        mesh = arr.to_mesh(Nmesh=nmesh, BoxSize=box_size, window=None, resampler='tsc', interlaced=True, compensated=True)
        
        logging.info("Calculating power spectrum")
        k_min = calculate_k_range(box_size, nmesh)[0]
        power_spectrum = FFTPower(mesh, mode='1d', kmin=k_min)
        power_spectrum.shot_noise = True
        k = power_spectrum.power['k']
        print(f'computed k min: {k.min()}, difference between computed and theoretical k min: {k.min() - k_min_theory}')
        pk = power_spectrum.power['power'].real
        shot_noise = power_spectrum.attrs['volume'] / power_spectrum.attrs['N1']
        
        print(f'computed dk: {power_spectrum.attrs["dk"]}')
        # Create directory for power spectrum data if it doesn't exist
        os.makedirs('power_spectrum', exist_ok=True)
        
        # Generate output filename based on parameters
        sim_name = input_path.rstrip('/').split('/')[-1]
        output_file = os.path.join('power_spectrum', 
                                 f"power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.txt")
        
        # Save to a file
        logging.info(f"Saving power spectrum to {output_file}")
        np.savetxt(output_file, np.column_stack([k, pk, np.repeat(shot_noise, len(pk))]), 
                  header="k P(k) shot_noise", fmt='%.8e')
        
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
    logging.info("Creating power spectrum plot")
    
    # Calculate k_min and k_max
    k_min_theory, k_max_theory = calculate_k_range(box_size, nmesh)
    
    # Set colors using a modified Wes Anderson palette
    measured_color = '#0D324D'    # prussian blue
    linear_color = '#3F88C5'      # orange
    nonlinear_color = '#FFBA08'    # vermilion
    vertical_color = '#999999'    # gray
    
    # Create figure with two panels
    fig = plt.figure(figsize=(12, 10))
    gs = plt.GridSpec(2, 1, height_ratios=[2, 1], hspace=0)  # Set hspace to 0
    
    # Top panel for power spectra
    ax1 = plt.subplot(gs[0])
    
    # Plot measured power spectrum with clean line
    ax1.loglog(k, pk - shot_noise, color=measured_color, 
              linestyle='-', linewidth=3, 
              label="Measured P(k)", zorder=3)
    
    # Calculate and plot both linear and non-linear CAMB power spectra
    k_lin, pk_lin, pk_nonlin = get_camb_linear_pk(box_size)
    ax1.loglog(k_lin, pk_lin, '--', color=linear_color, 
              label="Linear P(k) (WMAP9)", linewidth=3, zorder=1)
    ax1.loglog(k_lin, pk_nonlin, ':', color=nonlinear_color, 
              label="Non-linear P(k) (WMAP9)", linewidth=4, zorder=2)
    
    # Add vertical lines for k_min and k_max
    ax1.axvline(x=k_min_theory, color=vertical_color, linestyle=':', linewidth=1, alpha=0.7)
    ax1.axvline(x=k_max_theory, color=vertical_color, linestyle=':', linewidth=1, alpha=0.7)
    
    # Add text annotations for k_min and k_max
    
    
    # Bottom panel for differences
    ax2 = plt.subplot(gs[1])
    
    # Interpolate theoretical power spectra to measured k values
    f_lin = interp1d(k_lin, pk_lin, bounds_error=False, fill_value='extrapolate')
    f_nonlin = interp1d(k_lin, pk_nonlin, bounds_error=False, fill_value='extrapolate')
    
    # Calculate differences
    diff_lin = (pk - shot_noise) - f_lin(k)
    diff_nonlin = (pk - shot_noise) - f_nonlin(k)
    
    # Plot differences with log scale on x-axis only
    ax2.plot(k, diff_lin, '--', color=linear_color, 
            label="Measured - Linear", linewidth=2)
    ax2.plot(k, diff_nonlin, ':', color=nonlinear_color, 
            label="Measured - Non-linear", linewidth=3)
    ax2.plot(k, np.zeros_like(k), color=measured_color, linestyle='-', linewidth=3)
    ax2.axvline(x=k_min_theory, color=vertical_color, linestyle=':', linewidth=1, alpha=0.7)
    ax2.axvline(x=k_max_theory, color=vertical_color, linestyle=':', linewidth=1, alpha=0.7)
       
    # Set x-axis to log scale after plotting
    ax2.set_xscale('log')

    ymin = ax2.get_ylim()[0]
    ax2.text(k_min_theory*1.1   , ymin*1.5, 'k_min', rotation=90, 
            color=vertical_color, fontsize=12, alpha=0.8)
    ax2.text(k_max_theory*1.1, ymin*1.5, 'k_max', rotation=90, 
            color=vertical_color, fontsize=12, alpha=0.8)
    
    # Customize both panels
    for ax in [ax1, ax2]:
        ax.set_xlim(k_min_theory/1.5, k_max_theory*1.5)
        #ax.grid(True, which='both', linestyle=':', alpha=0.2)
        ax.tick_params(axis='both', which='major', labelsize=12, length=8, width=1)
        ax.tick_params(axis='both', which='minor', labelsize=10, length=4, width=1)
        # ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        # ax.yaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        # ax.ticklabel_format(style='sci', scilimits=(-2,2), axis='both')
    
    # Hide x-axis labels for top panel
    ax1.tick_params(labelbottom=False)
    
    # Labels and title
    ax2.set_xlabel('k [h/Mpc]', fontsize=14, labelpad=10)
    ax1.set_ylabel('P(k) [(Mpc/h)³]', fontsize=14, labelpad=10)
    ax2.set_ylabel('ΔP(k) [(Mpc/h)³]', fontsize=14, labelpad=10)
    ax1.set_title(f'Matter Power Spectrum\n(Box={box_size} Mpc/h, Nmesh={nmesh})', 
                 fontsize=16, pad=20)
    
    # Legends
    ax1.legend(fontsize=12, framealpha=0.9, loc='lower left', 
             bbox_to_anchor=(0.02, 0.02), borderaxespad=0.)
    ax2.legend(fontsize=12, framealpha=0.9, loc='upper right')
    
    # Add some padding around the plot
    plt.tight_layout()
    
    # Create plots directory if it doesn't exist
    os.makedirs('plots', exist_ok=True)
    
    # Generate plot filename based on parameters
    plot_file = os.path.join('plots', 
                           f"power_spectrum_box{box_size}_nmesh{nmesh}_{sim_name}.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logging.info(f"Plot saved to {plot_file}")
    
    if show:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Calculate and plot power spectrum from BigFile')
    parser.add_argument('input_path', help='Path to the input BigFile')
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
