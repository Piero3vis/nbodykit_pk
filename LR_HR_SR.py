import matplotlib.pyplot as plt
import numpy as np
import camb
from camb import model, initialpower
import os
import argparse

from scipy.interpolate import interp1d

def load_data(path):
    """Load power spectrum data from file."""
    data = np.loadtxt(path)
    k = data[:, 0]  # First column is k
    pk = data[:, 1]
    shot_noise = data[:, 2]
    return k, pk, shot_noise

def get_camb_linear_pk(box_size, kmin=1e-4, kmax=10.0, npoints=4000):
    """Calculate linear and non-linear matter power spectrum using CAMB with WMAP9 cosmology"""
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
    
    return k_h, pk_lin[0], pk_nonlin[0] 

def create_comparison_plot(box_size=100.0, nmesh=256):
    # File paths
    lr_path = "power_spectrum/power_spectrum_box100.0_nmesh256_lr_ds_from_sim_64.txt"
    hr_path = "power_spectrum/power_spectrum_box100.0_nmesh256_PART_010.txt"
    sr_path = "power_spectrum/power_spectrum_box100.0_nmesh256_sr_from_sim_32_x2.txt"
    
    # Load measured power spectra
    k_lr, pk_lr, shot_noise_lr = load_data(lr_path)
    k_hr, pk_hr, shot_noise_hr = load_data(hr_path)
    k_sr, pk_sr, shot_noise_sr = load_data(sr_path)
    
    # Get theoretical power spectra
    k_theory, pk_lin, pk_nonlin = get_camb_linear_pk(box_size)
    
    # Create interpolation functions for theoretical spectra
    pk_lin_interp = interp1d(k_theory, pk_lin, bounds_error=False, fill_value='extrapolate')
    pk_nonlin_interp = interp1d(k_theory, pk_nonlin, bounds_error=False, fill_value='extrapolate')
    
    # Create figure with 3x3 grid, more rectangular shape
    fig = plt.figure(figsize=(15, 10))
    gs = plt.GridSpec(3, 3, hspace=0, wspace=0.3, height_ratios=[1.5, 1, 1])
    
    # Calculate k limits based on data
    k_min = min(k_lr[0], k_hr[0], k_sr[0])/1.2
    k_max = max(k_lr[-1], k_hr[-1], k_sr[-1])*1.1
    
    # Lists for easy iteration
    k_data = [k_lr, k_hr, k_sr]
    pk_data = [pk_lr, pk_hr, pk_sr]
    shot_noise_data = [shot_noise_lr, shot_noise_hr, shot_noise_sr]
    measured_pk = [pk_data[i] - shot_noise_data[i] for i in range(len(pk_data))]
    titles = ['LR', 'HR', 'SR']
    colors = ['blue', 'red', 'green']
    
    # Create all 9 panels
    for i, (k, pk, title, color) in enumerate(zip(k_data, measured_pk, titles, colors)):
        # Top row: P(k)
        ax1 = plt.subplot(gs[0, i])
        ax1.loglog(k, pk, color=color, label=f'Measured {title}')
        ax1.loglog(k_theory, pk_lin, '--', color='orange', label='Linear Theory')
        ax1.loglog(k_theory, pk_nonlin, ':', color='brown', label='Non-linear Theory')
        ax1.set_xlim(k_min, k_max)
        
        # Add column titles
        ax1.set_title(title, pad=10, fontsize=14)
        
        if i == 0:
            ax1.set_ylabel('P(k) [(Mpc/h)³]')
        if i == 1:
            ax1.text(0.5, 1.15, f'Matter Power Spectrum\n(Box={box_size} Mpc/h, Nmesh={nmesh})', 
                    horizontalalignment='center', transform=ax1.transAxes)
        if i == 0:
            ax1.legend()
        
        # Middle row: ΔP(k)
        ax2 = plt.subplot(gs[1, i])
        pk_lin_interp_vals = pk_lin_interp(k)
        pk_nonlin_interp_vals = pk_nonlin_interp(k)
        ax2.semilogx(k, pk - pk_lin_interp_vals, '--', color='orange', label='Measured - Linear')
        ax2.semilogx(k, pk - pk_nonlin_interp_vals, ':', color='brown', label='Measured - Non-linear')
        ax2.set_xlim(k_min, k_max)
        ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        if i == 0:
            ax2.set_ylabel('ΔP(k) [(Mpc/h)³]')
        if i == 0:
            ax2.legend()
        
        # Bottom row: P(k)_measured/P(k)_theory
        ax3 = plt.subplot(gs[2, i])
        ax3.semilogx(k, pk/pk_lin_interp_vals, '--', color='orange', label='Measured/Linear')
        ax3.semilogx(k, pk/pk_nonlin_interp_vals, ':', color='brown', label='Measured/Non-linear')
        ax3.set_xlim(k_min, k_max)
        ax3.axhline(y=1, color='black', linestyle='-', alpha=0.3)
        if i == 0:
            ax3.set_ylabel('P(k)_measured/P(k)_theory')
        if i == 0:
            ax3.legend()
        
        # Only show x-axis labels on bottom row
        if ax3 != plt.gca():
            plt.setp(plt.gca().get_xticklabels(), visible=False)
        else:
            ax3.set_xlabel('k [h/Mpc]')
        
        # Remove grid
        for ax in [ax1, ax2, ax3]:
            ax.grid(False)
    
    plt.savefig('power_spectrum_comparison.png', bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    create_comparison_plot()