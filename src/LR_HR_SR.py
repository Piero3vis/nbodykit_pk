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

def calculate_k_range(box_size, nmesh):
    """Calculate theoretical k_min and k_max given box size and nmesh"""
    k_min = 2 * np.pi / box_size
    k_max = np.pi * nmesh / box_size
    return k_min, k_max

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

def plot_simulation_comparison(box_size=100.0, nmesh=256, lr_path=None, hr_path=None, sr_path=None):
    """Create a three-panel plot comparing LR, HR, and SR simulations with theory."""
    
    # Use provided paths or defaults
    lr_path = lr_path or "outputs/power_spectrum/power_spectrum_box100.0_nmesh256_lr_ds_from_sim_64.txt"
    hr_path = hr_path or "outputs/power_spectrum/power_spectrum_box100.0_nmesh256_PART_010.txt"
    sr_path = sr_path or "outputs/power_spectrum/power_spectrum_box100.0_nmesh256_sr_from_sim_64_8x_ds_32.txt"
    
    # Load measured power spectra
    k_lr, pk_lr, shot_noise_lr = load_data(lr_path)
    k_hr, pk_hr, shot_noise_hr = load_data(hr_path)
    k_sr, pk_sr, shot_noise_sr = load_data(sr_path)
    
    # Get theoretical power spectra
    k_theory, pk_lin, pk_nonlin = get_camb_linear_pk(box_size)
    
    # Create interpolation function for theoretical spectrum
    f_nonlin = interp1d(k_theory, pk_nonlin, bounds_error=False, fill_value='extrapolate')
    
    # Calculate measured P(k) minus shot noise
    pk_lr_measured = pk_lr - shot_noise_lr
    pk_hr_measured = pk_hr - shot_noise_hr
    pk_sr_measured = pk_sr - shot_noise_sr
    
    # Calculate ratios and differences with non-linear theory
    pk_nonlin_lr = f_nonlin(k_lr)
    pk_nonlin_hr = f_nonlin(k_hr)
    pk_nonlin_sr = f_nonlin(k_sr)
    
    ratio_lr = pk_lr_measured / pk_nonlin_lr
    ratio_hr = pk_hr_measured / pk_nonlin_hr
    ratio_sr = pk_sr_measured / pk_nonlin_sr
    
    diff_lr = pk_lr_measured - pk_nonlin_lr
    diff_hr = pk_hr_measured - pk_nonlin_hr
    diff_sr = pk_sr_measured - pk_nonlin_sr
    
    # Create figure with three panels
    fig = plt.figure(figsize=(10, 12))
    gs = plt.GridSpec(3, 1, height_ratios=[1.2, 1, 1], hspace=0)
    
    # Calculate k range
    k_min_theory, k_max_theory = calculate_k_range(box_size, nmesh)
    x_min = k_min_theory/1.2
    x_max = k_max_theory*1.1

    
    # Colors

    palette = ["#084c61","#ffba08","#890620","#032b43"]
    colors = {'LR': palette[0], 'HR': palette[1], 'SR': palette[2], 'Theory': palette[-1]}
    
    # Top panel: P(k)
    ax1 = plt.subplot(gs[0])
    ax1.loglog(k_lr, pk_lr_measured, '-', color=colors['LR'], label='LR (Ng = 32)', linewidth=2)
    ax1.loglog(k_hr, pk_hr_measured, '-', color=colors['HR'], label='HR (Ng = 64)', linewidth=2)
    ax1.loglog(k_sr, pk_sr_measured, '-', color=colors['SR'], label='SR (Ng = 32 (upsampling 8x, random sampled to 32))', linewidth=2)
    ax1.loglog(k_theory, pk_nonlin, '--', color=colors['Theory'], label='Non-linear theory (WMAP9)', linewidth=2)
    ax1.axvline(x=k_max_theory, color='black', linestyle='--', alpha=0.4)
    ax1.set_ylabel('P(k) [(Mpc/h)³]', fontsize=14)
    ax1.legend(fontsize=12)
    
    # Middle panel: Ratio
    ax2 = plt.subplot(gs[1])
    ax2.semilogx(k_lr, ratio_lr, '-', color=colors['LR'], label='LR/Theory', linewidth=2)
    ax2.semilogx(k_hr, ratio_hr, '-', color=colors['HR'], label='HR/Theory', linewidth=2)
    ax2.semilogx(k_sr, ratio_sr, '-', color=colors['SR'], label='SR/Theory', linewidth=2)
    ax2.axhline(y=1, color=colors['Theory'], linestyle='-', alpha=0.8)
    ax2.set_ylabel('P(k)$_{measured}$/P(k)$_{theory}$', fontsize=14)
    ax2.axvline(x=k_max_theory, color='black', linestyle='--', alpha=0.4)
    ax2.legend(fontsize=12)
    
    # Bottom panel: Difference
    ax3 = plt.subplot(gs[2])
    ax3.semilogx(k_lr, diff_lr, '-', color=colors['LR'], label='LR - Theory', linewidth=2)
    ax3.semilogx(k_hr, diff_hr, '-', color=colors['HR'], label='HR - Theory', linewidth=2)
    ax3.semilogx(k_sr, diff_sr, '-', color=colors['SR'], label='SR - Theory', linewidth=2)
    ax3.axhline(y=0, color=colors['Theory'], linestyle='-', alpha=0.8)
    ax3.axvline(x=k_max_theory, color='black', linestyle='--', alpha=0.4)
    ax3.set_ylabel('diff [(Mpc/h)³]', fontsize=14)

    ax3.set_xlabel('k [h/Mpc]', fontsize=14)
    ax3.legend(fontsize=12)
    
    # Format all panels
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim(x_min, x_max)
        ax.grid(False)
        ax.tick_params(axis='both', which='major', labelsize=12)
        
        # Add vertical line for k_nyq with annotation
        ax.axvline(x=k_max_theory, color='black', linestyle='--', alpha=0.4)
        if ax == ax1:  # Only add text annotation on top panel
            ax.text(k_max_theory*0.8, ax.get_ylim()[1]*0.5, '$k_{nyq}$', 
                   rotation=90, verticalalignment='center', fontsize=12)
        
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
    
    # Add title
    plt.suptitle(f'Matter Power Spectrum Comparison\n(Box={box_size} Mpc/h, Nmesh={nmesh})', 
                y=0.95, fontsize=16)
    
    
    os.makedirs('outputs/plots/comparison_plots', exist_ok=True)
    plt.savefig('outputs/plots/comparison_plots/power_spectrum_comparison.png', bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()




def create_comparison_plot(box_size=100.0, nmesh=256):
    # File paths
    lr_path = "power_spectrum/power_spectrum_box100.0_nmesh256_lr_ds_from_sim_64.txt"
    hr_path = "power_spectrum/power_spectrum_box100.0_nmesh256_PART_010.txt"
    sr_path = "power_spectrum/power_spectrum_box100.0_nmesh64_sr_from_sim_64_8x.txt"
    
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
        ax1.set_title(title, pad=10, fontsize=16)
        
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

def main():
    parser = argparse.ArgumentParser(description='Generate comparison plots for LR, HR, and SR power spectra')
    parser.add_argument('--box-size', type=float, default=100.0,
                       help='Box size in Mpc/h (default: 100.0)')
    parser.add_argument('--nmesh', type=int, default=256,
                       help='Number of mesh cells (default: 256)')
    parser.add_argument('--lr-path', 
                       default="outputs/power_spectrum/power_spectrum_box100.0_nmesh256_lr_ds_from_sim_64.txt",
                       help='Path to LR power spectrum file (default: outputs/power_spectrum/power_spectrum_box100.0_nmesh256_lr_ds_from_sim_64.txt)')
    parser.add_argument('--hr-path',
                       default="outputs/power_spectrum/power_spectrum_box100.0_nmesh256_PART_010.txt",
                       help='Path to HR power spectrum file (default: outputs/power_spectrum/power_spectrum_box100.0_nmesh256_PART_010.txt)')
    parser.add_argument('--sr-path',
                       default="outputs/power_spectrum/power_spectrum_box100.0_nmesh256_sr_from_sim_64_8x_ds_32.txt",
                       help='Path to SR power spectrum file (default: outputs/power_spectrum/power_spectrum_box100.0_nmesh256_sr_from_sim_64_8x_ds_32.txt)')
    
    args = parser.parse_args()
    
    plot_simulation_comparison(
        box_size=args.box_size,
        nmesh=args.nmesh,
        lr_path=args.lr_path,
        hr_path=args.hr_path,
        sr_path=args.sr_path
    )

if __name__ == "__main__":
    main()