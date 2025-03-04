
import os
import argparse
import logging
import pk
from pk import setup_logging, calculate_power_spectrum, plot_power_spectrum
#from LR_HR_SR import plot_simulation_comparison


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

    # Create output directories
    print(f"Creating output directories: outputs/plots, outputs/power_spectrum")
    os.makedirs(f"outputs/plots", exist_ok=True)
    os.makedirs(f"outputs/power_spectrum", exist_ok=True)
    
    k, pk, shot_noise, sim_name = calculate_power_spectrum(
        args.input_path, 
        nmesh=args.nmesh, 
        box_size=args.box_size
    )
    
    plot_power_spectrum(k, pk, shot_noise, args.box_size, args.nmesh, sim_name, show=not args.no_show)

    

if __name__ == "__main__":
    main() 