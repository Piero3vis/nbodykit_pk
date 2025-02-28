import os
import sys
import logging
import numpy as np
from bigfile import BigFile, File

def setup_logging():
    """Set up logging configuration."""
    formatter = logging.Formatter('%(message)s')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logging.root.handlers = []
    logging.root.addHandler(handler)
    logging.root.setLevel(logging.INFO)

def is_cube(n):
    """Check if a number is a perfect cube and return its cube root."""
    cube_root = n**(1./3.)
    is_perfect_cube = abs(round(cube_root) ** 3 - n) < 1e-10
    if is_perfect_cube:
        logging.info(f"Particle number {n:,} is a perfect cube!")
        logging.info(f"Cube root is {round(cube_root):,} (i.e., {round(cube_root):,}³ particles)")
        return round(cube_root)
    else:
        logging.info(f"Particle number {n:,} is not a perfect cube")
        logging.info(f"Cube root would be approximately {cube_root:.2f}")
        return None

def read_simple_format(path):
    """Read data from simple Position/Velocity format."""
    bf = BigFile(path)
    positions = bf['Position'][:]
    return positions, None, None, None

def read_header_format(path):
    """Read data from format with Header and numbered subgroups."""
    bigf = File(path)
    header = bigf.open('Header')
    boxsize = header.attrs['BoxSize'][0]
    redshift = 1./header.attrs['Time'][0] - 1
    
    Ng = header.attrs['TotNumPart'][1] ** (1/3)
    Ng = int(np.rint(Ng))
    
    cellsize = boxsize / Ng
    
    positions = bigf.open('1/Position')[:]
    return positions, boxsize, redshift, cellsize

def check_simulations():
    """Check for LR, HR, SR simulations and their data."""
    sim_types = ['LR', 'HR', 'SR']
    sim_info = {}
    
    logging.info("\nChecking for simulation directories...")
    
    if not os.path.exists('data'):
        logging.error("'data' directory not found!")
        sys.exit(1)
    
    for sim_type in sim_types:
        base_path = os.path.join('data', sim_type)
        if os.path.exists(base_path):
            logging.info(f"\nFound {sim_type} directory:")
            logging.info(f"Base path: {base_path}")
            
            for subfolder in os.listdir(base_path):
                full_path = os.path.join(base_path, subfolder)
                if os.path.isdir(full_path):
                    logging.info(f"\nChecking subfolder: {subfolder}")
                    
                    try:
                        # Try both formats
                        if os.path.exists(os.path.join(full_path, 'Position')):
                            logging.info("Found simple Position format")
                            positions, boxsize, redshift, cellsize = read_simple_format(full_path)
                            format_type = "simple"
                        else:
                            logging.info("Trying Header format")
                            positions, boxsize, redshift, cellsize = read_header_format(full_path)
                            format_type = "header"
                        
                        logging.info("Successfully read Position data:")
                        logging.info(f"Array shape: {positions.shape}")
                        
                        # Check if number of particles is a perfect cube
                        logging.info("\nChecking particle count...")
                        n_particles = positions.shape[0]
                        cube_root = is_cube(n_particles)
                        
                        logging.info(f"\nData type: {positions.dtype}")
                        logging.info(f"Memory size: {positions.nbytes / 1e9:.2f} GB")
                        
                        # Basic statistics
                        logging.info("\nPosition ranges:")
                        for i, axis in enumerate(['x', 'y', 'z']):
                            min_val = positions[:,i].min()
                            max_val = positions[:,i].max()
                            mean_val = positions[:,i].mean()
                            logging.info(f"{axis}: [{min_val:.2f}, {max_val:.2f}], mean: {mean_val:.2f}")
                        
                        sim_info_dict = {
                            'path': full_path,
                            'format': format_type,
                            'shape': positions.shape,
                            'n_particles': n_particles,
                            'cube_root': cube_root,
                            'dtype': str(positions.dtype),
                            'size_gb': positions.nbytes / 1e9,
                            'ranges': {
                                'x': [float(positions[:,0].min()), float(positions[:,0].max())],
                                'y': [float(positions[:,1].min()), float(positions[:,1].max())],
                                'z': [float(positions[:,2].min()), float(positions[:,2].max())]
                            }
                        }
                        
                        # Add Header information if available
                        if format_type == "header":
                            sim_info_dict.update({
                                'boxsize': boxsize,
                                'redshift': redshift,
                                'cellsize': cellsize
                            })
                        
                        sim_key = f"{sim_type}/{subfolder}"
                        sim_info[sim_key] = sim_info_dict
                        
                    except Exception as e:
                        logging.warning(f"Error reading simulation data: {str(e)}")
        else:
            logging.warning(f"{sim_type} directory not found in data/")
    
    if not sim_info:
        logging.error("No valid simulation data found")
        sys.exit(1)
    
    return sim_info

def check_shapes_and_ranges(sim_info):
    """Check relationships between shapes and validate ranges."""
    logging.info("\nValidating shapes and ranges...")
    
    # Get shapes for each resolution
    hr_shape = None
    sr_shape = None
    lr_shape = None
    
    for sim_key, info in sim_info.items():
        if sim_key.startswith('HR'):
            hr_shape = info['n_particles']
        elif sim_key.startswith('SR'):
            sr_shape = info['n_particles']
        elif sim_key.startswith('LR'):
            lr_shape = info['n_particles']
    
    # Check shapes
    if hr_shape and sr_shape:
        if hr_shape == sr_shape:
            logging.info(f"✓ HR and SR have same number of particles: {hr_shape:,}")
        else:
            logging.warning(f"! HR ({hr_shape:,}) and SR ({sr_shape:,}) have different numbers of particles")
    
    if lr_shape and hr_shape:
        ratio = hr_shape / lr_shape
        if ratio.is_integer():
            # Check if the ratio is a perfect cube
            ratio_cube_root = round(ratio**(1./3.))
            if abs(ratio_cube_root**3 - ratio) < 1e-10:
                logging.info(f"✓ HR/SR has {ratio:.0f}x more particles than LR")
                logging.info(f"  LR: {lr_shape:,} particles (Ng_LR)")
                logging.info(f"  HR: {hr_shape:,} particles = {ratio:.0f} × LR")
                logging.info(f"  Ratio is {ratio_cube_root}³ = {ratio:.0f}")
                logging.info(f"  Therefore: Ng_HR = Ng_SR = {ratio_cube_root} × Ng_LR")
            else:
                logging.warning(f"! Ratio {ratio:.0f} is not a perfect cube!")
                logging.warning(f"  This means the resolution increase is not uniform across dimensions")
        else:
            logging.warning(f"! HR/SR particle count ({hr_shape:,}) is not a multiple of LR ({lr_shape:,})")
    
    # Check ranges (in kpc)
    logging.info("\nChecking position ranges (in kpc)...")
    for sim_key, info in sim_info.items():
        logging.info(f"\n{sim_key}:")
        for axis, (min_val, max_val) in info['ranges'].items():
            if 0 <= min_val <= 100000 and 0 <= max_val <= 100000:
                logging.info(f"✓ {axis}: [{min_val:.2f}, {max_val:.2f}] kpc - in valid range")
            else:
                logging.warning(f"! {axis}: [{min_val:.2f}, {max_val:.2f}] kpc - outside [0, 100000]")
    
    logging.info("\nIMPORTANT REMINDER:")
    logging.info("Position data is in kpc and needs to be divided by 1000 for Pk calculations (Mpc)!")

def main():
    setup_logging()
    sim_info = check_simulations()
    
    # Print summary
    logging.info("\nSummary of available simulations:")
    for sim_key, info in sim_info.items():
        logging.info(f"\n{sim_key}:")
        logging.info(f"  Format: {info['format']}")
        logging.info(f"  Shape: {info['shape']} ({info['size_gb']:.2f} GB)")
        if info['cube_root']:
            logging.info(f"  Particles: {info['n_particles']:,} = {info['cube_root']:,}³")
        else:
            logging.info(f"  Particles: {info['n_particles']:,} (not a perfect cube)")
        
        if info['format'] == 'header':
            logging.info(f"  Box size: {info['boxsize']}")
            logging.info(f"  Redshift: {info['redshift']:.2f}")
            logging.info(f"  Cell size: {info['cellsize']:.6f}")
        
        logging.info(f"  Data type: {info['dtype']}")
        logging.info("  Position ranges:")
        for axis, (min_val, max_val) in info['ranges'].items():
            logging.info(f"    {axis}: [{min_val:.2f}, {max_val:.2f}]")

    check_shapes_and_ranges(sim_info)

if __name__ == "__main__":
    main() 