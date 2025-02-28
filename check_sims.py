import os
import sys
import logging
from bigfile import BigFile

def setup_logging():
    """Set up logging configuration."""
    formatter = logging.Formatter('%(message)s')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logging.root.handlers = []
    logging.root.addHandler(handler)
    logging.root.setLevel(logging.INFO)

def check_simulations():
    """Check for LR, HR, SR simulations and their Position data."""
    sim_types = ['LR', 'HR', 'SR']
    sim_info = {}
    
    logging.info("\nChecking for simulation directories...")
    
    # First check if data directory exists
    if not os.path.exists('data'):
        logging.error("'data' directory not found!")
        sys.exit(1)
    
    for sim_type in sim_types:
        path = os.path.join('data', sim_type)
        if os.path.exists(path):
            logging.info(f"\nFound {sim_type} simulation directory:")
            logging.info(f"Path: {path}")
            
            # Check for Position directory
            pos_path = os.path.join(path, 'Position')
            if os.path.exists(pos_path):
                logging.info("Found Position directory")
                try:
                    # Try to read Position data with BigFile
                    bf = BigFile(path)
                    pos = bf['Position']
                    shape = pos.shape if hasattr(pos, 'shape') else None
                    
                    logging.info(f"Successfully read Position data:")
                    logging.info(f"Shape: {shape}")
                    logging.info(f"Data type: {pos.dtype}")
                    
                    sim_info[sim_type] = {
                        'path': path,
                        'shape': shape,
                        'dtype': str(pos.dtype)
                    }
                except Exception as e:
                    logging.warning(f"Error reading Position data: {str(e)}")
            else:
                logging.warning(f"No Position directory found in {path}")
        else:
            logging.warning(f"{sim_type} simulation directory not found in data/")
    
    if not sim_info:
        logging.error("No valid simulation data found")
        sys.exit(1)
    
    return sim_info

def main():
    setup_logging()
    sim_info = check_simulations()
    
    # Print summary
    logging.info("\nSummary of available simulations:")
    for sim_type, info in sim_info.items():
        logging.info(f"\n{sim_type}:")
        for key, value in info.items():
            if key != 'path':
                logging.info(f"  {key}: {value}")

if __name__ == "__main__":
    main() 