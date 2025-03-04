import unittest
import numpy as np
from src.pk import calculate_k_range, get_camb_linear_pk

class TestPowerSpectrum(unittest.TestCase):
    def setUp(self):
        """Set up test parameters"""
        self.box_size = 100.0
        self.nmesh = 256
    
    def test_k_range_calculation(self):
        """Test if k_min and k_max are calculated correctly"""
        k_min, k_max = calculate_k_range(self.box_size, self.nmesh)
        
        # Test k_min = 2π/L
        expected_k_min = 2 * np.pi / self.box_size
        self.assertAlmostEqual(k_min, expected_k_min, 
            msg="k_min calculation incorrect")
        
        # Test k_max = πN/L
        expected_k_max = np.pi * self.nmesh / self.box_size
        self.assertAlmostEqual(k_max, expected_k_max, 
            msg="k_max calculation incorrect")
    
    def test_camb_power_spectrum(self):
        """Test if CAMB power spectrum calculation returns valid results"""
        k, pk_lin, pk_nonlin = get_camb_linear_pk(self.box_size)
        
        # Test that arrays are not empty
        self.assertTrue(len(k) > 0, "k array is empty")
        self.assertTrue(len(pk_lin) > 0, "linear P(k) array is empty")
        self.assertTrue(len(pk_nonlin) > 0, "non-linear P(k) array is empty")
        
        # Test that values are positive
        self.assertTrue(np.all(k > 0), "k values should be positive")
        self.assertTrue(np.all(pk_lin > 0), "linear P(k) values should be positive")
        self.assertTrue(np.all(pk_nonlin > 0), "non-linear P(k) values should be positive")
        
        # Test that non-linear power is greater than linear at small scales (large k)
        high_k_mask = k > 0.1  # k > 0.1 h/Mpc
        self.assertTrue(np.all(pk_nonlin[high_k_mask] >= pk_lin[high_k_mask]), 
            "Non-linear power should be >= linear power at small scales")

class TestFileOperations(unittest.TestCase):
    def test_output_directories(self):
        """Test if output directories exist"""
        import os
        
        # Test plots directory
        self.assertTrue(os.path.exists("outputs/plots"), 
            "Plots directory does not exist")
        
        # Test power spectra directory
        self.assertTrue(os.path.exists("outputs/power_spectrum"), 
            "Power spectra directory does not exist")

if __name__ == '__main__':
    unittest.main() 