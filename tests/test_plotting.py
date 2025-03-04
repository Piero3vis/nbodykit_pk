import unittest
import matplotlib.pyplot as plt
import numpy as np
from src.pk import plot_power_spectrum

class TestPlotting(unittest.TestCase):
    def setUp(self):
        """Set up test data"""
        self.k = np.logspace(-2, 1, 100)
        self.pk = 1000 * self.k**(-1.5)  # Simple power law
        self.shot_noise = 100.0
        self.box_size = 100.0
        self.nmesh = 256
        self.sim_name = "test_sim"
    
    def test_plot_creation(self):
        """Test if plotting function creates and saves a file"""
        import os
        
        # Create plot
        plot_power_spectrum(
            self.k, self.pk, self.shot_noise,
            self.box_size, self.nmesh, self.sim_name,
            show=False
        )
        
        # Check if file was created
        expected_file = f'outputs/plots/power_spectrum_box{self.box_size}_nmesh{self.nmesh}_{self.sim_name}.png'
        self.assertTrue(os.path.exists(expected_file),
            f"Plot file not created: {expected_file}")
        
        # Clean up
        if os.path.exists(expected_file):
            os.remove(expected_file) 