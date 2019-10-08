# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:25:10 2019

@authors: Alexander
"""

class Instrument:
    def __init__(self, Instr_name):
        
        
        
        self.transmission = 
        self.efficiency = 
        self.readout_noise = 
        self.filter_wavelengths =

# How to read our data files:
from astropy.io import ascii
data = ascii.read("test.data")
values = data["col1"]

rogers_mom = values[0]
hasan = values[1]

print(rogers_mom, hasan)
