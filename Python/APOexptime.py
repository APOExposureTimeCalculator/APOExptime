# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:25:10 2019

@authors: Alexander
"""

from astropy.io import ascii
from scipy import interpolate

class Instrument:
    def __init__(self, Instr_name):
        
        para = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_param.data')
        
        efficiency = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_eff.data')
        
      
        for row in para['Filters']:
            filt = ascii.read('data/apo3_5m/'+Instr_name +"/"+row)
            name = row.split('.data')[0]
            setattr(Instrument, name, filt)
        
        efficiency_wavelength=efficiency["col1"]
        efficiency_percent=efficiency["col2"]/100  #divided by 100 to turn into decimal
        filt_wavelength=f["col1"]
        filt_throughput=f["col2"]
        
        efficiency_interpolated = interpolate.InterpolatedUnivariateSpline(
                efficiency_wavelength, efficiency_percent)
        filt_interpolated = interpolate.InterpolatedUnivariateSpline(
                filt_wavelength, filt_throughput)

        
        self.transmission = 
        self.efficiency = efficiency_interpolated
        self.readout_noise = 
        self.filter_wavelengths = filt_interpolated

# How to read our data files:


