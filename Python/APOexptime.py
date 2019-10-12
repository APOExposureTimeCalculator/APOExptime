# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:25:10 2019

@authors: Alexander
"""

from astropy.io import ascii
import glob
from scipy import interpolate

class Instrument:
    def __init__(self, Instr_name):
        
        para = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_param.data')
        
        efficiency = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_eff.data')
        
        f = glob.glob('data/apo3_5m/'+Instr_name+ '/' + '*filter.data')
        for i in len(int(para['FilterNum'])):
            filt = ascii.read(row)
            name = row[row.find('filter')-1]+'filter'
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

data = ascii.read("test.data")
values = data["col1"]

rogers_mom = values[0]
hasan = values[1]

print(rogers_mom, hasan)
