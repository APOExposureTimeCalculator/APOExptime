# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:25:10 2019

@authors: Alexander
"""

from astropy.io import ascii
import glob

class Instrument:
    def __init__(self, Instr_name):
        
        para = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_param.data')
        
        efficiency = ascii.read('data/apo3_5m/'+Instr_name +"/"+Instr_name + '_eff.data')
        
        f = glob.glob('data/apo3_5m/'+Instr_name+ '/' + '*filter.data')
        for i in len(int(para['FilterNum'])):
            filt = ascii.read(row)
            name = row[row.find('filter')-1]+'filter'
            setattr(Instrument, name, filt)
        
        self.transmission = 
        self.efficiency = 
        self.readout_noise = 
        self.filter_wavelengths =

# How to read our data files:

data = ascii.read("test.data")
values = data["col1"]

rogers_mom = values[0]
hasan = values[1]

print(rogers_mom, hasan)
