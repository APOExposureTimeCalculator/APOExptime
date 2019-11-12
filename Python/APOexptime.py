# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:25:10 2019

@authors: Alexander, Bryson, Roger, Hassan, Manny
"""

from astropy.io import ascii
from scipy import interpolate

class Instrument:
    """Instantiates an object of the Instrument class.
    
    Args:
        Instr_name(str): Name of Instrument to be used.
        
    Attributes:
        efficiency(object): UnivariateInterpolatedSpline of instrument efficiency.
        readout_noise(float): Value of instrument readout noise
        filter_num(int): Number of filters for instrument
        gain(float): Gain of instrument
        """
    
    def __init__(self, Instr_name):
        
        para = ascii.read('../data/apo3_5m/'+Instr_name +"/"+Instr_name + '_param.data')
        
        efficiency = ascii.read('../data/apo3_5m/'+Instr_name +"/"+Instr_name + '_qe.data')
        
      
        for row in para['Filters']:
            filt = ascii.read('../data/apo3_5m/'+Instr_name +"/"+row)
            data_name = row.split('.dat')[0]
            filt_wavelength=filt["col1"]
            filt_throughput=filt["col2"]
            filt_interpolated = interpolate.InterpolatedUnivariateSpline(
            filt_wavelength, filt_throughput)
            setattr(Instrument, data_name, filt_interpolated)
            
            filt_range = [filt["col1"][0], filt["col1"][len(filt["col1"])-1]]
            range_name = row.split('filter.dat')[0]+'range'
            setattr(Instrument, range_name, filt_range)
        
        qefficiency_wavelength=efficiency["col1"]*10 #multiplied by 10 to turn to angstroms
        qefficiency_percent=efficiency["col2"]/100  #divided by 100 to turn into decimal

        
        efficiency_interpolated = interpolate.InterpolatedUnivariateSpline(
                qefficiency_wavelength, qefficiency_percent)
 
        if para['isImager'][0] == 0:
            self.dispersion = para['Dispersion']
        
        
        self.efficiency = efficiency_interpolated
        self.readout_noise = para['readoutnoise[electrons]'][0]
        self.filter_num = para['FilterNum'][0]
        self.gain = para['gain'][0]
        self.scale = para['plate_scale[arcsec/pix]'][0]
# How to read our data files:


