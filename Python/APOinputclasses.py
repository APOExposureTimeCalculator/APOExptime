# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:13:34 2019

@author: Alexander
"""
from synphot.models import BlackBody1D
from synphot import units
from astropy import units as u
import numpy as np
from synphot import SourceSpectrum
from scipy import interpolate

class Sky:
    def __init__(self, lunar_phase = 0, seeing = 1):
        self.lunarphase = lunar_phase
        self.seeing = seeing
        self.skySED = None #Temporary. This will be SED object that contains SED of sky from file
        
        
class Target:
    def __init__(self, mag, magsystem, filtRange, SED = None, Temp = 5778, location=None):
        if magsystem == 'VEGAMAG':
            sys = units.VEGAMAG
        elif magsystem == 'stmag':
            sys = u.STmag
        elif magsystem == 'abnu':
            sys = u.ABmag
        
        vega = SourceSpectrum.from_vega()
        
        self.mag = mag
        self.SED = SED
        self.temp = Temp
        self.inputFlux = units.convert_flux(filtRange, mag*sys, units.FLAM, vegaspec = vega)
        self.range = filtRange
        
        
    def starF_lambda(self):
        sp = SourceSpectrum(BlackBody1D, temperature=self.temp*u.K)
        sp_new = sp/np.mean(sp(self.range * u.AA, flux_unit= units.FLAM)/self.inputFlux)
        x = sp_new(range(1000,9000) * u.AA, flux_unit= units.FLAM)
        self.F_lambda = interpolate.InterpolatedUnivariateSpline(range(1000,9000), x)
        
class Counts:
        """Instantiates an object of the Instrument class.
    
    Args:
        Instr_name(str): Name of Instrument to be used.
        
    Attributes:
        efficiency(object): UnivariateInterpolatedSpline of instrument efficiency.
        readout_noise(float): Value of instrument readout noise
        filter_num(int): Number of filters for instrument
        gain(float): Gain of instrument
        """
        def __init__(self,  target, instrument, telescope=None, airmass = 1):
            self.integrate_range = instrument.gprime_range
            self.filter_profile = instrument.gprime_filter
            self.detector_qe = instrument.efficiency
           # telescope_transm = telescope.transmission
            self. telescope_area = (175**2)*3.14
            self.source = target.F_lambda
            
            self. interpolationrange = range(self.integrate_range[0], self.integrate_range[1])
            h= 6.626*10**(-27) #ergs*s
            c=2.9979*10**(18) #A/s
            s_integrade = s_integradeInterpolate([self.source, self.detector_qe, self.filter_profile], self.interpolationrange)
            
            s_prime = self.telescope_area*(1/(h*c))*s_integrade.integral(self.integrate_range[0], self.integrate_range[1])
            
            self.s_prime = s_prime
            
    
def s_integradeInterpolate(functions, interpolation_range):
    for i,f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range)*x
        
    return  interpolate.InterpolatedUnivariateSpline(interpolation_range, (x*interpolation_range))
            
    
        
        
        
                 
