# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:13:34 2019

@author: Alexander, Bryson, Roger, Hassan, Manny
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
        x = sp_new(range(1000,20000) * u.AA, flux_unit= units.FLAM)
        self.F_lambda = interpolate.InterpolatedUnivariateSpline(range(1000,20000), x)
        
class Observation:

   
        def __init__(self,  target, telescope=None, airmass = 1):
            
           # telescope_transm = telescope.transmission
            self. telescope_area = (175**2)*3.14
            self.source = target.F_lambda
            
        def counts(self, instrument):
            att = dir(instrument)
            self.detector_qe = instrument.efficiency
            for row in att:
                if row.find('filter') > 0:
                    
                    filter_profile = getattr(instrument, row)
                    integrate_range = getattr(instrument,row.replace('filter', 'range'))
                    interpolationrange = range(integrate_range[0], integrate_range[1])
                    h= 6.626*10**(-27) #ergs*s
                    c=2.9979*10**(18) #A/s
                    s_integrade = s_integradeInterpolate([self.source, self.detector_qe, filter_profile], interpolationrange)
                    
                    s_prime = self.telescope_area*(1/(h*c))*s_integrade.integral(integrate_range[0], integrate_range[1])
                    count_name = row.replace('_filter', '') +'_countrate'
                    setattr(Observation, count_name, s_prime)
                    
            
#        def SNfromTime(self, exptime):
#            
#            
#        def TimefromSN(self, SN):
            
            
            
            
def s_integradeInterpolate(functions, interpolation_range):
    for i,f in enumerate(functions):
        if i == 0:
            x = np.ones(len(interpolation_range))
        x = f(interpolation_range)*x
        
    return  interpolate.InterpolatedUnivariateSpline(interpolation_range, (x*interpolation_range))
            
    
        
        
        
                 
