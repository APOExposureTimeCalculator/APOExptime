# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 15:13:34 2019

@author: Alexander
"""

class Sky:
    def __init__(self, lunar_phase = 0, altitude = 2780, seeing = 1):
        self.lunarphase = lunar_phase
        self.altitude = altitude
        self.seeing = seeing
        self.skySED = None #Temporary. This will be SED object that contains SED of sky from file
        
        
class Target:
    def __init__(self, mag, SED = None, Temp = 5778, location, magsystem):
        self.mag = mag
        self.SED = SED
        self.temp = Temp
        
    def starF_lambda(self):
        if self.SED != None:
            #scale SED with magnitude
        else:
            self.temp #use temperature in planck function then scale by magnitude
        
class Counts:
    def __init__(self, airmass = 1, target, instrument, telescope):
        integrate_range = instrument.gprime_range
        filter_profile = instrument.gprime_filter
        detector_qe = instrument.efficiency
        telescope_transm = telescope.transmission
        source = target.F_lambda
        
        interpolationrange = range(integrate_range[0], integrate_range[1])
        h= 6.626*10**(-27) #ergs/s
        c=2.9979*10**(18) #cm/s
        s_integrade = s_integradeInterpolate([filter_profile, detector_qe, telescope_transm, source],  interpolationrange)
        
        s_prime = (1/(h*c))*s_integrade.integral(integrate_range[0], integrate_range[1])
        
        
    def s_integradeInterpolate(self, functions, interpolation_range):
        for i,f in enumerate(functions):
            if i = 0:
                x = np.ones(len(interpolation_range))
            x = f(interpolation_range)*x
            
        return  interpolate.InterpolatedUnivariateSpline(interpolation_range, (x*interpolation_range))
                
        
        
        
        
                 
