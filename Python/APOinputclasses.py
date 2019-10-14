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
        
class Observation:
    def __init__(self, SN, time, airmass = None, target, sky, instrument, telescope):
        lambda1 = instrument.something
        lambda2 = instrument.something
        h= 6.626*10**(-27) #ergs/s
        c=2.9979*10**(10) #cm/s
        s_prime = target.SED.integral(lambda1,lambda2)/(h*c)*(lambda2**2-lambda1**2)/2*sky.skySED.integral(lambda1,lambda2)*
                    instrument.efficiency.integral(lambda1,lambda2)*instrument.filter_wavelengths.integral(lambda1,lambda2)*
                    telescope.transmission.integral(lambda1,lambda2)
        
        
        
                 
