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
    def __init__(self, mag, SED = None, Temp = 5778, location):
        self.mag = mag
        self.SED = SED
        self.temp = Temp
        
    def starSED(self):
        if self.SED != None:
            #scale SED with magnitude
        else:
            self.temp #use temperature in planck function then scale by magnitude
        
class Observation:
    def __init__(self, SN, )