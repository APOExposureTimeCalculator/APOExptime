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
    def __init__(self, mag, SED = None, Temp = None)
        