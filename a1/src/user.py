import numpy as np
from miscTools import movingAverage

class rotor:
    
    #Class Methods
    
    #Initiation
    def __init__(self, dia, tsr, bladeQuantity, bladeRootLoc, bladeTipLoc, bladePitch, bladeTwistDist, bladeChordDist, afoilPolar , name = "Rotor-default"):
        
        self.dia = dia
        self.tsr = tsr
        self.bladeQuantity = bladeQuantity
        self.bladeRootLoc = bladeRootLoc
        self.bladeTipLoc = bladeTipLoc
        self.bladePitch = bladePitch
        self.bladeTwistDist = bladeTwistDist
        self.bladeChordDist = bladeChordDist
        self.afoilPolar = afoilPolar
        self.name = name
        
        self.radius = 0.5 * dia
        self.area = np.pi * self.radius**2
        
class wind:
        
    def __init__(self, fstreamVelc, fstreamRho, fstreamPres):
        
        self.fstreamVelc = fstreamVelc
        self.fstreamRho = fstreamRho
        self.fstreamPres = fstreamPres
        
        
        
        
        
        
        
        
        
    
        
        
        
        
        
        
        
        