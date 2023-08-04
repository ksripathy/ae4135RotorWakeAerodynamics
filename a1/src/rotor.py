import numpy as np
from miscTools import movingAverage

class rotor:
    
    #Class Methods
    
    #Initiation
    def __init__(self, rotorDia, bladeQuantity, bladeRootLoc, bladeTipLoc, bladeTwistDist, bladeChordDist, afoilPolar, annuliQuantity=50, name = "Rotor-default"):
        
        self.rotorDia = rotorDia
        self.bladeQuantity = bladeQuantity
        self.bladeRootLoc = bladeRootLoc
        self.bladeTipLoc = bladeTipLoc
        self.bladeTwistDist = bladeTwistDist
        self.bladeChordDist = bladeChordDist
        self.afoilPolar = afoilPolar
        self.annuliQuantity = annuliQuantity
        self.name = name
        
        self.rotorRadius = 0.5 * rotorDia
        self.rotorArea = np.pi * self.rotorRadius**2
        
        self.annulusID = np.arange(self.annuliQuantity)
        radialLoc = np.linspace(self.bladeRootLoc, self.bladeTipLoc, self.annuliQuantity + 1)
        #Annulus location is defined wrt non-dimensionalized radial location as the mid point of radial boundaries
        self.annulusLoc = movingAverage(radialLoc)
        self.annulusSpan = np.diff(radialLoc) * self.rotorRadius
        self.annulusArea = 2 * np.pi * (self.annulusLoc * self.rotorRadius) * self.annulusSpan
        
        #Initialising annulus blade element parameters
        self.annulusBladeTwist = self.bladeTwistDist(self.annulusLoc)
        self.annulusBladeChord = self.bladeChordDist(self.annulusLoc)
        
        print(f"{name} has been created and discretized into {self.annuliQuantity} sections. Annuli blade element parameters initialized!")
        
        
        
        
        
        
        
        
        
    
        
        
        
        
        
        
        
        