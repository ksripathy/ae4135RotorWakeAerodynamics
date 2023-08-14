import numpy as np
from miscTools import movingAverage

class mesh:
    
    def __init__(self, rotorObj, windObj, annuliQuantity=50):
        
        self.rotor = rotorObj
        self.wind = windObj
        self.annuliQuantity = annuliQuantity
        
        self.annulusID = np.arange(self.annuliQuantity)
        radialLoc = np.linspace(self.rotor.bladeRootLoc, self.rotor.bladeTipLoc, self.annuliQuantity + 1)
        self.annulusLocRoot = radialLoc[:-1]
        self.annulusLocTip = radialLoc[1:]
        #Annulus location is defined wrt non-dimensionalized radial location as the mid point of radial boundaries
        #self.annulusLoc = movingAverage(radialLoc)
        self.annulusLoc = 0.5 * (self.annulusLocRoot + self.annulusLocTip)
        self.annulusSpan = np.diff(radialLoc) * self.rotor.radius
        #self.annulusArea = 2 * np.pi * (self.annulusLoc * self.radius) * self.annulusSpan
        self.annulusArea = np.pi * ((self.annulusLoc * self.rotor.radius + self.annulusSpan)**2 - (self.annulusLoc * self.rotor.radius)**2) #More accurate annulus area definition
        
        #Initialising annulus blade element parameters
        self.annulusBladeTwistRoot = self.rotor.bladeTwistDist(self.annulusLocRoot) + self.rotor.bladePitch
        self.annulusBladeTwistTip = self.rotor.bladeTwistDist(self.annulusLocTip) + self.rotor.bladePitch
        self.annulusBladeTwist = self.rotor.bladeTwistDist(self.annulusLoc) + self.rotor.bladePitch
        
        self.annulusBladeChordRoot = self.rotor.bladeChordDist(self.annulusLocRoot)
        self.annulusBladeChordTip = self.rotor.bladeChordDist(self.annulusLocTip)
        self.annulusBladeChord = self.rotor.bladeChordDist(self.annulusLoc)
        
        self.rotorOmega = self.rotor.tsr * self.wind.fstreamVelc / self.rotor.radius
        self.annulusTSR = self.rotor.tsr * self.annulusLoc
        
        self.annulusAxialInd = np.zeros(self.annuliQuantity)
        self.annulusAzimInd = np.zeros(self.annuliQuantity)
               
        #Initializing dependent attributes    
        self.annulusAxialVelc = np.zeros(self.annuliQuantity)
        self.annulusAzimVelc = np.zeros(self.annuliQuantity)
        self.annulusResVelc = np.zeros(self.annuliQuantity)
        self.annulusInflow = np.zeros(self.annuliQuantity)
        self.annulusBladeAoA = np.zeros(self.annuliQuantity)
        self.annulusBladeCl = np.zeros(self.annuliQuantity)
        self.annulusBladeCd = np.zeros(self.annuliQuantity)
        self.annulusBladeLift2D = np.zeros(self.annuliQuantity)
        self.annulusBladeDrag2D = np.zeros(self.annuliQuantity)
        self.annulusBladeAxialLoad2D = np.zeros(self.annuliQuantity)
        self.annulusBladeAzimLoad2D = np.zeros(self.annuliQuantity)
        self.annulusBladeAxialLoad3D = np.zeros(self.annuliQuantity)
        self.annulusBladeAzimLoad3D = np.zeros(self.annuliQuantity)
        self.annulusAeroPower = np.zeros(self.annuliQuantity)
        self.annulusMechPower = np.zeros(self.annuliQuantity)
        self.annulusBladeCT = np.zeros(self.annuliQuantity)
        self.annulusBladeCQ = np.zeros(self.annuliQuantity)
        self.annulusBladeAeroCP = np.zeros(self.annuliQuantity)
        self.annulusBladeMechCP = np.zeros(self.annuliQuantity)
        self.annulusAxialIndCorrected = np.zeros(self.annuliQuantity)
        self.annulusAzimIndCorrected = np.zeros(self.annuliQuantity)
        self.annulusMomCT = np.zeros(self.annuliQuantity)
        self.annulusMomCQ = np.zeros(self.annuliQuantity)
        self.annulusMomCP = np.zeros(self.annuliQuantity)
        
        print(f"{self.rotor.name} has been created, discretized into {self.annuliQuantity} sections and initialized!")