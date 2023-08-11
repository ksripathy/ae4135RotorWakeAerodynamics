from dataclasses import dataclass
from dataclasses_json import dataclass_json
import numpy as np

@dataclass_json
@dataclass
class solution:
    
    #Parameters specified by user
    rotorName: str = "Rotor result object not initialized" #[]
    airfoilName: str = "DU95W180"
    rotorDia: float = 0.0 #[m]
    bladeQuanity: int = 0 #[]
    bladePitch: float = 0.0 #[degrees]
    fstreamVelc: float = 0.0 #[ms^-1]
    fstreamRho: float = 0.0 #[kgm^-3]
    fstreamPres: float = 0.0 #[Pa]
    tsr: float = 0.0 #[]
    bladeRootLoc: float = 0.0 #[]
    bladeTipLoc: float = 0.0 #[]
    bladeTwistDist: str = "Rotor result object not initialized" #[]
    bladeChordDist: str = "Rotor result object not initialized" #[]
    annuliQuantity: int = 0 #[]
    
    #Code Output
    rotorCT : float = 0.0#[]
    rotorCP : float = 0.0 #[]
    annulusLoc : np.float64 = np.zeros(annuliQuantity) #[]
    annulusArea : np.float64 = np.zeros(annuliQuantity) #[m^2]
    annulusBladeAoA: np.float64 = np.zeros(annuliQuantity) #[degrees]
    annulusInflow: np.float64 = np.zeros(annuliQuantity) #[degrees]
    annulusAxialInd : np.float64 = np.zeros(annuliQuantity) #[]
    annulusAzimInd : np.float64 = np.zeros(annuliQuantity) #[]
    annulusCT : np.float64 = np.zeros(annuliQuantity) #[]
    annulusCQ : np.float64 = np.zeros(annuliQuantity) #[]
    annulusCP : np.float64 = np.zeros(annuliQuantity) #[]
    
    def saveAsJson(self, fileDir):
        
        with open(f"{fileDir}/results_{self.rotorName}_vinf={self.fstreamVelc}_tsr={self.tsr}_pitch=" + str(self.bladePitch).replace(".","p") + ".json", 'w') as jsonFile:
            
            jsonFile.write(self.to_json(indent=4))
    
    
    