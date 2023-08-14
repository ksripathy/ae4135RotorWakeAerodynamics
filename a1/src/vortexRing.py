import numpy as np
class vortexRing:
    
    def __init__(self, rootLoc, tipLoc, rootChord, tipChord, rootTwist, tipTwist, rotorRadius, bladeAzimInit):
        
        #Translating degrees to radians
        bladeAzimInit = bladeAzimInit * np.pi/180
        rootTwist = rootTwist * np.pi/180
        tipTwist = tipTwist * np.pi/180
        
        #Initiating vertices of vortex ring with resepct to global reference frame
        self.leadRootLoc = np.array([0.0, rootLoc * rotorRadius * np.cos(bladeAzimInit), rootLoc * rotorRadius * np.sin(bladeAzimInit)])
        self.leadTipLoc = np.array([0.0, tipLoc * rotorRadius * np.cos(bladeAzimInit), tipLoc * rotorRadius * np.sin(bladeAzimInit)])
        self.controlPointLoc = 0.5 * (self.leadRootLoc + self.leadTipLoc)
        self.trailRootLoc = self.leadRootLoc + np.array([rootChord * np.abs(np.sin(rootTwist)), rootChord * np.cos(rootTwist) * np.cos(bladeAzimInit), rootChord * np.cos(rootTwist) * np.sin(bladeAzimInit)])
        self.trailTipLoc = self.leadTipLoc + np.array([tipChord * np.abs(np.sin(tipTwist)), tipChord * np.cos(tipTwist) * np.cos(bladeAzimInit), tipChord * np.cos(tipTwist) * np.sin(bladeAzimInit)])
        
        