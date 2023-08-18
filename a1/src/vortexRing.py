import numpy as np

class vortexRing:
    
    def __init__(self, rootLoc=0, tipLoc=0, rootChord=0, tipChord=0, rootTwist=0, tipTwist=0, rotorRadius=0, bladeAzimInit=0):
        
        #Translating degrees to radians
        bladeAzimInit = bladeAzimInit * np.pi/180
        rootTwist = rootTwist * np.pi/180
        tipTwist = tipTwist * np.pi/180
        
        #Define local blade geometry
        self.rootChord = rootChord
        self.tipChord = tipChord
        self.rootTwist = rootTwist
        self.tipTwist = tipTwist
        
        #Define radius
        self.rootRad = rootLoc * rotorRadius
        self.tipRad = tipLoc * rotorRadius
        
        #Define blade azimuthal position
        self.bladeAzimPos = bladeAzimInit
        
        #Define convection displacement. For bound vortex ring it is 0
        self.convDisp = np.zeros(3)
        
        #Initiating vertices of vortex ring with resepct to global reference frame
        self.computeLoc()
        
        
    #Compute all coordinates of vortex ring
    def computeLoc(self):
        
        self.leadRootLoc = np.array([0, self.rootRad * np.sin(self.bladeAzimPos), self.rootRad * np.cos(self.bladeAzimPos)]) + self.convDisp
        self.leadTipLoc = np.array([0, self.tipRad * np.sin(self.bladeAzimPos), self.tipRad * np.cos(self.bladeAzimPos)]) + self.convDisp
        self.controlPointLoc = 0.5 * (self.leadRootLoc + self.leadTipLoc) + self.convDisp
        self.trailRootLoc = self.leadRootLoc + np.array([self.rootChord * np.abs(np.sin(self.rootTwist)), self.rootChord * np.cos(self.rootTwist) * np.sin(self.bladeAzimPos), self.rootChord * np.cos(self.rootTwist) * np.cos(self.bladeAzimPos)]) + self.convDisp
        self.trailTipLoc = self.leadTipLoc + np.array([self.tipChord * np.abs(np.sin(self.tipTwist)), self.tipChord * np.cos(self.tipTwist) * np.sin(self.bladeAzimPos), self.tipChord * np.cos(self.tipTwist) * np.cos(self.bladeAzimPos)]) + self.convDisp
    
    #Rotate bound vortex ring due to rotor rotation 
    def rotate(self, rotorOmega, deltaTime):
        
        self.bladeAzimPos = self.bladeAzimPos + rotorOmega * deltaTime
        self.computeLoc()
        
    #Convect bound vortex ring as trailing vortex rings
    def convect(self, destRing, wakeVelc, deltaTime):
        
        destRing.convDisp = wakeVelc * deltaTime
        
        destRing.leadRootLoc = self.leadRootLoc + destRing.convDisp
        destRing.leadTipLoc = self.leadTipLoc + destRing.convDisp
        destRing.controlPointLoc = self.controlPointLoc + destRing.convDisp
        destRing.trailRootLoc = self.trailRootLoc + destRing.convDisp
        destRing.trailTipLoc = self.trailTipLoc + destRing.convDisp
        
    def plotRings(self, pltAxisObj):
        
        line1 = np.array([self.leadRootLoc, self.leadTipLoc])
        line2 = np.array([self.leadTipLoc, self.trailTipLoc])
        line3 = np.array([self.trailTipLoc, self.trailRootLoc])
        line4 = np.array([self.trailRootLoc, self.leadRootLoc])
        
        pltAxisObj.plot(line1[:,0],line1[:,1],line1[:,2], 'r')
        pltAxisObj.plot(line2[:,0],line2[:,1],line2[:,2], 'r')
        pltAxisObj.plot(line3[:,0],line3[:,1],line3[:,2], 'r')
        pltAxisObj.plot(line4[:,0],line4[:,1],line4[:,2], 'r')
        
    def plotControlPoint(self, pltAxisObj, markerColor):
        
        pltAxisObj.scatter(self.controlPointLoc[0], self.controlPointLoc[1], self.controlPointLoc[2], marker='.', color=markerColor)
        