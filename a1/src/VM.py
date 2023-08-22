from vortexRing import vortexRing
import numpy as np
import copy

class VM:
    
    def __init__(self, meshObj, wakeCoverage, wakeVelc, deltaTime):
        
        self.mesh = meshObj
        self.wakeCoverage = wakeCoverage
        self.wakeVelc = wakeVelc
        self.deltaTime = deltaTime
        totalTime = self.wakeCoverage/self.wakeVelc[0]
        self.totalTimeSteps = int(totalTime/deltaTime) 
        
        for i in range(self.mesh.rotor.bladeQuantity):
            
             bladeObj = getattr(self.mesh, f"annulusBlade{i}")
             setattr(bladeObj, "rootVortHistory", np.empty(self.totalTimeSteps,dtype=np.ndarray))
             setattr(bladeObj, "tipVortHistory", np.empty(self.totalTimeSteps,dtype=np.ndarray))
             setattr(bladeObj, "QChordHistory", np.empty(self.totalTimeSteps,dtype=np.ndarray))
             
             #Initalising each time step of Qchord history with all annuli
             for j in range(self.totalTimeSteps):
                 
                 bladeObj.QChordHistory[j] = np.empty(self.mesh.annuliQuantity, dtype=np.ndarray)
                 
                 
             
        #self.rootVortHistory = np.empty(self.totalTimeSteps,dtype=np.ndarray)
        #self.tipVortHistory = np.empty(self.totalTimeSteps,dtype=np.ndarray)
        
        
    def buildVortRings(self):
        
        #Iterating over blades
        for i in range(self.mesh.rotor.bladeQuantity):
            
            #Getattr needs to be used iterate over different runtime generated blade objects
            bladeObj = getattr(self.mesh, f"annulusBlade{i}")
            
            #Iterating over each annuli of blade
            for j in range(self.mesh.annuliQuantity):
        
                #For each annulus generate vortexSystem object
                bladeObj.vortSys[j] = np.empty(self.totalTimeSteps+1, dtype=np.ndarray)
            
                #Initialize bound vortex ring
                bladeObj.vortSys[j][0] = vortexRing(self.mesh.annulusLocRoot[j], self.mesh.annulusLocTip[j], self.mesh.annulusBladeChordRoot[j], self.mesh.annulusBladeChordTip[j], self.mesh.annulusBladeTwistRoot[j], self.mesh.annulusBladeTwistTip[j], self.mesh.rotor.radius, bladeObj.azimInit)
                
                #Initialize trailing vortex rings
                for k in range(1,self.totalTimeSteps+1):
                    
                    bladeObj.vortSys[j][k] = vortexRing()
                    
                #Time marching to propagate bound vortex rings to wake
                for k in range(self.totalTimeSteps):
                    
                    #Convect existing vortex rings. Values are updated beginning with rings located in the fathermost wake region
                    for n in range(k,-1,-1):
                        
                        bladeObj.vortSys[j][n].convect(bladeObj.vortSys[j][n+1],self.wakeVelc,self.deltaTime)                        
                        
                    #Rotate the rotor by one time step thereby the locations of bound vortex as well
                    bladeObj.vortSys[j][0].rotate(self.mesh.rotorOmega, self.deltaTime)
                    
                    #Recording root and tip trailing vortex history
                    if j == 0:
                        
                        bladeObj.rootVortHistory[k] = copy.deepcopy(bladeObj.vortSys[j][:k+2])
                        
                    if j == self.mesh.annuliQuantity - 1:
                        
                        bladeObj.tipVortHistory[k] = copy.deepcopy(bladeObj.vortSys[j][:k+2])
                        
                    #Recording bound vortex history
                    bladeObj.QChordHistory[k][j] = copy.deepcopy(bladeObj.vortSys[j][0])
            
            
                        
                        
                    
                    
                    
                
                
                
        
        
        
