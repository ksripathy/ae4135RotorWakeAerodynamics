import numpy as np
import scipy.optimize as optz
from glauertsCorrection import glauertsCorrectionSingleDOF
from glauertsCorrection import glauertsCorrectionSingleDOFInverse
from prandtlsCorrection import prandtlsCorrectionSingleDOFv2


class BEM:
    
    def __init__(self, meshObj, axialIndInit=0.3, glauertToggle = True, prandtlToggle = True):
        
        self.mesh = meshObj
        
        self.mesh.annulusAxialInd[:] = axialIndInit
        
        self.glauertToggle = glauertToggle
        self.prandtlToggle = prandtlToggle
        
        self.annulusPrandtlCorrectionFactor = np.ones(self.mesh.annuliQuantity)
        
    def indToLoadSingleDOF(self, ID):
        
        #Velocities perceived at blade element
        self.mesh.annulusAxialVelc[ID] = self.mesh.wind.fstreamVelc * (1 - self.mesh.annulusAxialInd[ID])
        self.mesh.annulusAzimVelc[ID] = self.mesh.rotorOmega * (self.mesh.annulusLoc[ID] * self.mesh.rotor.radius) * (1 + self.mesh.annulusAzimInd[ID])
        self.mesh.annulusResVelc[ID] = np.sqrt (self.mesh.annulusAxialVelc[ID]**2 + self.mesh.annulusAzimVelc[ID]**2)
        
        #Blade element angles
        self.mesh.annulusInflow[ID] = np.arctan2(self.mesh.annulusAxialVelc[ID], self.mesh.annulusAzimVelc[ID]) * (180/np.pi)
        self.mesh.annulusBladeAoA[ID] = self.mesh.annulusInflow[ID] - self.mesh.annulusBladeTwist[ID]
         
        #Blade lift and drag
        self.mesh.annulusBladeCl[ID] = np.interp(self.mesh.annulusBladeAoA[ID], self.mesh.rotor.afoilPolar[:,0], self.mesh.rotor.afoilPolar[:,1])
        self.mesh.annulusBladeCd[ID] = np.interp(self.mesh.annulusBladeAoA[ID], self.mesh.rotor.afoilPolar[:,0], self.mesh.rotor.afoilPolar[:,2])
        
        self.mesh.annulusBladeLift2D[ID] = 0.5 * self.mesh.annulusBladeChord[ID] * self.mesh.wind.fstreamRho * (self.mesh.annulusResVelc[ID]**2) * self.mesh.annulusBladeCl[ID]
        self.mesh.annulusBladeDrag2D[ID] = 0.5 * self.mesh.annulusBladeChord[ID] * self.mesh.wind.fstreamRho * (self.mesh.annulusResVelc[ID]**2) * self.mesh.annulusBladeCd[ID]
        
        #Annulus normal and tangential laoding
        self.mesh.annulusBladeAxialLoad2D[ID] = (self.mesh.annulusBladeLift2D[ID] * np.cos(self.mesh.annulusInflow[ID] * (np.pi/180)) + self.mesh.annulusBladeDrag2D[ID] * np.sin(self.mesh.annulusInflow[ID] * (np.pi/180))) * self.mesh.rotor.bladeQuantity
        self.mesh.annulusBladeAzimLoad2D[ID] = (self.mesh.annulusBladeLift2D[ID] * np.sin(self.mesh.annulusInflow[ID] * (np.pi/180)) - self.mesh.annulusBladeDrag2D[ID] * np.cos(self.mesh.annulusInflow[ID] * (np.pi/180))) * self.mesh.rotor.bladeQuantity
        
        self.mesh.annulusBladeAxialLoad3D[ID] = self.mesh.annulusBladeAxialLoad2D[ID] * self.mesh.annulusSpan[ID]
        self.mesh.annulusAeroPower[ID] = self.mesh.annulusBladeAxialLoad3D[ID] * self.mesh.annulusAxialVelc[ID]
        
        self.mesh.annulusBladeAzimLoad3D[ID] = self.mesh.annulusBladeAzimLoad2D[ID] * self.mesh.annulusSpan[ID]
        self.mesh.annulusMechPower[ID] = self.mesh.annulusBladeAzimLoad3D[ID] * self.mesh.annulusLoc[ID] * self.mesh.rotor.radius * self.mesh.rotorOmega
        
        #Annulus thrust coefficieint
        self.mesh.annulusBladeCT[ID] = self.mesh.annulusBladeAxialLoad3D[ID]/(0.5 * self.mesh.wind.fstreamRho * self.mesh.wind.fstreamVelc**2 * self.mesh.annulusArea[ID])
        
        #Annulus torque coefficient
        self.mesh.annulusBladeCQ[ID] = (self.mesh.annulusLoc[ID] * self.mesh.annulusBladeAzimLoad3D[ID])/(0.5 * self.mesh.wind.fstreamRho * self.mesh.wind.fstreamVelc**2 * self.mesh.annulusArea[ID])
        
        #Annulus power coefficient
        self.mesh.annulusBladeAeroCP[ID] = self.mesh.annulusAeroPower[ID]/(0.5 * self.mesh.wind.fstreamRho * self.mesh.wind.fstreamVelc**3 * self.mesh.annulusArea[ID])
        self.mesh.annulusBladeMechCP[ID] = self.mesh.annulusMechPower[ID]/(0.5 * self.mesh.wind.fstreamRho * self.mesh.wind.fstreamVelc**3 * self.mesh.annulusArea[ID])
        
    def loadToIndSingleDOF(self, ID):
        
        if self.glauertToggle:
            
            self.mesh.annulusAxialIndCorrected[ID] = glauertsCorrectionSingleDOF(self.mesh.annulusBladeCT[ID])
            
        else:
            self.mesh.annulusAxialIndCorrected[ID] = 0.5 * (1 - np.sqrt(1 - self.mesh.annulusBladeCT[ID]))
            
        self.mesh.annulusMomCT[ID] = glauertsCorrectionSingleDOFInverse(self.mesh.annulusAxialIndCorrected[ID])
        
        if self.prandtlToggle:
            
            self.annulusPrandtlCorrectionFactor[ID] = prandtlsCorrectionSingleDOFv2(self.mesh.rotor.bladeQuantity, self.mesh.rotor.bladeRootLoc, self.mesh.annulusLoc[ID], self.mesh.annulusTSR[ID], self.mesh.annulusAxialIndCorrected[ID])
            self.mesh.annulusAxialIndCorrected[ID] = self.mesh.annulusAxialIndCorrected[ID] / self.annulusPrandtlCorrectionFactor[ID]
        
        self.mesh.annulusAzimIndCorrected[ID] = self.mesh.annulusBladeAzimLoad2D[ID] / (2 * self.mesh.wind.fstreamRho * (2 * np.pi* self.mesh.annulusLoc[ID] * self.mesh.rotor.radius) * self.mesh.wind.fstreamVelc**2 * (1 - self.mesh.annulusAxialIndCorrected[ID]) * self.mesh.rotor.tsr * self.mesh.annulusLoc[ID])
        
        if self.prandtlToggle:
            
            self.mesh.annulusAzimIndCorrected[ID] = self.mesh.annulusAzimIndCorrected[ID] / self.annulusPrandtlCorrectionFactor[ID]
        
        #As per Sorenson when Tip correction is on, either the blade element load coefficieints needed to be divided by prandtl correction factor or the momentum theory coefficients need to be multiplied for both of them to match
        self.mesh.annulusMomCT[ID] = glauertsCorrectionSingleDOFInverse(self.mesh.annulusAxialIndCorrected[ID] * self.annulusPrandtlCorrectionFactor[ID])
        self.mesh.annulusMomCQ[ID] = 4 * (self.mesh.annulusAzimIndCorrected[ID] * self.annulusPrandtlCorrectionFactor[ID]) * (1 - self.mesh.annulusAxialIndCorrected[ID]) * self.mesh.annulusTSR[ID] * self.mesh.annulusLoc[ID]
        self.mesh.annulusMomCP[ID] =  self.mesh.annulusMomCT[ID] * (1 - self.mesh.annulusAxialIndCorrected[ID])
        
    def classicSolver(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3, printIter = False, solverLog = False):
        
        for annulusID in range(self.mesh.annuliQuantity):
            
            print(f"Solving for annulus ID = {annulusID} ...")
            
            for i in range(maxIterations):
                
                self.indToLoadSingleDOF(annulusID)
                self.loadToIndSingleDOF(annulusID)
                
                iterationError = np.absolute(self.mesh.annulusAxialInd[annulusID] - self.mesh.annulusAxialIndCorrected[annulusID])
                
                if printIter:
                    
                    print("")
                    print(f"Iteration-{i}")
                    print(f"Axial Induction = {self.mesh.annulusAxialIndCorrected[annulusID]}")
                    print(f"Azimuthal Induction = {self.mesh.annulusAzimIndCorrected[annulusID]}")
                
                if np.absolute(iterationError) < minError or i == maxIterations - 1:
                    break
                
                self.mesh.annulusAxialInd[annulusID] = relaxationFactor * self.mesh.annulusAxialIndCorrected[annulusID] + (1 - relaxationFactor) * self.mesh.annulusAxialInd[annulusID]
                self.mesh.annulusAzimInd[annulusID] = relaxationFactor * self.mesh.annulusAzimIndCorrected[annulusID] + (1 - relaxationFactor) * self.mesh.annulusAzimInd[annulusID]
                
            #Compute final induction errors
            axialIndError = np.absolute(self.mesh.annulusAxialInd[annulusID] - self.mesh.annulusAxialIndCorrected[annulusID])
            azimIndError = np.absolute(self.mesh.annulusAzimInd[annulusID] - self.mesh.annulusAzimIndCorrected[annulusID])
            
            print(f"Annulus {annulusID} solved after {i + 1} iterations!")
            
            if solverLog:
                
                self.solverLog(annulusID, axialIndError, azimIndError)
        
    #Optimization function for the optimization solver. Have to use static method since the scipy optimization requires the first argument of objective function to have the objective variables and with a class method the first argument needs to be referenced to the object
    @staticmethod
    def optzFunction(optzArg, BEMObject, annulusID):
    
        #Since optimizer chooses the input argument randomly, azimuthal induction needs to be adjusted accordingly apriori in order to represent the inflow accurately
        BEMObject.mesh.annulusAxialInd[annulusID] = optzArg
        BEMObject.mesh.annulusAzimInd[annulusID] = BEMObject.mesh.annulusBladeAzimLoad2D[annulusID] / (2 * BEMObject.mesh.fstreamRho * (2 * np.pi* BEMObject.mesh.rotor.annulusLoc[annulusID] * BEMObject.mesh.rotor.rotorRadius) * BEMObject.mesh.fstreamVelc**2 * (1 - BEMObject.mesh.annulusAxialIndCorrected[annulusID]) * BEMObject.mesh.tsr * BEMObject.mesh.rotor.annulusLoc[annulusID])
        #Prandtl correction if any
        BEMObject.mesh.annulusAzimInd[annulusID] = BEMObject.mesh.annulusAzimInd[annulusID]/BEMObject.annulusPrandtlCorrectionFactor[annulusID]
        
        BEMObject.indToLoadSingleDOF(annulusID)
        BEMObject.loadToIndSingleDOF(annulusID)
        
        #axialIndError = BEMObject.mesh.annulusAxialInd[annulusID] - BEMObject.mesh.annulusAxialIndCorrected[annulusID]
        axialIndError = optzArg - BEMObject.mesh.annulusAxialIndCorrected[annulusID]
        
        iterationError = np.absolute(axialIndError)
        
        return iterationError
            
    def optzSolver(self, solverLog = False):
        
        for annulusID in range(self.mesh.annuliQuantity):
        
            print(f"Solving for annulus ID = {annulusID} ...")
            
            res = optz.direct(self.optzFunction, bounds = [(0,0.95),], args = (self, annulusID), maxfun=10000, maxiter=10000, len_tol=1e-10)
            print(res.message)
            
            #Compute final induction errors
            axialIndError = np.absolute(self.mesh.annulusAxialInd[annulusID] - self.mesh.annulusAxialIndCorrected[annulusID])
            azimIndError = np.absolute(self.mesh.annulusAzimInd[annulusID] - self.mesh.annulusAzimIndCorrected[annulusID])
            
            if solverLog:
                
                self.solverLog(annulusID, axialIndError, azimIndError)
        
    def solverLog(self, annulusID, axialIndError, azimIndError):
        
        print("")
        print(f"Prandtl correction factor = {self.annulusPrandtlCorrectionFactor[annulusID]}")
        print("")
        print(f"Annulus CT from blade element theory = {self.mesh.annulusBladeCT[annulusID]}")
        print(f"Annulus CQ from blade element theory = {self.mesh.annulusBladeCQ[annulusID]}")
        print(f"Annulus aerodynamic CP from blade element theory = {self.mesh.annulusBladeAeroCP[annulusID]}")
        print(f"Annulus mechanical CP from blade element theory = {self.mesh.annulusBladeMechCP[annulusID]}")
        print("")
        print(f"Annulus CT from momentum theory = {self.mesh.annulusMomCT[annulusID]}")
        print(f"Annulus CQ from momentum theory = {self.mesh.annulusMomCQ[annulusID]}")
        print(f"Annulus CP from momentum theory = {self.mesh.annulusMomCP[annulusID]}")
        print("")
        print(f"Final annulus axial induction = {self.mesh.annulusAxialIndCorrected[annulusID]}")
        print(f"Final annulus azimuthal induction = {self.mesh.annulusAzimIndCorrected[annulusID]}")
        print("")
        print(f"Final axial induction error = {axialIndError}")
        print(f"Final azimuthal induction error = {azimIndError}")
        print("===================================================================")
        
        
    
        
            
            
                    
            
            
        
        
        
           
        
        
    
        
        
        
        
        
        
        