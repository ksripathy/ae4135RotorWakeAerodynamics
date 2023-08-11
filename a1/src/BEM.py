import numpy as np
import scipy.optimize as optz
from glauertsCorrection import glauertsCorrectionSingleDOF
from glauertsCorrection import glauertsCorrectionSingleDOFInverse
from prandtlsCorrection import prandtlsCorrectionSingleDOFv2
from dataGen import rotorResults
import inspect


class BEM:
    
    def __init__(self, rotorObject, bladePitch, fstreamVelc, fstreamRho, fstreamPres, tsr, axialIndInit = 0.3, glauertToggle = True, prandtlToggle = True):
        
        self.rotor = rotorObject
        self.bladePitch = bladePitch
        self.fstreamVelc = fstreamVelc
        self.fstreamRho = fstreamRho
        self.fstreamPres = fstreamPres
        self.tsr = tsr
        self.axialIndInit = axialIndInit
        
        self.rotorOmega = tsr * fstreamVelc / self.rotor.rotorRadius
        self.annulusTSR = tsr * self.rotor.annulusLoc
        
        self.annulusAxialInd = axialIndInit * np.ones(self.rotor.annuliQuantity)
        self.annulusAzimInd = np.zeros(self.rotor.annuliQuantity)
        
        self.glauertToggle = glauertToggle
        self.prandtlToggle = prandtlToggle
               
        #Initializing dependent attributes    
        self.annulusAxialVelc = np.zeros(self.rotor.annuliQuantity)
        self.annulusAzimVelc = np.zeros(self.rotor.annuliQuantity)
        self.annulusResVelc = np.zeros(self.rotor.annuliQuantity)
        self.annulusInflow = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAoA = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeCl = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeCd = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeLift2D = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeDrag2D = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAxialLoad2D = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAzimLoad2D = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAxialLoad3D = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAzimLoad3D = np.zeros(self.rotor.annuliQuantity)
        self.annulusAeroPower = np.zeros(self.rotor.annuliQuantity)
        self.annulusMechPower = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeCT = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeCQ = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeAeroCP = np.zeros(self.rotor.annuliQuantity)
        self.annulusBladeMechCP = np.zeros(self.rotor.annuliQuantity)
        self.annulusAxialIndCorrected = np.zeros(self.rotor.annuliQuantity)
        self.annulusAzimIndCorrected = np.zeros(self.rotor.annuliQuantity)
        self.annulusMomCT = np.zeros(self.rotor.annuliQuantity)
        self.annulusMomCQ = np.zeros(self.rotor.annuliQuantity)
        self.annulusMomCP = np.zeros(self.rotor.annuliQuantity)
        self.annulusPrandtlCorrectionFactor = np.ones(self.rotor.annuliQuantity)
        self.rotorCT = 0
        self.rotorAeroCP = 0
        self.rotorMechCP = 0
        
        print(f"Momentum parameters initialized for {self.rotor.name}!")
        
    def indToBladeLoadSingleDOF(self, ID):
        
        #Velocities perceived at blade element
        self.annulusAxialVelc[ID] = self.fstreamVelc * (1 - self.annulusAxialInd[ID])
        self.annulusAzimVelc[ID] = self.rotorOmega * (self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * (1 + self.annulusAzimInd[ID])
        self.annulusResVelc[ID] = np.sqrt (self.annulusAxialVelc[ID]**2 + self.annulusAzimVelc[ID]**2)
        
        #Blade element angles
        self.annulusInflow[ID] = np.arctan2(self.annulusAxialVelc[ID], self.annulusAzimVelc[ID]) * (180/np.pi)
        self.annulusBladeAoA[ID] = self.annulusInflow[ID] - (self.rotor.annulusBladeTwist[ID] + self.bladePitch)
         
        #Blade lift and drag
        self.annulusBladeCl[ID] = np.interp(self.annulusBladeAoA[ID], self.rotor.afoilPolar[:,0], self.rotor.afoilPolar[:,1])
        self.annulusBladeCd[ID] = np.interp(self.annulusBladeAoA[ID], self.rotor.afoilPolar[:,0], self.rotor.afoilPolar[:,2])
        
        self.annulusBladeLift2D[ID] = 0.5 * self.rotor.annulusBladeChord[ID] * self.fstreamRho * (self.annulusResVelc[ID]**2) * self.annulusBladeCl[ID]
        self.annulusBladeDrag2D[ID] = 0.5 * self.rotor.annulusBladeChord[ID] * self.fstreamRho * (self.annulusResVelc[ID]**2) * self.annulusBladeCd[ID]
        
        #Annulus normal and tangential laoding
        self.annulusBladeAxialLoad2D[ID] = (self.annulusBladeLift2D[ID] * np.cos(self.annulusInflow[ID] * (np.pi/180)) + self.annulusBladeDrag2D[ID] * np.sin(self.annulusInflow[ID] * (np.pi/180))) * self.rotor.bladeQuantity
        self.annulusBladeAzimLoad2D[ID] = (self.annulusBladeLift2D[ID] * np.sin(self.annulusInflow[ID] * (np.pi/180)) - self.annulusBladeDrag2D[ID] * np.cos(self.annulusInflow[ID] * (np.pi/180))) * self.rotor.bladeQuantity
        
        self.annulusBladeAxialLoad3D[ID] = self.annulusBladeAxialLoad2D[ID] * self.rotor.annulusSpan[ID]
        self.annulusAeroPower[ID] = self.annulusBladeAxialLoad3D[ID] * self.annulusAxialVelc[ID]
        
        self.annulusBladeAzimLoad3D[ID] = self.annulusBladeAzimLoad2D[ID] * self.rotor.annulusSpan[ID]
        self.annulusMechPower[ID] = self.annulusBladeAzimLoad3D[ID] * self.rotor.annulusLoc[ID] * self.rotor.rotorRadius * self.rotorOmega
        
        #Annulus thrust coefficieint
        self.annulusBladeCT[ID] = self.annulusBladeAxialLoad3D[ID]/(0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.annulusArea[ID])
        
        #Annulus torque coefficient
        self.annulusBladeCQ[ID] = (self.rotor.annulusLoc[ID] * self.annulusMomAzimLoad3D[ID])/(0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.annulusArea[ID])
        
        #Annulus power coefficient
        self.annulusBladeAeroCP[ID] = self.annulusAeroPower[ID]/(0.5 * self.fstreamRho * self.fstreamVelc**3 * self.rotor.annulusArea[ID])
        self.annulusBladeMechCP[ID] = self.annulusMechPower[ID]/(0.5 * self.fstreamRho * self.fstreamVelc**3 * self.rotor.annulusArea[ID])
        
    def momLoadToIndSingleDOF(self, ID):
        
        if self.glauertToggle:
            
            self.annulusAxialIndCorrected[ID] = glauertsCorrectionSingleDOF(self.annulusBladeCT[ID])
            
        else:
            self.annulusAxialIndCorrected[ID] = 0.5 * (1 - np.sqrt(1 - self.annulusBladeCT[ID]))
            
        self.annulusMomCT[ID] = glauertsCorrectionSingleDOFInverse(self.annulusAxialIndCorrected[ID])
        
        if self.prandtlToggle:
            
            self.annulusPrandtlCorrectionFactor[ID] = prandtlsCorrectionSingleDOFv2(self.rotor.bladeQuantity, self.rotor.bladeRootLoc, self.rotor.annulusLoc[ID], self.annulusTSR[ID], self.annulusAxialIndCorrected[ID])
            self.annulusAxialIndCorrected[ID] = self.annulusAxialIndCorrected[ID] / self.annulusPrandtlCorrectionFactor[ID]
        
        self.annulusAzimIndCorrected[ID] = self.annulusBladeAzimLoad2D[ID] / (2 * self.fstreamRho * (2 * np.pi* self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * self.fstreamVelc**2 * (1 - self.annulusAxialIndCorrected[ID]) * self.tsr * self.rotor.annulusLoc[ID])
        
        if self.prandtlToggle:
            
            self.annulusAzimIndCorrected[ID] = self.annulusAzimIndCorrected[ID] / self.annulusPrandtlCorrectionFactor[ID]
        
        #For some reason induction input for CT from momentum theory needs to have the prandtl correction factor scaled off for it to match CT from blade element theory
        self.annulusMomCT[ID] = glauertsCorrectionSingleDOFInverse(self.annulusAxialIndCorrected[ID] * self.annulusPrandtlCorrectionFactor[ID])
        self.annulusMomCQ[ID] = 4 * self.annulusAzimIndCorrected[ID] * (1 - self.annulusAxialIndCorrected[ID]) * self.annulusTSR[ID] * self.rotor.annulusLoc[ID]
        self.annulusMomCP[ID] =  self.annulusMomCT[ID] * (1 - self.annulusAxialIndCorrected[ID])
        
    def classicSolver(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3, printIter = False, solverLog = False):
        
        for annulusID in range(self.rotor.annuliQuantity):
            
            print(f"Solving for annulus ID = {annulusID} ...")
            
            for i in range(maxIterations):
                
                self.indToBladeLoadSingleDOF(annulusID)
                self.momLoadToIndSingleDOF(annulusID)
                
                iterationError = np.absolute(self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID])
                
                if printIter:
                    
                    print("")
                    print(f"Iteration-{i}")
                    print(f"Axial Induction = {self.annulusAxialIndCorrected[annulusID]}")
                    print(f"Azimuthal Induction = {self.annulusAzimIndCorrected[annulusID]}")
                
                if np.absolute(iterationError) < minError or i == maxIterations - 1:
                    break
                
                self.annulusAxialInd[annulusID] = relaxationFactor * self.annulusAxialIndCorrected[annulusID] + (1 - relaxationFactor) * self.annulusAxialInd[annulusID]
                self.annulusAzimInd[annulusID] = relaxationFactor * self.annulusAzimIndCorrected[annulusID] + (1 - relaxationFactor) * self.annulusAzimInd[annulusID]
                
            #Compute final induction errors
            axialIndError = np.absolute(self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID])
            azimIndError = np.absolute(self.annulusAzimInd[annulusID] - self.annulusAzimIndCorrected[annulusID])
            
            print(f"Annulus {annulusID} solved after {i + 1} iterations!")
            
            if solverLog:
                
                self.solverLog(annulusID, axialIndError, azimIndError)
            
        self.rotorPerformance()
        print(f"{self.rotor.name} CT = {self.rotorCT}")
        print(f"{self.rotor.name} AeroCP = {self.rotorAeroCP}")
        print(f"{self.rotor.name} MechCP = {self.rotorMechCP}")
        
    #Optimization function for the optimization solver. Have to use static method since the scipy optimization requires the first argument of objective function to have the objective variables and with a class method the first argument needs to be referenced to the object
    @staticmethod
    def optzFunction(optzArg, BEMObject, annulusID):
    
        #Since optimizer chooses the input argument randomly, azimuthal induction needs to be adjusted accordingly apriori in order to represent the inflow accurately
        BEMObject.annulusAxialInd[annulusID] = optzArg
        BEMObject.annulusAzimInd[annulusID] = BEMObject.annulusBladeAzimLoad2D[annulusID] / (2 * BEMObject.fstreamRho * (2 * np.pi* BEMObject.rotor.annulusLoc[annulusID] * BEMObject.rotor.rotorRadius) * BEMObject.fstreamVelc**2 * (1 - BEMObject.annulusAxialIndCorrected[annulusID]) * BEMObject.tsr * BEMObject.rotor.annulusLoc[annulusID])
        #Prandtl correction if any
        BEMObject.annulusAzimInd[annulusID] = BEMObject.annulusAzimInd[annulusID]/BEMObject.annulusPrandtlCorrectionFactor[annulusID]
        
        BEMObject.indToBladeLoadSingleDOF(annulusID)
        BEMObject.momLoadToIndSingleDOF(annulusID)
        
        #axialIndError = BEMObject.annulusAxialInd[annulusID] - BEMObject.annulusAxialIndCorrected[annulusID]
        axialIndError = optzArg - BEMObject.annulusAxialIndCorrected[annulusID]
        
        iterationError = np.absolute(axialIndError)
        
        return iterationError
            
    def optzSolver(self, solverLog = False):
        
        for annulusID in range(self.rotor.annuliQuantity):
        
            print(f"Solving for annulus ID = {annulusID} ...")
            
            res = optz.direct(self.optzFunction, bounds = [(0,0.95),], args = (self, annulusID), maxfun=10000, maxiter=10000, len_tol=1e-10)
            print(res.message)
            
            #Compute final induction errors
            axialIndError = np.absolute(self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID])
            azimIndError = np.absolute(self.annulusAzimInd[annulusID] - self.annulusAzimIndCorrected[annulusID])
            
            if solverLog:
                
                self.solverLog(annulusID, axialIndError, azimIndError)
            
        self.rotorPerformance()
        print(f"{self.rotor.name} CT = {self.rotorCT}")
        print(f"{self.rotor.name} AeroCP = {self.rotorAeroCP}")
        print(f"{self.rotor.name} MechCP = {self.rotorMechCP}")        
            
    def rotorPerformance(self):
    
        rotorThrust = np.sum(self.annulusBladeAxialLoad3D)
        rotorAeroPower = np.sum(self.annulusAeroPower)
        rotorMechPower = np.sum(self.annulusMechPower)
        
        self.rotorCT = rotorThrust / (0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.rotorArea)
        self.rotorAeroCP = rotorAeroPower / (0.5 * self.fstreamRho * self.fstreamVelc**3 * self.rotor.rotorArea)
        self.rotorMechCP = rotorMechPower / (0.5 * self.fstreamRho * self.fstreamVelc**3 * self.rotor.rotorArea)
        
    def solverLog(self, annulusID, axialIndError, azimIndError):
        
        print("")
        print(f"Prandtl correction factor = {self.annulusPrandtlCorrectionFactor[annulusID]}")
        print("")
        print(f"Annulus CT from blade element theory = {self.annulusBladeCT[annulusID]}")
        print(f"Annulus CQ from blade element theory = {self.annulusBladeCQ[annulusID]}")
        print(f"Annulus aerodynamic CP from blade element theory = {self.annulusBladeAeroCP[annulusID]}")
        print(f"Annulus mechanical CP from blade element theory = {self.annulusBladeMechCP[annulusID]}")
        print("")
        print(f"Annulus CT from momentum theory = {self.annulusMomCT[annulusID]}")
        print(f"Annulus CQ from momentum theory = {self.annulusMomCQ[annulusID]}")
        print(f"Annulus CP from momentum theory = {self.annulusMomCP[annulusID]}")
        print("")
        print(f"Final annulus axial induction = {self.annulusAxialIndCorrected[annulusID]}")
        print(f"Final annulus azimuthal induction = {self.annulusAzimIndCorrected[annulusID]}")
        print("")
        print(f"Final axial induction error = {axialIndError}")
        print(f"Final azimuthal induction error = {azimIndError}")
        print("===================================================================")
        
    def saveOutput(self):
        
        result = rotorResults()
        
        #Save metadata
        result.rotorName = self.rotor.name
        result.rotorDia = self.rotor.rotorDia
        result.bladeQuanity = self.rotor.bladeQuantity
        result.bladePitch = self.bladePitch
        result.fstreamVelc = self.fstreamVelc
        result.fstreamRho = self.fstreamRho
        result.fstreamPres = self.fstreamPres
        result.tsr = self.tsr
        result.bladeRootLoc = self.rotor.bladeRootLoc
        result.bladeTipLoc = self.rotor.bladeTipLoc
        result.bladeTwistDist = inspect.getsource(self.rotor.bladeTwistDist)[:-1].split(":")[1]#translate lambda expression into string
        result.bladeChordDist = inspect.getsource(self.rotor.bladeChordDist)[:-1].split(":")[1]
        result.annuliQuantity = self.rotor.annuliQuantity
        
        #Save solution
        result.annulusLoc = self.rotor.annulusLoc
        result.annulusArea = self.rotor.annulusArea
        result.annulusBladeAoA = self.annulusBladeAoA
        result.annulusInflow = self.annulusInflow
        result.annulusAxialInd = self.annulusAxialInd
        result.annulusAzimInd = self.annulusAzimInd
        result.annulusCT = self.annulusBladeCT
        result.annulusCQ = self.annulusBladeCQ
        result.annulusCP = self.annulusBladeMechCP
        result.rotorCT = self.rotorCT
        result.rotorCP = self.rotorMechCP
        
        return result
        
        
    
        
            
            
                    
            
            
        
        
        
           
        
        
    
        
        
        
        
        
        
        