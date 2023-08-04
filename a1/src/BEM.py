import numpy as np
import scipy.optimize as optz
from glauertsCorrection import glauertsCorrectionSingleDOF
from glauertsCorrection import glauertsCorrectionSingleDOFInverse
from prandtlsCorrection import prandtlsCorrectionSingleDOF
from prandtlsCorrection import prandtlsCorrectionSingleDOFv2
from miscTools import optzFunction
from miscTools import optzFunction2Args


class BEM:
    
    def __init__(self, rotor, bladePitch, fstreamVelc, fstreamRho, fstreamPres, tsr, axialIndInit = 0.3, glauertToggle = True, prandtlToggle = True):
        
        self.rotor = rotor
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
        self.annulusMomAxialLoad3D = np.zeros(self.rotor.annuliQuantity)
        self.annulusMomAzimLoad3D = np.zeros(self.rotor.annuliQuantity)
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
        self.rotorCP = 0
        
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
        
        self.annulusMomAxialLoad3D[ID] = 2 * self.fstreamRho * self.fstreamVelc**2 * self.annulusAxialInd[ID] * (1 - self.annulusAxialInd[ID]) * self.rotor.annulusArea[ID]
        
        self.annulusMomAzimLoad3D[ID] = 2 * self.fstreamRho * self.fstreamVelc * (1 - self.annulusAxialInd[ID]) * self.rotorOmega * self.annulusAzimInd[ID] * (self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * self.rotor.annulusArea[ID]
        
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
        
        self.annulusMomCT[ID] = glauertsCorrectionSingleDOFInverse(self.annulusAxialIndCorrected[ID])
        self.annulusMomCQ[ID] = 4 * self.annulusAzimIndCorrected[ID] * (1 - self.annulusAxialIndCorrected[ID]) * self.annulusTSR[ID] * self.rotor.annulusLoc[ID]
        self.annulusMomCP[ID] = 4 * self.annulusAxialIndCorrected[ID] * (1 - self.annulusAxialIndCorrected[ID])**2
        
        '''if self.prandtlToggle:
            
            self.annulusAxialIndCorrected[ID], self.annulusAzimIndCorrected[ID] = prandtlsCorrectionSingleDOF(self.rotor.bladeQuantity, self.rotor.bladeRootLoc, self.rotor.annulusLoc[ID], self.annulusTSR[ID], self.annulusAxialIndCorrected[ID], self.annulusAzimIndCorrected[ID])'''   
        
    def classicSolver(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3, printInd = False, adaptiveRel = False):
        
        for annulusID in range(self.rotor.annuliQuantity):
            
            print(f"Solving for annulus ID = {annulusID} ...")
            
            for i in range(maxIterations):
                
                self.indToBladeLoadSingleDOF(annulusID)
                self.momLoadToIndSingleDOF(annulusID)
                
                iterationError = np.absolute(self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID])
                
                if printInd:
                    
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
            print("")    
            '''print(f"Annulus axial loading from blade element theory = {np.round(self.annulusBladeAxialLoad3D[annulusID],2)} N")
            print(f"Annnulus azimuthal loading from blade element theory = {np.round(self.annulusBladeAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Annulus axial loading from momentum theory = {np.round(self.annulusMomAxialLoad3D[annulusID],2)} N")
            print(f"Annulus azimuthal loading from momentum theory = {np.round(self.annulusMomAzimLoad3D[annulusID],2)} N")'''
            print(f"Annulus CT from blade element theory = {self.annulusBladeCT[annulusID]}")
            print(f"Annulus CQ from blade element theory = {self.annulusBladeCQ[annulusID]}")
            print(f"Annulus CP from blade element theory = {self.annulusBladeCP[annulusID]}")
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
            
    def optzSolver(self):
        
        for annulusID in range(self.rotor.annuliQuantity):
        
            print(f"Solving for annulus ID = {annulusID} ...")
            
            res = optz.direct(optzFunction, bounds = [(0,0.95),], args = (self, annulusID), maxfun=10000, maxiter=10000, len_tol=1e-10)
            #res = optz.direct(optzFunction2Args, bounds = [(0,0.95), (0,0.95)], args = (self, annulusID))
            #res = optz.shgo(optzFunction, bounds = [(0,0.95),], args = (self, annulusID))
            print(res.message)
            
            #Compute final induction errors
            axialIndError = np.absolute(self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID])
            azimIndError = np.absolute(self.annulusAzimInd[annulusID] - self.annulusAzimIndCorrected[annulusID])
                                    
            '''#Pass on correct induction values to determine loads accurately
            self.annulusAxialInd[annulusID] = self.annulusAxialIndCorrected[annulusID]
            self.annulusAzimInd[annulusID] = self.annulusAzimIndCorrected[annulusID]
            
            #Update load computation
            self.indToBladeLoadSingleDOF(annulusID)
            self.momLoadToIndSingleDOF(annulusID)'''
            
            '''print("")    
            print(f"Annulus axial loading from blade element theory = {np.round(self.annulusBladeAxialLoad3D[annulusID],2)} N")
            print(f"Annnulus azimuthal loading from blade element theory = {np.round(self.annulusBladeAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Annulus axial loading from momentum theory = {np.round(self.annulusMomAxialLoad3D[annulusID],2)} N")
            print(f"Annulus azimuthal loading from momentum theory = {np.round(self.annulusMomAzimLoad3D[annulusID],2)} N")'''
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
            
        self.rotorPerformance()
        print(f"{self.rotor.name} CT = {self.rotorCT}")
        print(f"{self.rotor.name} CP = {self.rotorCP}")
        
            
    def rotorPerformance(self):
    
        rotorThrust = np.sum(self.annulusBladeAxialLoad3D)
        rotorPower = np.sum(self.annulusMechPower)
        
        self.rotorCT = rotorThrust / (0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.rotorArea)
        self.rotorCP = rotorPower / (0.5 * self.fstreamRho * self.fstreamVelc**3 * self.rotor.rotorArea)
    
        
            
            
                    
            
            
        
        
        
           
        
        
    
        
        
        
        
        
        
        