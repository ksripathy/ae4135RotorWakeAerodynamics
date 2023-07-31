import numpy as np
import scipy.optimize as optz
from glauertsCorrection import glauertsCorrectionSingleDOF
from miscTools import optzFunction


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
        
        self.annulusAxialInd = axialIndInit * np.ones(rotor.annuliQuantity)
        self.annulusAzimInd = np.zeros(rotor.annuliQuantity)
        
        self.glauertToggle = glauertToggle
        self.prandtlToggle = prandtlToggle
               
        #Initializing dependent attributes    
        self.annulusAxialVelc = np.zeros(rotor.annuliQuantity)
        self.annulusAzimVelc = np.zeros(rotor.annuliQuantity)
        self.annulusResVelc = np.zeros(rotor.annuliQuantity)
        self.annulusInflow = np.zeros(rotor.annuliQuantity)
        self.annulusBladeAoA = np.zeros(rotor.annuliQuantity)
        self.annulusBladeCl = np.zeros(rotor.annuliQuantity)
        self.annulusBladeCd = np.zeros(rotor.annuliQuantity)
        self.annulusBladeLift2D = np.zeros(rotor.annuliQuantity)
        self.annulusBladeDrag2D = np.zeros(rotor.annuliQuantity)
        self.annulusBladeAxialLoad2D = np.zeros(rotor.annuliQuantity)
        self.annulusBladeAzimLoad2D = np.zeros(rotor.annuliQuantity)
        self.annulusBladeAxialLoad3D = np.zeros(rotor.annuliQuantity)
        self.annulusBladeAzimLoad3D = np.zeros(rotor.annuliQuantity)
        self.annulusMomAxialLoad3D = np.zeros(rotor.annuliQuantity)
        self.annulusMomAzimLoad3D = np.zeros(rotor.annuliQuantity)
        self.annulusCT = np.zeros(rotor.annuliQuantity)
        self.annulusAxialIndCorrected = np.zeros(rotor.annuliQuantity)
        self.annulusAzimIndCorrected = np.zeros(rotor.annuliQuantity)
        
        print(f"Momentum parameters initialized for {self.rotor.name}!")
        
    def indToBladeLoadSingleDOF(self, ID):
        
        #Velocities perceived at blade element
        self.annulusAxialVelc[ID] = self.fstreamVelc * (1 - self.annulusAxialInd[ID])
        self.annulusAzimVelc[ID] = self.rotorOmega * (self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * (1 + self.annulusAzimInd[ID])
        self.annulusResVelc[ID] = np.sqrt (self.annulusAxialVelc[ID]**2 + self.annulusAzimVelc[ID]**2)
        
        #Blade element angles
        self.annulusInflow[ID] = np.round(np.arctan2(self.annulusAxialVelc[ID], self.annulusAzimVelc[ID]) * (180/np.pi), 2)
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
        self.annulusBladeAzimLoad3D[ID] = self.annulusBladeAzimLoad2D[ID] * self.rotor.annulusSpan[ID]
        
        self.annulusMomAxialLoad3D[ID] = 2 * self.fstreamRho * self.fstreamVelc**2 * self.annulusAxialInd[ID] * (1 - self.annulusAxialInd[ID]) * self.rotor.annulusArea[ID]
        self.annulusMomAzimLoad3D[ID] = 2 * self.fstreamRho * self.fstreamVelc * (1 - self.annulusAxialInd[ID]) * self.rotorOmega * self.annulusAzimInd[ID] * (self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * self.rotor.annulusArea[ID]
        
        #Annulus thrust coefficieint
        self.annulusCT[ID] = self.annulusBladeAxialLoad3D[ID]/(0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.annulusArea[ID])
        
    def momLoadToIndSingleDOF(self, ID):
        
        if self.glauertToggle:
            
            self.annulusAxialIndCorrected[ID] = glauertsCorrectionSingleDOF(self.annulusCT[ID])
            
        else:
            self.annulusAxialIndCorrected[ID] = 0.5 * (1 - np.sqrt(1 - self.annulusCT[ID]))        
        
        self.annulusAzimIndCorrected[ID] = self.annulusBladeAzimLoad2D[ID] / (2 * self.fstreamRho * (2 * np.pi* self.rotor.annulusLoc[ID] * self.rotor.rotorRadius) * self.fstreamVelc**2 * (1 - self.annulusAxialIndCorrected[ID]) * self.tsr * self.rotor.annulusLoc[ID])
                
    def classicSolver(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3):
        
        for annulusID in range(self.rotor.annuliQuantity):
            
            print(f"Solving for annulus ID = {annulusID} ...")
            
            for i in range(maxIterations):
                
                self.indToBladeLoadSingleDOF(annulusID)
                self.momLoadToIndSingleDOF(annulusID)
                
                inductionError = self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID]
                
                self.annulusAxialInd[annulusID] = relaxationFactor * self.annulusAxialIndCorrected[annulusID] + (1 - relaxationFactor) * self.annulusAxialInd[annulusID]
                self.annulusAzimInd[annulusID] = relaxationFactor * self.annulusAzimIndCorrected[annulusID] + (1 - relaxationFactor) * self.annulusAzimInd[annulusID]
                
                if np.absolute(inductionError) < minError:
                    break
            
            print(f"Annulus {annulusID} solved after {i + 1} iterations!")
            print("")    
            print(f"Annulus axial loading from blade element theory = {np.round(self.annulusBladeAxialLoad3D[annulusID],2)} N")
            print(f"Annnulus azimuthal loading from blade element theory = {np.round(self.annulusBladeAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Annulus axial loading from momentum theory = {np.round(self.annulusMomAxialLoad3D[annulusID],2)} N")
            print(f"Annulus azimuthal loading from momentum theory = {np.round(self.annulusMomAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Final annulus axial induction = {self.annulusAxialInd[annulusID]}")
            print(f"Final annulus azimuthal induction = {self.annulusAzimInd[annulusID]}")
            print("")
            print(f"Final axial induction error = {self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID]}")
            print(f"Final azimuthal induction error = {self.annulusAzimInd[annulusID] - self.annulusAzimIndCorrected[annulusID]}")
            print("===================================================================")
            
    def optzSolver(self):
        
        for annulusID in range(self.rotor.annuliQuantity):
        
            print(f"Solving for annulus ID = {annulusID} ...")
            
            optz.direct(optzFunction, bounds = [(0,0.95), (0,0.95)], args = (self, annulusID))
            
            print("")    
            print(f"Annulus axial loading from blade element theory = {np.round(self.annulusBladeAxialLoad3D[annulusID],2)} N")
            print(f"Annnulus azimuthal loading from blade element theory = {np.round(self.annulusBladeAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Annulus axial loading from momentum theory = {np.round(self.annulusMomAxialLoad3D[annulusID],2)} N")
            print(f"Annulus azimuthal loading from momentum theory = {np.round(self.annulusMomAzimLoad3D[annulusID],2)} N")
            print("")
            print(f"Final annulus axial induction = {self.annulusAxialInd[annulusID]}")
            print(f"Final annulus azimuthal induction = {self.annulusAzimInd[annulusID]}")
            print("")
            print(f"Final axial induction error = {self.annulusAxialInd[annulusID] - self.annulusAxialIndCorrected[annulusID]}")
            print(f"Final azimuthal induction error = {self.annulusAzimInd[annulusID] - self.annulusAzimIndCorrected[annulusID]}")
            print("===================================================================")
    
        
            
            
                    
            
            
        
        
        
           
        
        
    
        
        
        
        
        
        
        