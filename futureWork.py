import numpy as np
import scipy.optimize as optz
from glauertsCorrection import glauertsCorrection
from glauertsCorrection import glauertsCorrectionSingleDOF
from glauertsCorrection import glauertsCorrectionSingleDOFv2
from miscTools import optzFunction


class BEM:
    
    def __init__(self, rotor, bladePitch, fstreamVelc, fstreamRho, fstreamPres, tsr, axialIndInit = 0.3, glauertToggle = True, prandtlToggle = True, singleDOFMode = False):
        
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
        self.singleDOFMode = singleDOFMode
        
        if singleDOFMode:
            
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
        
    def indToBladeLoad(self):
        
        #Velocities perceived at blade element
        self.annulusAxialVelc = self.fstreamVelc * (1 - self.annulusAxialInd)
        self.annulusAzimVelc = self.rotorOmega * (self.rotor.annulusLoc * self.rotor.rotorRadius) * (1 + self.annulusAzimInd)
        self.annulusResVelc = np.sqrt (self.annulusAxialVelc**2 + self.annulusAzimVelc**2)
        
        #Blade element angles
        self.annulusInflow = np.arctan2(self.annulusAxialVelc, self.annulusAzimVelc) * (180/np.pi)
        self.annulusBladeAoA = self.annulusInflow - (self.rotor.annulusBladeTwist + self.bladePitch)
         
        #Blade lift and drag
        self.annulusBladeCl = np.interp(self.annulusBladeAoA, self.rotor.afoilPolar[:,0], self.rotor.afoilPolar[:,1])
        self.annulusBladeCd = np.interp(self.annulusBladeAoA, self.rotor.afoilPolar[:,0], self.rotor.afoilPolar[:,2])
        
        self.annulusBladeLift2D = 0.5 * self.rotor.annulusBladeChord * self.fstreamRho * (self.annulusResVelc**2) * self.annulusBladeCl
        self.annulusBladeDrag2D = 0.5 * self.rotor.annulusBladeChord * self.fstreamRho * (self.annulusResVelc**2) * self.annulusBladeCd
        
        #Annulus normal and tangential laoding
        self.annulusBladeAxialLoad2D = (self.annulusBladeLift2D * np.cos(self.annulusInflow * (np.pi/180)) + self.annulusBladeDrag2D * np.sin(self.annulusInflow * (np.pi/180))) * self.rotor.bladeQuantity
        self.annulusBladeAzimLoad2D = (self.annulusBladeLift2D * np.sin(self.annulusInflow * (np.pi/180)) - self.annulusBladeDrag2D * np.cos(self.annulusInflow * (np.pi/180))) * self.rotor.bladeQuantity
        
        self.annulusBladeAxialLoad3D = self.annulusBladeAxialLoad2D * self.rotor.annulusSpan
        self.annulusBladeAzimLoad3D = self.annulusBladeAzimLoad2D * self.rotor.annulusSpan
        
        self.annulusMomAxialLoad3D = 2 * self.fstreamRho * self.fstreamVelc**2 * self.annulusAxialInd * (1 - self.annulusAxialInd) * self.rotor.annulusArea
        self.annulusMomAzimLoad3D = 2 * self.fstreamRho * self.fstreamVelc * (1 - self.annulusAxialInd) * self.rotorOmega * self.annulusAzimInd * (self.rotor.annulusLoc * self.rotor.rotorRadius) * self.rotor.annulusArea
        
        #Annulus thrust coefficieint
        self.annulusCT = self.annulusBladeAxialLoad3D/(0.5 * self.fstreamRho * self.fstreamVelc**2 * self.rotor.annulusArea)
        
    def momLoadToInd(self):
        
        if self.glauertToggle:
            
            self.annulusAxialIndCorrected = glauertsCorrection(self.annulusCT)
            
        else:
            self.annulusAxialIndCorrected = 0.5 * (1 - np.sqrt(1 - self.annulusCT))        
        
        self.annulusAzimIndCorrected = self.annulusBladeAzimLoad2D / (2 * self.fstreamRho * (2 * np.pi* self.rotor.annulusLoc * self.rotor.rotorRadius) * self.fstreamVelc**2 * (1 - self.annulusAxialIndCorrected) * self.tsr * self.rotor.annulusLoc)
        
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
        
        
    def optimizationSolver(self):
        
        optzBounds = []
        
        #Setting bound tuples for axial induction
        for i in range(self.rotor.annuliQuantity):
            
            optzBounds.append((0, 0.95))
            
        #Setting bound tuples for azimuthal induction
        for i in range(self.rotor.annuliQuantity):
            
            optzBounds.append((0, 0.95))
        
        '''axialIndLbound = np.zeros(self.rotor.annuliQuantity)
        azimIndLbound = np.zeros(self.rotor.annuliQuantity)
        
        axialIndUbound = 0.95 * np.ones(self.rotor.annuliQuantity)
        azimIndUbound = 0.95 * np.ones(self.rotor.annuliQuantity)'''
        
        optz.direct(optzFunction, bounds = optzBounds, args = (self,))
        
        print("\n")
        print(f"Rotor axial loading from blade element theory = {np.round(np.sum(self.annulusBladeAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from blade element theory = {np.round(np.sum(self.annulusBladeAzimLoad3D),2)} N")
        print("\n")
        print(f"Rotor axial loading from momentum theory = {np.round(np.sum(self.annulusMomAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from momentum theory = {np.round(np.sum(self.annulusMomAzimLoad3D),2)} N")
        print("\n")
        print(f"Final axial induction distribution = {self.annulusAxialInd}")
        print(f"Final azimuthal induction distribution = {self.annulusAzimInd}")
        print("\n")
        print("BEM solved :D")
        
    def classicSolver(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3):
        
        for i in range(maxIterations):
        
            print(f"Iteration-{i}:")
            
            self.indToBladeLoad()
            
            print(f"Rotor axial loading from blade element theory = {np.round(np.sum(self.annulusBladeAxialLoad3D),2)} N")
            print(f"Rotor azimuthal loading from blade element theory = {np.round(np.sum(self.annulusBladeAzimLoad3D),2)} N")
            print("\n")
            print(f"Rotor axial loading from momentum theory = {np.round(np.sum(self.annulusMomAxialLoad3D),2)} N")
            print(f"Rotor azimuthal loading from momentum theory = {np.round(np.sum(self.annulusMomAzimLoad3D),2)} N")
            
            self.momLoadToInd()
            
            print("\n")
            
            iterationError = self.annulusAxialInd - self.annulusAxialIndCorrected
            normError = np.linalg.norm(iterationError)/self.rotor.annuliQuantity
            
            print(f"Averaged norm of iteration error = {normError}")
            print(f"Max error in induction is {np.max(np.absolute(iterationError))} at annulus ID {np.argmax(np.absolute(iterationError))}")
            
            #Induction values for next iteration
            self.annulusAxialInd = relaxationFactor * self.annulusAxialIndCorrected + (1 - relaxationFactor) * self.annulusAxialInd
            self.annulusAzimInd = relaxationFactor * self.annulusAzimIndCorrected + (1 - relaxationFactor) * self.annulusAzimInd
            
            if np.max(np.absolute(iterationError)) < minError:
                break
            
            print("=================================================================================================")

        print("End of iterations!")
        print("\n")
        print(f"Rotor axial loading from blade element theory = {np.round(np.sum(self.annulusBladeAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from blade element theory = {np.round(np.sum(self.annulusBladeAzimLoad3D),2)} N")
        print("\n")
        print(f"Rotor axial loading from momentum theory = {np.round(np.sum(self.annulusMomAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from momentum theory = {np.round(np.sum(self.annulusMomAzimLoad3D),2)} N")
        print("\n")
        print(f"Final axial induction distribution = {self.annulusAxialInd}")
        print(f"Final azimuthal induction distribution = {self.annulusAzimInd}")
        print("\n")
        print("BEM solved :D")
        
    def optimizationFunction(self, optimizationArg):
        
        print(f"New Iteration:")
        print(f"Optimization Argument = {optimizationArg}")
        
        self.annulusAxialInd = optimizationArg
            
        self.indToBladeLoad()
        
        '''print(f"Rotor axial loading from blade element theory = {np.round(np.sum(self.annulusBladeAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from blade element theory = {np.round(np.sum(self.annulusBladeAzimLoad3D),2)} N")
        print("\n")
        print(f"Rotor axial loading from momentum theory = {np.round(np.sum(self.annulusMomAxialLoad3D),2)} N")
        print(f"Rotor azimuthal loading from momentum theory = {np.round(np.sum(self.annulusMomAzimLoad3D),2)} N")'''
        
        self.momLoadToInd()
        
        '''print("\n")'''
        
        iterationError = self.annulusAxialInd - self.annulusAxialIndCorrected
        normError = np.linalg.norm(iterationError)/self.rotor.annuliQuantity
        
        '''print(f"Averaged norm of iteration error = {normError}")
        print(f"Max error in induction is {np.max(np.absolute(iterationError))} at annulus ID {np.argmax(np.absolute(iterationError))}")'''
        print(f"Max iteration error = {np.max(np.absolute(iterationError))}")
        print("====================================================================")
        
        return np.max(np.absolute(iterationError))
        
    def solveBEMv2(self):
        
        optz.minimize(self.optimizationFunction, 0.75 * np.ones(self.rotor.annuliQuantity))
        
    def solveBEMv3(self):
        
        optz.shgo(self.optimizationFunction, bounds=optz.Bounds(0.01 * np.ones(self.rotor.annuliQuantity),0.95 * np.ones(self.rotor.annuliQuantity)))
        
    def solveBEMv4(self, minError = 1e-5, maxIterations = 100, relaxationFactor = 0.3):
        
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
            print("===================================================================")
            
    def optzSolver(self):
        
        for annulusID in range(self.rotor.annuliQuantity):
        
            print(f"Solving for annulus ID = {annulusID} ...")
            
            #optz.shgo(optzFunctionv2, bounds = [(0,0.95), [0,0.95]], args = (self, annulusID), options = {"f_tol" : 1e-5, "disp" : False})
            #optz.basinhopping(optzFunctionv2, [self.axialIndInit, 0], minimizer_kwargs = {"args" : (self, annulusID)})
            optz.direct(optzFunction, bounds = [(0,0.95), (0,0.95)], args = (self, annulusID))
            #optz.dual_annealing(optzFunctionv2, bounds = [(0,0.95), [0,0.95]], args = (self, annulusID))
            
            '''#Update final induction to calculate right loads
            self.annulusAxialInd = self.annulusAxialIndCorrected
            self.annulusAzimInd = self.annulusAzimIndCorrected
            self.indToBladeLoadSingleDOF(annulusID)'''
            
            
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
    
        
            
            
                    
            
            
        
        
        
           
        
        
    
        
        
        
        
        
        
        