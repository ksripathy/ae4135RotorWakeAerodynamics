import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

from src.rotor import rotor
from src.BEM import BEM
import src.dataGen as dataGen

dia1 = 100
numberOfBlades1 = 3
bladePitch1 = -2.5
fstreamVelc1 = 10
fstreamRho1 = 1.225
fstreamPres1 = 1e5
tsr1 = 10
bladeRootLoc1 = 0.2
bladeTipLoc1 = 1.0
bladeTwistDist1 = lambda ndimRadius : 14*(1 - ndimRadius)
bladeChordDist1 = lambda ndimRadius : (3*(1 - ndimRadius)) + 1

#DU95W180 polar
afoilPolar1 = np.loadtxt(homeDir + "/polarDU95W180.csv", delimiter=",")

rotor1 = rotor(dia1, numberOfBlades1, bladeRootLoc1, bladeTipLoc1, bladeTwistDist1, bladeChordDist1, afoilPolar1, 80)
BEMrotor1 = BEM(rotor1, bladePitch1, fstreamVelc1, fstreamRho1, fstreamPres1, tsr1, glauertToggle=True, prandtlToggle=True)

BEMrotor1.classicSolver(printIter=False, solverLog=False)

#BEMrotor1.optzSolver()

BEMrotor1result = BEMrotor1.saveOutput()
BEMrotor1result.saveAsJson(dataDir)
BEMrotor1resultReimport = dataGen.readJSON(dataDir + "/results_Rotor-default_vinf=10_tsr=10_pitch=-2.json")









