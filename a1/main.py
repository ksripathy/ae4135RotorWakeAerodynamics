import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")

#Add src folder to python path
sys.path.append(srcDir)

from src.rotor import rotor
from src.BEM import BEM

dia1 = 50
numberOfBlades1 = 3
bladePitch1 = -2
fstreamVelc1 = 10
fstreamRho1 = 1.225
fstreamPres1 = 1e5
tsr1 = 6
bladeRootLoc1 = 0.2
bladeTipLoc1 = 1.0
bladeTwistDist1 = lambda ndimRadius : 14*(1 - ndimRadius)
bladeChordDist1 = lambda ndimRadius : (3*(1 - ndimRadius)) + 1

#DU95W180 polar
afoilPolar1 = np.loadtxt(homeDir + "/polarDU95W180.csv", delimiter=",")

rotor1 = rotor(dia1, numberOfBlades1, bladeRootLoc1, bladeTipLoc1, bladeTwistDist1, bladeChordDist1, afoilPolar1, 40)
BEMrotor1 = BEM(rotor1, bladePitch1, fstreamVelc1, fstreamRho1, fstreamPres1, tsr1, glauertToggle=True, prandtlToggle= True)

#BEMrotor1.classicSolver()

BEMrotor1.optzSolver()








