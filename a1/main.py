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

from src.user import rotor, wind
from src.mesh import mesh
from src.BEM import BEM
from src.ioUtils import saveData

dia1 = 100
numberOfBlades1 = 3
bladePitch1 = -8
fstreamVelc1 = 10
fstreamRho1 = 1.225
fstreamPres1 = 1e5
tsr1 = 8
bladeRootLoc1 = 0.2
bladeTipLoc1 = 1.0
bladeTwistDist1 = lambda ndimRadius : 14*(1 - ndimRadius)
bladeChordDist1 = lambda ndimRadius : (3*(1 - ndimRadius)) + 1

#DU95W180 polar
afoilPolar1 = np.loadtxt(homeDir + "/polarDU95W180.csv", delimiter=",")

rotor1 = rotor(dia1, tsr1, numberOfBlades1, bladeRootLoc1, bladeTipLoc1, bladePitch1, bladeTwistDist1, bladeChordDist1, afoilPolar1)
wind1 = wind(fstreamVelc1, fstreamRho1, fstreamPres1)
mesh11 = mesh(rotor1, wind1, 40)
bemMesh11 = BEM(mesh11)

bemMesh11.classicSolver(printIter=False, solverLog=True)

#BEMrotor1.optzSolver(solverLog=True)

output11 = saveData(mesh11,dataDir)

'''BEMrotor1result = BEMrotor1.saveOutput()
BEMrotor1result.saveAsJson(dataDir)
BEMrotor1resultReimport = solution.readJSON(dataDir + "/results_Rotor-default_vinf=10_tsr=10_pitch=-2.json")'''









