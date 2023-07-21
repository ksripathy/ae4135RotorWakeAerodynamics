import os
import sys

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")

#Add src folder to python path
sys.path.append(srcDir)

from src.rotor import rotor

dia1 = 50
numberOfBlades1 = 3
bladePitch1 = -2
fstreamVelc1 = 10
tsr1 = 8
bladeTwistDist1 = lambda ndimRadius : 14*(1 - ndimRadius)
bladeChordDist1 = lambda ndimRadius : (3*(1 - ndimRadius)) + 1

rotor1 = rotor(dia1, numberOfBlades1, fstreamVelc1, tsr1, bladePitch1, bladeTwistDist1, bladeChordDist1)


