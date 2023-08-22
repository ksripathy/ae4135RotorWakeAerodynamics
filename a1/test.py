import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

from src.user import rotor, wind
from src.mesh import mesh
from src.vortexRing import vortexRing
from src.miscTools import attributeArray
from src.VM import VM
from src.plot import animVortRings

dia1 = 100
numberOfBlades1 = 3
bladePitch1 = -2
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
mesh11 = mesh(rotor1, wind1, 20)

#Build blades for vortex model
mesh11.buildBlades()

#Parameters for vortex model
wakeCoverage = 2 * dia1
wakeVelc = np.array([fstreamVelc1, 0, 0])
totalTime = wakeCoverage/wakeVelc[0]
deltaTime = 0.1
totalTimeSteps = int(totalTime/deltaTime)

#Initialize vortex model simulation
VMmesh11 = VM(mesh11, wakeCoverage, wakeVelc, deltaTime)

#Build vortex rings for all blades
VMmesh11.buildVortRings()


plots = animVortRings(VMmesh11)
anim = plots.run()

'''fig1, ax1 = plt.subplots(subplot_kw={"projection":"3d"})

b0Qchord = np.empty((mesh11.annuliQuantity,3), dtype=np.ndarray)
b1Qchord = np.empty((mesh11.annuliQuantity,3), dtype=np.ndarray)
b2Qchord = np.empty((mesh11.annuliQuantity,3), dtype=np.ndarray)

for i in range(0,mesh11.annuliQuantity):
    b0Points = attributeArray(mesh11.annulusBlade0.vortSys[i],'controlPointLoc')
    b1Points = attributeArray(mesh11.annulusBlade1.vortSys[i],'controlPointLoc')
    b2Points = attributeArray(mesh11.annulusBlade2.vortSys[i],'controlPointLoc')
    
    b0Qchord[i,0], b0Qchord[i,1], b0Qchord[i,2] = mesh11.annulusBlade0.vortSys[i][0].controlPointLoc
    b1Qchord[i,0], b1Qchord[i,1], b1Qchord[i,2] = mesh11.annulusBlade1.vortSys[i][0].controlPointLoc
    b2Qchord[i,0], b2Qchord[i,1], b2Qchord[i,2] = mesh11.annulusBlade2.vortSys[i][0].controlPointLoc
    
    #Plot vortex filament for root and tip annulus
    if i in (0, 19, 39, 59, 79):
        ax1.plot(b0Points[:,0],b0Points[:,1],b0Points[:,2], color='r')
        ax1.plot(b1Points[:,0],b1Points[:,1],b1Points[:,2], color='g')
        ax1.plot(b2Points[:,0],b2Points[:,1],b2Points[:,2], color='b')
    
ax1.plot(b0Qchord[:,0], b0Qchord[:,1], b0Qchord[:,2], color='k')
ax1.plot(b1Qchord[:,0], b1Qchord[:,1], b1Qchord[:,2], color='k')
ax1.plot(b2Qchord[:,0], b2Qchord[:,1], b2Qchord[:,2], color='k')
        
plt.show()'''

'''#Initialize numpy array of vortex rings
blade1VortexRings = np.empty(totalTimeSteps + 1, dtype=np.ndarray)
#Bound vortex ring definition
blade1VortexRings[0] = vortexRing(mesh11.annulusLocRoot[19], mesh11.annulusLocTip[19], mesh11.annulusBladeChordRoot[19], mesh11.annulusBladeChordTip[19], mesh11.annulusBladeTwistRoot[19], mesh11.annulusBladeTwistTip[19], mesh11.rotor.radius, 0)
#Trailing vortex ring definition
for i in range(1,len(blade1VortexRings)):
    blade1VortexRings[i] = vortexRing()

#Time marching
for i in range(totalTimeSteps):

    #Convect vortex ring
    for j in range(i,-1,-1):
        
        blade1VortexRings[j].convect(blade1VortexRings[j+1], wakeVelc, deltaTime)
        
    #Rotate turbine and hence the bound vortex ring locations
    blade1VortexRings[0].rotate(mesh11.rotorOmega, deltaTime)
    
controlPoints = attributeArray(blade1VortexRings, 'controlPointLoc')
leadRoot = attributeArray(blade1VortexRings, 'leadRootLoc')
leadTip = attributeArray(blade1VortexRings, 'leadTipLoc')
trailRoot = attributeArray(blade1VortexRings, 'trailRootLoc')
trailTip = attributeArray(blade1VortexRings, 'trailTipLoc')

#plotting
fig1, ax1 = plt.subplots(subplot_kw={"projection":"3d"})
#attributeArray(blade1VortexRings, 'plotRings', (ax1,))
#ax1.plot(controlPoints[:,0],controlPoints[:,1],controlPoints[:,2])
ax1.scatter(controlPoints[:,0],controlPoints[:,1],controlPoints[:,2], marker='.')
#ax1.scatter(leadRoot[:,0], leadRoot[:,1], leadRoot[:,2], marker='.')
#ax1.scatter(trailRoot[:,0], trailRoot[:,1], trailRoot[:,2],marker='.')
#ax1.scatter(leadTip[:,0], leadTip[:,1], leadTip[:,2], marker='.')
#ax1.scatter(trailTip[:,0], trailTip[:,1], trailTip[:,2], marker='.')

plt.show()'''

        
    
    