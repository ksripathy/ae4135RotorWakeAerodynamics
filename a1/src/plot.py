import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from miscTools import attributeArray
from miscTools import animVortUpdate

class animVortRings:
    
    def __init__(self, VMObj):
        
        self.VMObj = VMObj
        self.animVortRingsArtists = []
        self.colorList = ['r', 'g', 'b', 'c', 'm', 'y']
        
        #Initalise figure artists for all blades
        self.QChordLines = np.empty(self.VMObj.mesh.rotor.bladeQuantity, dtype=np.ndarray)
        self.rootVortLines = np.empty(self.VMObj.mesh.rotor.bladeQuantity, dtype=np.ndarray)
        self.tipVortLines = np.empty(self.VMObj.mesh.rotor.bladeQuantity, dtype=np.ndarray)
                
    def run(self):
        
        #Empty figure
        self.animVortFig, self.animVortAx = plt.subplots(subplot_kw={"projection":"3d"})
        
        
        
        for j in range(self.VMObj.mesh.rotor.bladeQuantity):
            
            bladeObj = getattr(self.VMObj.mesh, f"annulusBlade{j}")
            
            #Retreive coordinate locations to plot
            QChordPos = attributeArray(bladeObj.QChordHistory[0], 'controlPointLoc')
            rootVortPos = attributeArray(bladeObj.rootVortHistory[0], 'controlPointLoc')
            tipVortPos = tipVortPos = attributeArray(bladeObj.tipVortHistory[0], 'controlPointLoc')
            
            #Set figure artists for each blade. Plot function generates a tuple, But the animation function update works only with an plot element. Hence the first element of a plot tuple is assigned
            self.QChordLines[j] = self.animVortAx.plot(QChordPos[:,0], QChordPos[:,1], QChordPos[:,2], color='k')[0]
            self.rootVortLines[j] = self.animVortAx.plot(rootVortPos[:,0], rootVortPos[:,1], rootVortPos[:,2], color=self.colorList[j])[0]
            self.tipVortLines[j] = self.animVortAx.plot(tipVortPos[:,0], tipVortPos[:,1], tipVortPos[:,2], color=self.colorList[j])[0]
            #Axis limits are not updated each. Hence set the wake limit in advance
            self.animVortAx.set(xlim=[0,self.VMObj.wakeCoverage])
            
        ani = animation.FuncAnimation(fig=self.animVortFig, func=animVortUpdate, frames=range(self.VMObj.totalTimeSteps), fargs=(self,), interval=30)
        
        return ani
            
            
            
            