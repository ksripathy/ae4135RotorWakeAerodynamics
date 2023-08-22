import numpy as np

def movingAverage(a, n=2):
    
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#Function for obtaining array of attribute of multiple objects of same class organized as an array
def attributeArray(objectArray, attribute, args=None):
    '''Returns the requested attribute of all the objects of numpy array'''
    
    if args == None:#Requested attribute is a variable
        
        return np.array([getattr(obj, attribute) for obj in objectArray])
    
    else:#Requested attribute is a method
        
        return np.array([getattr(obj, attribute)(*args) for obj in objectArray])
    
#Update function for vortRing animation
def animVortUpdate(frame, *fargs):
        
        plotObj = fargs[0]
        
        funcReturnList = []
        
        for j in range(plotObj.VMObj.mesh.rotor.bladeQuantity):
            
            bladeObj = getattr(plotObj.VMObj.mesh, f"annulusBlade{j}")
            
            #Retreive coordinate locations to plot
            QChordPos = attributeArray(bladeObj.QChordHistory[frame], 'controlPointLoc')
            rootVortPos = attributeArray(bladeObj.rootVortHistory[frame], 'controlPointLoc')
            tipVortPos = tipVortPos = attributeArray(bladeObj.tipVortHistory[frame], 'controlPointLoc')
            
            plotObj.QChordLines[j].set_xdata(QChordPos[:,0])
            plotObj.QChordLines[j].set_ydata(QChordPos[:,1])
            plotObj.QChordLines[j].set_3d_properties(QChordPos[:,2])
            funcReturnList.append(plotObj.QChordLines[j])
            
            plotObj.rootVortLines[j].set_xdata(rootVortPos[:,0])
            plotObj.rootVortLines[j].set_ydata(rootVortPos[:,1])
            plotObj.rootVortLines[j].set_3d_properties(rootVortPos[:,2])
            funcReturnList.append(plotObj.rootVortLines[j])
            
            plotObj.tipVortLines[j].set_xdata(tipVortPos[:,0])
            plotObj.tipVortLines[j].set_ydata(tipVortPos[:,1])
            plotObj.tipVortLines[j].set_3d_properties(tipVortPos[:,2])
            funcReturnList.append(plotObj.tipVortLines[j])
            
        return tuple(funcReturnList)
            
            
            
            
            
    
        
