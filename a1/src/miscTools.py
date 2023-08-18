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
    
