import numpy as np

def movingAverage(a, n=2):
    
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def optzFunction(optzArg, BEMobject, annulusID, relaxationFactor):
    
    BEMobject.annulusAxialInd[annulusID] = optzArg
    
    BEMobject.indToBladeLoadSingleDOF(annulusID)
    BEMobject.momLoadToIndSingleDOF(annulusID)
    
    BEMobject.annulusAzimInd = relaxationFactor * BEMobject.annulusAzimIndCorrected + (1 - relaxationFactor) * BEMobject.annulusAzimInd
    
    axialIndError = BEMobject.annulusAxialInd[annulusID] - BEMobject.annulusAxialIndCorrected[annulusID]
    
    return np.absolute(axialIndError)

def optzFunctionv2(optzArgs, BEMobject, annulusID):
    
    BEMobject.annulusAxialInd[annulusID] = optzArgs[0]
    BEMobject.annulusAzimInd[annulusID] = optzArgs[1]
    
    BEMobject.indToBladeLoadSingleDOF(annulusID)
    BEMobject.momLoadToIndSingleDOF(annulusID)
    
    #BEMobject.annulusAzimInd[annulusID] = BEMobject.annulusAzimIndCorrected[annulusID]
    
    axialIndError = BEMobject.annulusAxialInd[annulusID] - BEMobject.annulusAxialIndCorrected[annulusID]    
    azimIndError = BEMobject.annulusAzimInd[annulusID] - BEMobject.annulusAzimIndCorrected[annulusID]
    
    iterationError = np.sqrt(axialIndError**2 + azimIndError**2)
    
    return iterationError

def optzFunctionv3(optzArgs, BEMobject, annulusID):
    
    BEMobject.annulusAxialInd[annulusID] = optzArgs[0]
    BEMobject.annulusAzimInd[annulusID] = optzArgs[1]
    
    BEMobject.indToBladeLoadSingleDOF(annulusID)
    BEMobject.momLoadToIndSingleDOF(annulusID)
    
    axialLoadError = BEMobject.annulusBladeAxialLoad3D[annulusID] - BEMobject.annulusMomAxialLoad3D[annulusID]
    azimLoadError = BEMobject.annulusBladeAzimLoad3D[annulusID] - BEMobject.annulusMomAzimLoad3D[annulusID]
    
    iterationError = np.sqrt(axialLoadError**2 + azimLoadError**2)
    
    return iterationError