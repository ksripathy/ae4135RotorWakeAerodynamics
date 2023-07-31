import numpy as np

def movingAverage(a, n=2):
    
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def optzFunction(optzArgs, BEMObject, annulusID):
    
    BEMObject.annulusAxialInd[annulusID] = optzArgs[0]
    BEMObject.annulusAzimInd[annulusID] = optzArgs[1]
    
    BEMObject.indToBladeLoadSingleDOF(annulusID)
    BEMObject.momLoadToIndSingleDOF(annulusID)
    
    axialIndError = BEMObject.annulusAxialInd[annulusID] - BEMObject.annulusAxialIndCorrected[annulusID]    
    azimIndError = BEMObject.annulusAzimInd[annulusID] - BEMObject.annulusAzimIndCorrected[annulusID]
    
    iterationError = np.sqrt(axialIndError**2 + azimIndError**2)

    return iterationError