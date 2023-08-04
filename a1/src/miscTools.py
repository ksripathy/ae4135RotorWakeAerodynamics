import numpy as np

def movingAverage(a, n=2):
    
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def optzFunction(optzArg, BEMObject, annulusID):
    
    #Since optimizer chooses the input argument randomly, azimuthal induction needs to be adjusted accordingly apriori in order to represent the inflow accurately
    BEMObject.annulusAxialInd[annulusID] = optzArg
    BEMObject.annulusAzimInd[annulusID] = BEMObject.annulusBladeAzimLoad2D[annulusID] / (2 * BEMObject.fstreamRho * (2 * np.pi* BEMObject.rotor.annulusLoc[annulusID] * BEMObject.rotor.rotorRadius) * BEMObject.fstreamVelc**2 * (1 - BEMObject.annulusAxialIndCorrected[annulusID]) * BEMObject.tsr * BEMObject.rotor.annulusLoc[annulusID])
    #Prandtl correction if any
    BEMObject.annulusAzimInd[annulusID] = BEMObject.annulusAzimInd[annulusID]/BEMObject.annulusPrandtlCorrectionFactor[annulusID]
    
    BEMObject.indToBladeLoadSingleDOF(annulusID)
    BEMObject.momLoadToIndSingleDOF(annulusID)
    
    #axialIndError = BEMObject.annulusAxialInd[annulusID] - BEMObject.annulusAxialIndCorrected[annulusID]
    axialIndError = optzArg - BEMObject.annulusAxialIndCorrected[annulusID]
    
    iterationError = np.absolute(axialIndError)
    
    return iterationError

def optzFunction2Args(optzArgs, BEMObject, annulusID):
    
    BEMObject.annulusAxialInd[annulusID] = optzArgs[0]
    BEMObject.annulusAzimInd[annulusID] = optzArgs[1]
    
    BEMObject.indToBladeLoadSingleDOF(annulusID)
    BEMObject.momLoadToIndSingleDOF(annulusID)
    
    axialIndError = (BEMObject.annulusAxialInd[annulusID] - BEMObject.annulusAxialIndCorrected[annulusID]) / BEMObject.annulusAxialIndCorrected[annulusID]    
    azimIndError = (BEMObject.annulusAzimInd[annulusID] - BEMObject.annulusAzimIndCorrected[annulusID]) / BEMObject.annulusAzimIndCorrected[annulusID]
    
    iterationError = np.sqrt(axialIndError**2 + azimIndError**2)

    return iterationError