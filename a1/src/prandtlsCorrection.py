import numpy as np

def prandtlsCorrectionSingleDOF(bladeQuantity, bladeRootLoc, annulusLoc, annulusTSR, axialInd, azimInd):
    
    B = bladeQuantity
    muRoot = bladeRootLoc
    mu = annulusLoc
    tsr = annulusTSR
    a = axialInd
    
    #Prandtl correction parameters
    const = -0.5 * B * np.sqrt(1 + (tsr**2 * mu**2 / (1 - a)**2))
    fTip = 2/np.pi * np.arccos(np.exp((1 - mu) * const / mu))
    fRoot = 2/np.pi * np.arccos(np.exp((mu - muRoot) * const / mu))
    fTotal = fTip * fRoot
    
    axialIndCorrected = axialInd / fTotal
    azimIndCorrected = azimInd / fTotal
    
    return axialIndCorrected, azimIndCorrected
    