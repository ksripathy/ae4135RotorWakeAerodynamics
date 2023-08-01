import numpy as np

def glauertsCorrectionSingleDOF(annulusCT):
    
    #Glauert's empirical parameters
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    if annulusCT < CT2:
            
        annulusAxialIndCorrected = 0.5 * (1 - np.sqrt(1 - annulusCT))
        
    else:
        
        annulusAxialIndCorrected = 1 + ((annulusCT - CT1) / (4 * (np.sqrt(CT1) - 1)))
        
    return annulusAxialIndCorrected

def glauertsCorrectionSingleDOFInverse(annulusAxialInd):
    
    #Glauert's empirical parameters
    CT1 = 1.816
    filtParam = 1 - 0.5 * np.sqrt(CT1)
    
    if annulusAxialInd < filtParam:
        
        annulusCTCorrected = 4 * annulusAxialInd * (1 - annulusAxialInd)
        
    else:
        
        annulusCTCorrected = CT1 - 4 * (np.sqrt(CT1) - 1) * (1 - annulusAxialInd)
        
    return annulusCTCorrected