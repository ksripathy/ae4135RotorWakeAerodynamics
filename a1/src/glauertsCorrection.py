import numpy as np

def glauertsCorrection(annulusCT):
    
    #Glauert's empirical parameters
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    annulusAxialIndCorrected = np.zeros(np.shape(annulusCT)[0])
    
    for i in range(np.shape(annulusCT)[0]):
        
        if annulusCT[i] < CT2:
            
           annulusAxialIndCorrected[i] = 0.5 * (1 - np.sqrt(1 - annulusCT[i]))
           
        else:
            
            annulusAxialIndCorrected[i] = 1 + ((annulusCT[i] - CT1) / (4 * (np.sqrt(CT1) - 1)))
            
    return annulusAxialIndCorrected

def glauertsCorrectionSingleDOF(annulusCT):
    
    #Glauert's empirical parameters
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    
    if annulusCT < CT2:
            
        annulusAxialIndCorrected = 0.5 * (1 - np.sqrt(1 - annulusCT))
        
    else:
        
        annulusAxialIndCorrected = 1 + ((annulusCT - CT1) / (4 * (np.sqrt(CT1) - 1)))
        
    return annulusAxialIndCorrected

def glauertsCorrectionSingleDOFv2(annulusAxialInd):
    
    #Glauert's empirical parameters
    CT1 = 1.816
    CT2 = 2 * np.sqrt(CT1) - CT1
    a1 = 1 - 0.5 * np.sqrt(CT1)
    
    if annulusAxialInd < a1:
        
        annulusCTCorrected = 4 * annulusAxialInd * (1 - annulusAxialInd)
        
    else:
        
        annulusCTCorrected = CT1 - 4 * (np.sqrt(CT1) - 1) * (1 - a1)
        
    if annulusCTCorrected < CT2:
        
        annulusAxialIndCorrected = 0.5 * (1 - np.sqrt(1 - annulusCTCorrected))
        
    else:
        
        annulusAxialIndCorrected = 1 + (annulusCTCorrected - CT1)/(4 * (np.sqrt(CT1) - 1))
        
    return annulusAxialIndCorrected