from solution import solution
import inspect
import numpy as np

def saveData(meshObj, dataDir):
    
    output = solution()
    
    #Save metadata
    output.rotorName = meshObj.rotor.name
    output.rotorDia = meshObj.rotor.dia
    output.bladeQuanity = meshObj.rotor.bladeQuantity
    output.bladePitch = meshObj.rotor.bladePitch
    output.fstreamVelc = meshObj.wind.fstreamVelc
    output.fstreamRho = meshObj.wind.fstreamRho
    output.fstreamPres = meshObj.wind.fstreamPres
    output.tsr = meshObj.rotor.tsr
    output.bladeRootLoc = meshObj.rotor.bladeRootLoc
    output.bladeTipLoc = meshObj.rotor.bladeTipLoc
    output.bladeTwistDist = inspect.getsource(meshObj.rotor.bladeTwistDist)[:-1].split(":")[1]#translate lambda expression into string
    output.bladeChordDist = inspect.getsource(meshObj.rotor.bladeChordDist)[:-1].split(":")[1]
    output.annuliQuantity = meshObj.annuliQuantity
    
    #Save solution
    output.annulusLoc = meshObj.annulusLoc
    output.annulusArea = meshObj.annulusArea
    output.annulusBladeAoA = meshObj.annulusBladeAoA
    output.annulusInflow = meshObj.annulusInflow
    output.annulusAxialInd = meshObj.annulusAxialInd
    output.annulusAzimInd = meshObj.annulusAzimInd
    output.annulusCT = meshObj.annulusBladeCT
    output.annulusCQ = meshObj.annulusBladeCQ
    output.annulusCP = meshObj.annulusBladeMechCP
    output.rotorCT = np.sum(meshObj.annulusBladeAxialLoad3D) / (0.5 * meshObj.wind.fstreamRho * meshObj.wind.fstreamVelc**2 * meshObj.rotor.area)
    output.rotorCP = np.sum(meshObj.annulusMechPower) / (0.5 * meshObj.wind.fstreamRho * meshObj.wind.fstreamVelc**3 * meshObj.rotor.area)
    
    output.saveAsJson(dataDir)
    
    return output

def loadData(fileLoc):
    
    with open(fileLoc, 'r') as jsonFile:
        
        input = solution.from_json(jsonFile.read())
        #For JSON serialization numpy arrays are translated into list. Hence after deserialization lists are converted back to numpy arrays
        input.annulusLoc = np.asarray(input.annulusLoc)
        input.annulusArea = np.asarray(input.annulusArea)
        input.annulusBladeAoA = np.asarray(input.annulusBladeAoA)
        input.annulusInflow = np.asarray(input.annulusInflow)
        input.annulusAxialInd = np.asarray(input.annulusAxialInd)
        input.annulusAzimInd = np.asarray(input.annulusAzimInd)
        input.annulusCT = np.asarray(input.annulusCT)
        input.annulusCQ = np.asarray(input.annulusCQ)
        input.annulusCP = np.asarray(input.annulusCP)
        
    return input