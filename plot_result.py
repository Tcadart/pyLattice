import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

def Type_lattice(Lattice):
    if (Lattice == 0):
        Type = 'BCC'
    if (Lattice == 1):
        Type = 'Octet'
    if (Lattice == 2):
        Type = 'OctetExt'
    if (Lattice == 3):
        Type = 'OctetInt'
    if (Lattice == 4):
        Type = 'BCCZ'
    if (Lattice == 5):
        Type = 'Cubic'
    if (Lattice == 6):
        Type = 'OctetInt_without'
    if (Lattice == 7):
        Type = 'Octet_without'
    if (Lattice == 8):
        Type = 'OctahedronZ'
    if (Lattice == 9):
        Type = 'OctahedronZcross'
    if (Lattice == 10):
        Type = 'Kelvin'
    if (Lattice == 11):
        Type = 'CubicV2'
    if (Lattice == 12):
        Type = 'Octet_corr'
    return Type


def openFile(Lattice_Type,number_cell,AnalysisType,MethodSim):
    file_path = "D:/travail_Abaqus/Lattice/"+Type_lattice(Lattice_Type)+"_"+str(number_cell)+str(AnalysisType)+str(MethodSim)+".txt"
    with open(file_path, 'r') as file:
        dataRF = []
        dataU = []
        dataTime = []
        idx = 0
        for line in file:
            linet = np.fromstring(line.strip('[]\n'), sep=',')
            if idx == 0:
                dataRF.append(linet)
            elif idx == 1:
                dataU.append(linet)
            elif idx == 2:
                dataTime.append(linet)
            if ']' in line:
                idx+=1
    return np.concatenate(dataRF),np.concatenate(dataU),np.concatenate(dataTime)

def processDataStressStrain(dataRF, dataU,Number_cell):
    # Modifier les valeurs RF en les divisant par la surface
    dataRF = [rf / (Number_cell * Number_cell) for rf in dataRF]
    
    # Modifier les valeurs U en les multipliant par la longueur du lattice selon Z
    dataU = [-u / (Number_cell) for u in dataU]
    return dataRF, dataU

def plotData(dataRF, dataU,legend):
    
    plt.plot(dataU, dataRF,'--', label = legend) 

    # Ajoutez des titres et des étiquettes
    plt.title('Graphique des Contrainte-Déformation', fontsize=20)
    plt.xlabel('Déformation macroscopique', fontsize=20)
    plt.ylabel('Contrainte macroscopique', fontsize=20)
    
    # Affichez le tracé
    plt.grid(True)


def open_solid_file(chemin_fichier):
    # Initialiser les tableaux pour chaque ensemble de données
    x_values = []
    xydata_1_values = []
    xydata_2_values = []
    current_section = None

    with open(chemin_fichier, 'r') as fichier:
        for ligne in fichier:
            ligne = ligne.strip()
            if not ligne:
                continue
            
            # Identifier la section actuelle
            if "XYData-1" in ligne:
                current_section = 'xydata_1'
                continue
            elif "XYData-2" in ligne:
                current_section = 'xydata_2'
                continue

            # Extraire et convertir les données
            if current_section and ligne[0].isdigit():
                x, value = ligne.split()
                x = float(x)
                value = float(value)

                # Ajouter les valeurs aux tableaux appropriés
                if current_section == 'xydata_1':
                    x_values.append(x)
                    xydata_1_values.append(value)
                elif current_section == 'xydata_2':
                    xydata_2_values.append(value)

    # Convertir les listes en tableaux NumPy
    x_array = np.array(x_values)
    xydata_1_array = np.array(xydata_1_values)
    xydata_2_array = np.array(xydata_2_values)

    return x_array, xydata_1_array, xydata_2_array



number_cell = 4

Lattice_Type = 0
# 0 => BCC
# 1 => Octet 
# 2 => OctetExt 
# 3 => OctetInt 
# 4 => BCCZ 
# 5 => Cubic 
# 6 => OctahedronZ 
# 7 => OctahedronZcross 
# 8 => Kelvin
# 9 => Cubic formulation 2 (centered)
# 10 => Auxetic

AnalysisType = 1
# 0 Modelisation lattice only
# 1 Compression Z static
# 2 Compression Z Implicit
# 3 Compression BCC Solid Plastic Static

MethodSim = 0
# 0 No modification
# 1 Node Modification


# Beam
MethodSim = 0
dataRF,dataU,dataTime = openFile(Lattice_Type,number_cell,AnalysisType,MethodSim)
dataRF, dataU = processDataStressStrain(dataRF, dataU,number_cell)
plotData(dataRF, dataU, 'Beam')

MethodSim = 1
dataRF,dataU,dataTime = openFile(Lattice_Type,number_cell,AnalysisType,MethodSim)
dataRF, dataU = processDataStressStrain(dataRF, dataU,number_cell)
plotData(dataRF, dataU,'BeamMod')

dataTimeSolid, dataRF3Solid, dataU3Solid = open_solid_file("D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/Plasticity/result_solid_"+str(number_cell)+".txt")
dataRF3Solid, dataU3Solid = processDataStressStrain(dataRF3Solid, dataU3Solid,number_cell)
plotData(dataRF3Solid, dataU3Solid,'Solid')

plt.legend(fontsize=15)
plt.show()