import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pandas as pd
from Lattice_description import *


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

def openFile_all_data(Lattice_Type, number_cell, AnalysisType, MethodSim):
    file_path = f"D:/travail_Abaqus/Lattice/{Type_lattice(Lattice_Type)}_{number_cell}{AnalysisType}{MethodSim}.txt"
    dataRF = {'RF1': [], 'RF2': [], 'RF3': []}
    dataU = {'U1': [], 'U2': [], 'U3': []}
    dataTime = []
    current_data_list = []

    # Déterminer l'ordre dans lequel les données apparaissent
    data_order = ['RF1', 'RF2', 'RF3', 'U1', 'U2', 'U3', 'dataTime']
    current_index = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Accumuler les données numériques
            current_data_list.extend(np.fromstring(line.strip('[]\n'), sep=','))
            if ']' in line:  # Vérifier la fin de la catégorie actuelle
                # Attribuer les données accumulées à la catégorie correspondante
                if current_index < 3:  # RF1, RF2, RF3
                    dataRF[data_order[current_index]] = np.array(current_data_list)
                elif current_index < 6:  # U1, U2, U3
                    dataU[data_order[current_index]] = np.array(current_data_list)
                else:  # dataTime
                    dataTime = np.array(current_data_list)
                # Réinitialiser pour la prochaine catégorie de données
                current_data_list = []
                current_index += 1

    return dataRF, dataU, dataTime



def processDataStressStrain(dataRF, dataU,Number_cell,length_cell):
    # Modifier les valeurs RF en les divisant par la surface
    dataRF = [rf / (Number_cell * length_cell * Number_cell * length_cell) for rf in dataRF]
    
    # Modifier les valeurs U en les multipliant par la longueur du lattice selon Z
    dataU = [-u / (Number_cell * length_cell) for u in dataU]
    return dataRF, dataU


def processDataStressStrain_all_data(dataRF, dataU, Number_cell, length_cell):
    surface = (Number_cell * length_cell) ** 2
    length = Number_cell * length_cell

    for key in dataRF:
        dataRF[key] = dataRF[key] / surface

    for key in dataU:
        dataU[key] = -dataU[key] / length

    return dataRF, dataU


def plotData(dataRF, dataU,legend):
    
    plt.plot(dataU, dataRF,'--', label = legend)


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


def plotDataExperience(number_cell, length_cell):
    number_cell = 8
    liste_feuille = ['BCC_1','BCC_2','BCC_3']
    df = pd.read_excel('D:/Fichiers/4_Experimentation/Compression lattice/compression_lattice.xlsx', sheet_name=liste_feuille[0])
    Force = df.iloc[:, 1]
    Displacement = df.iloc[:, 2]
    dataRF, dataU = processDataStressStrain(Force, -Displacement, number_cell, length_cell)
    plotData(dataRF, dataU, 'Experiments')



number_cell = 3
length_cell = 7.5

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

MethodSim = 1
# 0 No modification
# 1 Node Modification


# Beam
# MethodSim = 0
# dataRF,dataU,dataTime = openFile(Lattice_Type,number_cell,AnalysisType,MethodSim)
# dataRF, dataU = processDataStressStrain(dataRF, dataU,number_cell)
# plotData(dataRF, dataU, 'Beam')

dataRF,dataU,dataTime = openFile_all_data(Lattice_Type,number_cell,AnalysisType,MethodSim)
print(dataRF)
print(dataU)
dataRF, dataU = processDataStressStrain_all_data(dataRF, dataU,number_cell, length_cell)

plotData(dataRF['RF1'], dataU['U3'],'RF1')
plotData(dataRF['RF2'], dataU['U3'],'RF2')
plotData(dataRF['RF3'], dataU['U3'],'RF3')


# dataTimeSolid, dataRF3Solid, dataU3Solid = open_solid_file("D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/Plasticity/result_solid_"+str(number_cell)+".txt")
# dataRF3Solid, dataU3Solid = processDataStressStrain(dataRF3Solid, dataU3Solid,number_cell)
# plotData(dataRF3Solid, dataU3Solid,'Solid')

# plotDataExperience(number_cell, length_cell)

# Ajoutez des titres et des étiquettes
plt.title('Graphique des Contrainte-Déformation', fontsize=20)
plt.xlabel('Déformation macroscopique', fontsize=20)
plt.ylabel('Contrainte macroscopique', fontsize=20)

# Affichez le tracé
plt.grid(True)
plt.legend(fontsize=15)
plt.show()