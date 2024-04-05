import odbAccess
from odbAccess import openOdb
import sys
import os

import numpy as np
import time
from Lattice import *
import re
from Materials import *

from Lattice_description import *


def get_result(name_Job, name_step):
    odb = openOdb(name_Job + '.odb')
    historyRegions = odb.steps[name_step].historyRegions
    dataRF = {'RF1': [], 'RF2': [], 'RF3': []}
    dataU = {'U1': [], 'U2': [], 'U3': []}
    dataTime = []
    for key in historyRegions.keys():
        historyRegion = historyRegions[key]
        # Traitement pour RF1, RF2, RF3
        for rf_key in dataRF.keys():
            try:
                rf_data = historyRegion.historyOutputs[rf_key].data
                if not dataRF[rf_key]:  # Si la liste est vide
                    dataRF[rf_key] = [value for _, value in rf_data]
                else:
                    for i, (_, value) in enumerate(rf_data):
                        dataRF[rf_key][i] += value
            except KeyError:
                pass
        # Traitement pour U1, U2, U3
        for u_key in dataU.keys():
            try:
                u_data = historyRegion.historyOutputs[u_key].data
                if not dataU[u_key]:  # Si la liste est vide
                    dataU[u_key] = [value for _, value in u_data]
                    if u_key == 'U1':  # Assurer que dataTime est seulement rempli une fois
                        dataTime = [time for time, _ in u_data]
            except KeyError:
                pass
    return dataRF, dataU, dataTime

def save_result(Lattice_Type, number_cell, AnalysisType, MethodSim, dataRF, dataU, dataTime):
    Radius = "0_5"
    filename = Type_lattice(Lattice_Type)+"_"+str(number_cell)+str(AnalysisType)+str(MethodSim)+Radius+".txt"
    print(filename)
    with open(filename, "w") as f:
        for rf_key in ['RF1', 'RF2', 'RF3']:
            rf_data_str = np.array2string(np.array(dataRF[rf_key]), separator=',')
            f.write(rf_data_str + "\n")
        for u_key in ['U1', 'U2', 'U3']:
            u_data_str = np.array2string(np.array(dataU[u_key]), separator=',')
            f.write(u_data_str + "\n")
        time_data_str = np.array2string(np.array(dataTime), separator=',')
        f.write(time_data_str + "\n")

name_Job = 'Job-beamsansmod'
# name_Job = 'Job_1'
# name_Job = "Job-1_5"
# name_Job = "Job-1_5mod"
name_Job = "Job-0_5"
# name_Job = "Job-0_5beam"
name_Job = 'BCC_R1_modMid'
name_step = 'Step-1'
dataRF, dataU, dataTime = get_result(name_Job,name_step)
save_result(0, 6, 1, 1, dataRF, dataU, dataTime)