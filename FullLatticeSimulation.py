# Version python 3.9.7
# ----------------------------------------------------------------------------
# Created By  : Thomas Cadart
# Created Date: 30/06/2023
# version ='1.0'
# ---------------------------------------------------------------------------
""" Code for constructing lattice""" 
# ---------------------------------------------------------------------------


# *******************************************************************************************************************
# *******************************************************************************************************************

#                                           Inputs

# *******************************************************************************************************************
# *******************************************************************************************************************
from abaqus import *
from abaqusConstants import *
import regionToolset
import __main__
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add
import sys
import os

sys.path.insert(8, r"D:/travail_Abaqus/MicroMechanics_v1.18/MicroMechanics")
import microMechanics
from microMechanics.mmpBackend import Interface
from microMechanics.mmpBackend.mmpInterface.mmpRVEConstants import *
from microMechanics.mmpBackend.mmpKernel.mmpLibrary import *
from microMechanics.mmpBackend import mmpKernel as Kernel
from odbAccess import openOdb

import numpy as np
import time
from Lattice import *
import re
from Materials import *
# *******************************************************************************************************************
# *******************************************************************************************************************

#                                           Functions

# *******************************************************************************************************************
# *******************************************************************************************************************

def CreateModel(name_model):
    #Create Model
    mdb.Model(name=name_model, modelType=STANDARD_EXPLICIT)

def CreatePart(name_model,name_Part):
    #Create Part
    s1 = mdb.models[name_model].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Line(point1=(0.0, 0.0), point2=(0.0, 5.0))
    s1.VerticalConstraint(entity=g[2], addUndoState=False)
    p = mdb.models[name_model].Part(name=name_Part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[name_model].parts[name_Part]
    p.BaseWire(sketch=s1)
    s1.unsetPrimaryObject()
    p = mdb.models[name_model].parts[name_Part]
    del mdb.models[name_model].sketches['__profile__']
    # Delete Part
    del p.features['Wire-1']

def CreateNodes(name_model,name_Part):
    p = mdb.models[name_model].parts[name_Part]
    for i in range(len(lattice.nodes)):
        p.DatumPointByCoordinate(coords=(lattice.nodes[i][1], lattice.nodes[i][2], lattice.nodes[i][3])) 

def CreateBeams(name_model,name_Part):
    p = mdb.models[name_model].parts[name_Part]
    d2 = p.datums
    idx_default = 2
    beam_points = []
    for i in range(len(lattice.beams)):
        idx1 = 3 if int(lattice.beams[i][1]) + 2 > 999 else idx_default
        idx2 = 3 if int(lattice.beams[i][2]) + 2 > 999 else idx_default

        datum1 = d2[int(lattice.beams[i][1]) + idx1]
        datum2 = d2[int(lattice.beams[i][2]) + idx2]

        # Ajoutez les points a la liste
        beam_points.append((datum1, datum2))

    # Utilisez la liste des points pour creer les poutres en une seule fois
    p.WirePolyLine(points=beam_points, mergeType=IMPRINT, meshable=ON)

def SetAbaqusWindows(name_model,name_Part):
    p = mdb.models[name_model].parts[name_Part]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)

def Assembly_beam(name_model,name_Part,name_Assembly):
    a1 = mdb.models[name_model].rootAssembly
    a1.DatumCsysByDefault(CARTESIAN)
    p = mdb.models[name_model].parts[name_Part]
    a1.Instance(name=name_Assembly, part=p, dependent=ON)

def Assembly_beam_Surface(name_model,name_Part,name_surface):
    # Assembly 2 parts and move surface at good position
    a1 = mdb.models[name_model].rootAssembly
    p = mdb.models[name_model].parts[name_Part]
    name_Part_Assembly = name_Part+'-1'
    a1.Instance(name=name_Part_Assembly, part=p, dependent=ON)
    p = mdb.models[name_model].parts[name_surface]
    name_surface_Assembly = name_surface+'-1'
    a1.Instance(name=name_surface_Assembly, part=p, dependent=ON)
    a1.translate(instanceList=(name_surface_Assembly, ), vector=(0.0, 0.0, lattice.zMax+0.1*lattice.zMax))
    return name_Part_Assembly, name_surface_Assembly

def CreateSet(name_model,lattice,name_Set,Position,node_data,name_Assembly):
    #Traitement Position
    if 'x' in Position:
        Direction = 1
        if '+' in Position:
            valselected =lattice.xMax
        else:
            valselected = lattice.xMin
    elif 'y' in Position:
        Direction = 2
        if '+' in Position:
            valselected =lattice.yMax
        else:
            valselected = lattice.yMin
    elif 'z' in Position:
        Direction = 3
        if '+' in Position:
            valselected = lattice.zMax
        else:
            valselected = lattice.zMin
    m = mdb.models[name_model]
    a = m.rootAssembly
    v1 = a.instances[name_Assembly].vertices
    verts1 = []
    for i in range(len(node_data)):
        if node_data[i][Direction]==valselected:
            verts1.append(v1.findAt(((node_data[i][1],node_data[i][2],node_data[i][3]),),))
    a.Set(vertices=verts1, name=name_Set)

def mesh_Part(name_model,name_Part,mesh_size):
    p = mdb.models[name_model].parts[name_Part]
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

def getUniqueRadius():
    uniqueRadius = list(set(lattice.radius))
    return uniqueRadius

def Create_Beam_Profile_Mod(name_model,name_Part,VectorOrientation,name_region,name_Beam_Profile):
    #Beam Section Orientation
    p = mdb.models[name_model].parts[name_Part]
    Region_mid,Region_Ext = selectBeamRegion(name_region)
    p.assignBeamSectionOrientation(region=Region_mid, method=N1_COSINES, n1=(VectorOrientation[0], VectorOrientation[1], VectorOrientation[2]))
    p.assignBeamSectionOrientation(region=Region_Ext, method=N1_COSINES, n1=(VectorOrientation[0], VectorOrientation[1], VectorOrientation[2]))
    # Beam Profile creation
    uniqueRadius = getUniqueRadius()
    name = []
    for i in range(len(uniqueRadius)):
        mdb.models[name_model].CircularProfile(name=name_Beam_Profile+'_'+str(uniqueRadius[i]), r=uniqueRadius[i])
        name.append(name_Beam_Profile+'_'+str(uniqueRadius[i]))
    return name, Region_mid, Region_Ext

def Create_Beam_Profile(name_model,name_Part,VectorOrientation,name_region,name_Beam_Profile):
    #Beam Section Orientation
    p = mdb.models[name_model].parts[name_Part]
    AllBeams = selectBeamRegion(name_region)
    p.assignBeamSectionOrientation(region=AllBeams, method=N1_COSINES, n1=(VectorOrientation[0], VectorOrientation[1], VectorOrientation[2]))
    # Beam Profile creation
    uniqueRadius = getUniqueRadius()
    name = []
    for i in range(len(uniqueRadius)):
        mdb.models[name_model].CircularProfile(name=name_Beam_Profile+'_'+str(uniqueRadius[i]), r=uniqueRadius[i])
        name.append(name_Beam_Profile+'_'+str(uniqueRadius[i]))
    return name, AllBeams

def selectBeamRegion(name_region):
    p = mdb.models[name_model].parts[name_Part]
    e = p.edges
    if name_region == 'AllBeams': #Traitement all beams
        edges = e.getByBoundingBox(xMin=lattice.xMin,yMin=lattice.yMin,zMin=lattice.zMin,xMax=lattice.xMax,yMax=lattice.yMax,zMax=lattice.zMax)
        Region = p.Set(edges=edges, name=name_region)
    elif name_region == 'BeamMod':
        edges = e.getByBoundingBox(xMin=lattice.xMin,yMin=lattice.yMin,zMin=lattice.zMin,xMax=lattice.xMax,yMax=lattice.yMax,zMax=lattice.zMax)
        Region = p.Set(edges=edges, name='AllBeams')
        # For normal beams
        edges_mid = []
        for i in range(len(lattice.beams)):
            if lattice.beams[i][3] == 0: #Check type of beam Center or modified
                edges_mid.append(e.findAt(((lattice.nodes[lattice.beams[i][1]][1]+(lattice.nodes[lattice.beams[i][2]][1]-lattice.nodes[lattice.beams[i][1]][1])/2,lattice.nodes[lattice.beams[i][1]][2]+(lattice.nodes[lattice.beams[i][2]][2]-lattice.nodes[lattice.beams[i][1]][2])/2,lattice.nodes[lattice.beams[i][1]][3]+(lattice.nodes[lattice.beams[i][2]][3]-lattice.nodes[lattice.beams[i][1]][3])/2),),))
        region_mid = p.Set(edges=edges_mid, name='BeamMid')
        # For Modified beams
        p = mdb.models[name_model].parts[name_Part]
        region_ext = p.SetByBoolean(name = 'BeamMod',operation=DIFFERENCE,sets=(p.sets['AllBeams'],p.sets['BeamMid'],))
        return region_mid, region_ext
    else: # Traitement couche par couche ################## Probleme lors de poutres qui apparraisent dans 2 couche differentes Il faut un traitement different
        couche_number = int(re.search(r'\d+', name_region).group())
        if 'X' in name_region:
            edges = e.getByBoundingBox(xMin=lattice.xMax-((lattice.num_cells_x+1-couche_number)*lattice.cell_size_x),yMin=lattice.yMin,zMin=lattice.zMin,xMax=lattice.xMax-((lattice.num_cells_x-couche_number)*lattice.cell_size_x),yMax=lattice.yMax,zMax=lattice.zMax)
            Region = p.Set(edges=edges, name=name_region)
        if 'Y' in name_region:
            edges = e.getByBoundingBox(xMin=lattice.xMin,yMin=lattice.yMax-((lattice.num_cells_y+1-couche_number)*lattice.cell_size_y),zMin=lattice.zMin,xMax=lattice.xMax,yMax=lattice.yMax-((lattice.num_cells_y-couche_number)*lattice.cell_size_y),zMax=lattice.zMax)        
            Region = p.Set(edges=edges, name=name_region)
        if 'Z' in name_region:
            edges = e.getByBoundingBox(xMin=lattice.xMin,yMin=lattice.yMin,zMin=lattice.zMax-((lattice.num_cells_z+1-couche_number)*lattice.cell_size_z),xMax=lattice.xMax,yMax=lattice.yMax,zMax=lattice.zMax-((lattice.num_cells_z-couche_number)*lattice.cell_size_z))
            Region = p.Set(edges=edges, name=name_region)
    return Region

def create_material(name_model,Material_Type,Material_Carac):
    m = MatProperties(Material_Type,Material_Carac)
    mdb.models[name_model].Material(name=m.name)
    mdb.models[name_model].materials[m.name].Density(table=((m.density, ), ))
    mdb.models[name_model].materials[m.name].Elastic(table=((m.elastic[0], m.elastic[1]), ))
    if Material_Carac == 1 and Material_Type == 0:
        mdb.models[name_model].materials[m.name].Plastic(scaleStress=None, table=m.plastic)
    elif Material_Carac == 1 and Material_Type == 3:
        mdb.models[name_model].materials[m.name].Plastic(hardening=JOHNSON_COOK, 
            scaleStress=None, table=(m.plastic, ))
    return m.name

def Create_Beam_Section(name_model,name_Section,name_material,name_Beam_Profile,Region):
    # Beam Section creation
    mdb.models[name_model].BeamSection(name=name_Section, integration=DURING_ANALYSIS, 
    poissonRatio=0.0, profile=name_Beam_Profile, material=name_material, 
    temperatureVar=LINEAR, consistentMassMatrix=False)
    # Beam assignment
    p = mdb.models[name_model].parts[name_Part]
    p.SectionAssignment(region=Region, sectionName=name_Section, offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

def Create_Job(name_Job,name_model):
    #Create a job
    mdb.Job(name=name_Job, model=name_model, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

def Submit_Job(name_Job,name_model):
    #Submit and wait
    job1 = mdb.Job(name=name_Job,model=name_model)
    job1.submit()
    job1.waitForCompletion()

def Create_Step(name_model,name_step,name_previous_step,AnalysisType):
    if AnalysisType == 1:
        # mdb.models[name_model].StaticStep(name=name_step, previous=name_previous_step)
        mdb.models[name_model].StaticStep(name=name_step, previous=name_previous_step, nlgeom=ON,
                                                stabilizationMagnitude=0.0002,
                                                stabilizationMethod=DISSIPATED_ENERGY_FRACTION,
                                                continueDampingFactors=False, adaptiveDampingRatio=0.05, initialInc=0.1)
    elif AnalysisType == 2:
        mdb.models[name_model].StaticStep(name=name_step, previous=name_previous_step, nlgeom=ON)

def Create_BC_Fixed(name_model,name_step,name_set):
    a = mdb.models[name_model].rootAssembly
    region = region = a.sets[name_set]
    mdb.models[name_model].EncastreBC(name='Fixed', createStepName=name_step, region=region, localCsys=None)

def Create_Loads(name_model,name_step,name_set,name_load,vector_load):
    a = mdb.models[name_model].rootAssembly
    region = a.sets[name_set]
    mdb.models[name_model].DisplacementBC(name=name_load, createStepName=name_step, 
    region=region, u1=vector_load[0], u2=vector_load[1], u3=vector_load[2], ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

def Create_Loads_Surface(name_model,name_step,name_surface_Assembly,name_load,vector_load):
    a = mdb.models[name_model].rootAssembly
    r1 = a.instances[name_surface_Assembly].referencePoints
    refPoints1=(r1[2], )
    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models[name_model].DisplacementBC(name=name_load, createStepName=name_step, 
    region=region, u1=vector_load[0], u2=vector_load[1], u3=vector_load[2], ur1=UNSET, ur2=UNSET, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


def Create_GetReactionForce(name_model,name_step,name_set):
    regionDef=mdb.models[name_model].rootAssembly.sets[name_set]
    mdb.models[name_model].HistoryOutputRequest(name='Reaction_Force', 
        createStepName=name_step, variables=('RF1', 'RF2', 'RF3'), region=regionDef, 
        sectionPoints=DEFAULT, rebar=EXCLUDE)

def Create_GetDisplacement(name_model,name_step,name_set):
    regionDef=mdb.models[name_model].rootAssembly.sets[name_set]
    mdb.models[name_model].HistoryOutputRequest(name='Displacement', 
        createStepName=name_step, variables=('U1', 'U2', 'U3'), region=regionDef, 
        sectionPoints=DEFAULT, rebar=EXCLUDE)
    
def Delete_output_default(name_model):
    del mdb.models[name_model].historyOutputRequests['H-Output-1']
    mdb.models[name_model].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'U'))


def getVectorOrientation(Lattice):
    if (Lattice == 0):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 1):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 2):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 3):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 4):
        VectorOrientation = [0,1,-1]
    elif (Lattice == 5):
        VectorOrientation = [0,1,-1]
    elif (Lattice == 6):
        VectorOrientation = [1,-1,1]
    elif (Lattice == 7):
        VectorOrientation = [1,-1,-1]
    elif (Lattice == 8):
        VectorOrientation = [0,0,-1]
    elif (Lattice == 9):
        VectorOrientation = [1,0,-1]
    elif (Lattice == 10):
        VectorOrientation = [1,0,-1]
    elif (Lattice == 11):
        VectorOrientation = [1,0,-1]
    else:
        VectorOrientation = [0,0,-1]
    return VectorOrientation

def visualizationSimulation(name_Job):
    o3 = session.openOdb(name='D:/travail_Abaqus/Lattice/'+name_Job+'.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=1 )


def constructLatticeAbaqus(name_model,name_Part,name_Assembly,Lattice_Type,MethodSim):
    CreateModel(name_model)
    CreatePart(name_model,name_Part)
    CreateNodes(name_model,name_Part)
    CreateBeams(name_model,name_Part)
    mesh_Part(name_model,name_Part,0.087)
    if MethodSim == 1:
        name_region = 'BeamMod'
        name_Beam_Profile = 'Circ'
        name_Beam_Profile,Region_mid,Region_Ext = Create_Beam_Profile_Mod(name_model,name_Part,getVectorOrientation(Lattice_Type),name_region,name_Beam_Profile)
    elif MethodSim == 0:
        name_region = 'AllBeams'
        name_Beam_Profile = 'Circ'
        name_Beam_Profile, AllBeam = Create_Beam_Profile(name_model,name_Part,getVectorOrientation(Lattice_Type),name_region,name_Beam_Profile)
    
    name_material = create_material(name_model,4,1)

    if MethodSim == 1:
        Create_Beam_Section(name_model,'CircBeam',name_material,name_Beam_Profile[1],Region_mid)
        Create_Beam_Section(name_model,'CircBeammod',name_material,name_Beam_Profile[0],Region_Ext)
    elif MethodSim == 0:
        Create_Beam_Section(name_model,'CircBeam',name_material,name_Beam_Profile[0],AllBeam)

# def constructLatticeAbaqusMultiMat(name_model,name_Part,name_Assembly,node_data,Beam_data,Radius,Lattice_Type):
#     CreateModel(name_model)
#     CreatePart(name_model,name_Part)
#     CreateNodes(name_model,name_Part)
#     CreateBeams(name_model,name_Part)
#     mesh_beam(name_model,name_Part)
#     name_region = ['Couche _Z1','Couche _Z2','Couche _Z3','Couche _Z4']
#     Material_type = [0,1,2,0]
#     Radius = [0.5,0.055,0.2,0.01]
#     for i in range(lattice.num_cells_z):
#         name_Beam_Profile = 'Circ'+str(i)
#         #region = Create_Beam_Profile(name_model,name_Part,getVectorOrientation(Lattice_Type),name_region[i],name_Beam_Profile)
#         name_material = create_material(name_model,Material_type[i],0)
#         Create_Beam_Section(name_model,'CircBeam_'+str(i),name_material,name_Beam_Profile,region)
#     Assembly_beam(name_model,name_Part,name_Assembly)

def getReactionForce(name_Job):
    # Create file with all reaction force data
    file_path = 'D:/travail_Abaqus/Lattice/abaqus.rpt'
    odb = session.odbs['D:/travail_Abaqus/Lattice/'+name_Job+'.odb']
    session.writeFieldReport(fileName=file_path, append=OFF, 
    sortItem='Node Label', odb=odb, step=0, frame=1, outputPosition=NODAL, 
    variable=(('RF', NODAL), ), stepFrame=SPECIFY)
    reactionforceData = []
    with open(file_path, 'r') as file:
        is_data_section = False  # Indicateur pour savoir si nous sommes dans la section des donnees
        for line in file:
            # Recherche de la ligne qui indique le debut des donnees des noeuds
            if "Node Label" in line:
                is_data_section = True
                continue
            if "Minimum" in line:
                is_data_section = False
                continue
            if is_data_section:
                # Separation des valeurs en utilisant les espaces comme separateurs
                values = line.split()
                if len(values) == 5 and float(values[1]) != 0:  # Nous nous attendons a 4 valeurs numeriques par ligne
                    node_label = int(values[0])
                    rf_magnitude = float(values[1])
                    rf_rf1 = float(values[2])
                    rf_rf2 = float(values[3])
                    rf_rf3 = float(values[4])
                    reactionforceData.append((node_label, rf_magnitude, rf_rf1, rf_rf2, rf_rf3))
    fixedNodeReactionForce = []
    reactionForceZ = 0
    for i in range(len(reactionforceData)):
        if reactionforceData[i][4] > 0:
            fixedNodeReactionForce.append(reactionforceData[i])
            reactionForceZ += reactionforceData[i][4]
    print(reactionForceZ)
    return reactionForceZ

def get_result(name_Job,name_step):
    odb = openOdb(name_Job+'.odb')
    NameOutput = odb.steps[name_step].historyRegions
    dataRF3 = []
    dataU3 = []
    dataTime = []
    for i in range(len(NameOutput)):
        try:
            if len(dataRF3) == 0:
                rf3_data = [odb.steps[name_step].historyRegions[NameOutput.keys()[i]].historyOutputs['RF3'].data][0]
                for i in range(len(rf3_data)):
                    dataRF3.append(rf3_data[i][1])
            else:
                rf3_data = [odb.steps[name_step].historyRegions[NameOutput.keys()[i]].historyOutputs['RF3'].data][0]
                for i in range(len(rf3_data)):
                    dataRF3[i] = dataRF3[i] + rf3_data[i][1]
        except:
            try:
                if len(dataU3) == 0:
                    U3_data = [odb.steps[name_step].historyRegions[NameOutput.keys()[i]].historyOutputs['U3'].data][0]
                    time_data = [odb.steps[name_step].historyRegions[NameOutput.keys()[i]].historyOutputs['U3'].data][0]
                    for i in range(len(U3_data)):
                        dataU3.append(U3_data[i][1])
                        dataTime.append(time_data[i][0])
            except:
                continue
    return dataRF3, dataU3, dataTime

def save_result(Lattice_Type, number_cell, AnalysisType, MethodSim, dataRF3, dataU3, dataTime):
    f = open(Type_lattice(Lattice_Type)+"_"+str(number_cell)+str(AnalysisType)+str(MethodSim)+".txt", "w")
    f.writelines([np.array2string(np.array(dataRF3), separator=',')+"\n"])
    f.writelines([np.array2string(np.array(dataU3), separator=',')+"\n"])
    f.writelines([np.array2string(np.array(dataTime), separator=',')+"\n"])
    f.close()

def CreateSurface(name_model,name_set_RF):
    # Create Rectangle Surface
    namePart = 'Surface-Top'
    s = mdb.models[name_model].ConstrainedSketch(name='__profile__',sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)
    xPoint1 = lattice.xMin-(0.1*lattice.xMax)
    yPoint1 = lattice.yMin-(0.1*lattice.yMax)
    xPoint2 = lattice.xMax+(0.1*lattice.xMax)
    yPoint2 = lattice.yMax+(0.1*lattice.yMax)
    s.rectangle(point1=(xPoint1, yPoint1), point2=(xPoint2, yPoint2))
    p = mdb.models[name_model].Part(name=namePart, dimensionality=THREE_D,type=DISCRETE_RIGID_SURFACE)
    p = mdb.models[name_model].parts[namePart]
    p.BaseShell(sketch=s)
    # Create Reference Point
    p = mdb.models[name_model].parts[namePart]
    v, e, d, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=p.InterestingPoint(edge=e[3], rule=MIDDLE))
    r = p.referencePoints
    refPoints=(r[2], )
    p.Set(referencePoints=refPoints, name=name_set_RF)
    # Assign section
    # name_section_surface = 'Shell_rigid'
    # numberIntegrationPoint = 5
    # Epaisseur = 1.0
    # mdb.models[name_model].HomogeneousShellSection(name=name_section_surface, 
    # preIntegrate=OFF, material=name_material, thicknessType=UNIFORM, 
    # thickness=Epaisseur, thicknessField='', nodalThicknessField='', 
    # idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
    # thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
    # integrationRule=SIMPSON, numIntPts=numberIntegrationPoint)
    return namePart

def createContactProperty(name_model,name_contact_prop,name_contact,name_step):
    mdb.models[name_model].ContactProperty(name_contact_prop)
    mdb.models[name_model].interactionProperties[name_contact_prop].TangentialBehavior(
        formulation=FRICTIONLESS)
    mdb.models[name_model].interactionProperties[name_contact_prop].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    #: The interaction property "Contact_prop" has been created.
    mdb.models[name_model].ContactStd(name=name_contact, 
        createStepName=name_step)
    mdb.models[name_model].interactions[name_contact].includedPairs.setValuesInStep(
        stepName=name_step, useAllstar=ON)
    mdb.models[name_model].interactions[name_contact].contactPropertyAssignments.appendInStep(
        stepName=name_step, assignments=((GLOBAL, SELF, name_contact_prop), ))
    #: The interaction "Surface to beam" has been created.

def gradSettings(rule,direction,parameters,number_cell_X,number_cell_Y, number_cell_Z):
    """
    Generate gradient settings based on the provided rule, direction, and parameters.

    :param rule: Gradient rule ('constant', 'linear', 'parabolic', 'sinusoide', 'exponential')
    :param direction: Direction of gradient for each axis (list of 3 integers [X,Y,Z])
    :param parameters: Gradient parameters for each axis (list of 3 floats [X,Y,Z])
    :param number_cell_X: Number of cells along X axis (int)
    :param number_cell_Y: Number of cells along Y axis (int)
    :param number_cell_Z: Number of cells along Z axis (int)
    :return: Generated gradient settings (list of lists)
    """
    def apply_rule(i, number_cell, dir_value, param_value, rule):
        if i < number_cell:
            if i >= 1 and dir_value == 1:
                if rule == 'constant':
                    return 1.0
                elif rule == 'linear':
                    return i * param_value
                elif rule == 'parabolic':
                    return i * param_value if i < number_cell / 2 else (number_cell - i - 1) * param_value
                elif rule == 'sinusoide':
                    if i < number_cell / 4:
                        return i * parameters
                    elif i < number_cell / 2:
                        return (number_cell / 2 - i) * parameters
                    elif i == number_cell / 2:
                        return 1.0
                    elif i < 3 / 4 * number_cell:
                        return (3 / 4 * number_cell - i) * (1 / parameters)
                elif rule == 'exponential':
                    return math.exp(i * param_value)
            return 1.0
        return 0.0  

    # Initialization matrix
    max_cells = max(number_cell_X, number_cell_Y, number_cell_Z)
    gradient = [[0.0, 0.0, 0.0] for _ in range(max_cells)]

    # Processing multiple rules
    for i in range(max_cells):
        number_cells = [number_cell_X, number_cell_Y, number_cell_Z]
        for dim_index in range(3):
            gradient[i][dim_index] = apply_rule(i, number_cells[dim_index], direction[dim_index], parameters[dim_index], rule)
    return gradient

def gradMaterialSetting(Multimat,direction,number_cell_X,number_cell_Y, number_cell_Z):
    if Multimat == 0: # Mono material
        gradMat = [[[1 for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    elif Multimat == 1: # Graded material
        if direction == 1:
            gradMat = [[[X for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
        if direction == 2:
            gradMat = [[[Y for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
        if direction == 3:
            gradMat = [[[Z for X in range(number_cell_X)] for Y in range(number_cell_Y)] for Z in range(number_cell_Z)]
    return gradMat

def delete_all_models():
    model_names = list(mdb.models.keys())
    for name in model_names:
        if name != 'Model-1':  # 'Model-1' basic model
            del mdb.models[name]


#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Fonction for

#*******************************************************************************************************************
#*******************************************************************************************************************

def Lattice_geometry_solid(Lattice,Radius_geom): #Tricks -0.0001 and 1.0001 for trigger cutting part outside cell
    BCC = [(-0.0001, -0.0001, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (0.5, 0.5, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (0.5, 0.5, 0.5, -0.0001, -0.0001, 1.0001, Radius_geom),
        (0.5, 0.5, 0.5, -0.0001, 1.0001, -0.0001, Radius_geom),
        (0.5, 0.5, 0.5, -0.0001, 1.0001, 1.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0001, -0.0001, -0.0001, Radius_geom)]
    Octet = [(-0.0001, -0.0001, -0.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (-0.0001, -0.0001, 1.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (1.0001, -0.0001, -0.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (-0.0001, -0.0001, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, 1.0001, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, -0.0001, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, 1.0001, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, -0.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (1.0001, 1.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (1.0001, -0.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (-0.0001, 1.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (-0.0001, -0.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 1.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (-0.0001, 1.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (1.0001, -0.0001, -0.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (1.0001, 0.5, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (-0.0001, 1.0001, -0.0001, 0.5, 1.0001, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (-0.0001, 1.0001, 1.0001, 0.5, 1.0001, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, -0.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 1.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 0.5, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, 1.0001, 0.5, 0.5, Radius_geom)]
    OctetExt = [(-0.0001, -0.0001, -0.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (-0.0001, -0.0001, 1.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (1.0001, -0.0001, -0.0001, 0.5, -0.0001, 0.5, Radius_geom),
        (-0.0001, -0.0001, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, 1.0001, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, -0.0001, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, 1.0001, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (-0.0001, -0.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (1.0001, 1.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (1.0001, -0.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (-0.0001, 1.0001, -0.0001, 0.5, 0.5, -0.0001, Radius_geom),
        (-0.0001, -0.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 1.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (-0.0001, 1.0001, 1.0001, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (1.0001, -0.0001, -0.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (1.0001, 0.5, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (-0.0001, 1.0001, -0.0001, 0.5, 1.0001, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (-0.0001, 1.0001, 1.0001, 0.5, 1.0001, 0.5, Radius_geom)]
    OctetInt = [(0.5, -0.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, -0.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 1.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 0.5, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, 1.0001, 0.5, 0.5, Radius_geom)]
    BCCZ = [(0.5, 0.5, 0.5, 1.0001, 1.0001, 1.0001, Radius_geom),
        (-0.0001, -0.0001, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 1.0001, 1.0001, -0.0001, Radius_geom),
        (-0.0001, -0.0001, 1.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, -0.0001, 1.0001, -0.0001, Radius_geom),
        (1.0001, -0.0001, 1.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, -0.0001, 1.0001, 1.0001, Radius_geom),
        (1.0001, -0.0001, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 0.5, 0.5, 0.5, 1.0001, Radius_geom)]
    Cubic = [(-0.0001, -0.0001, -0.0001, -0.0001, -0.0001, 1.0001, Radius_geom),
          (1.0001, -0.0001, -0.0001, 1.0001, -0.0001, 1.0001, Radius_geom),
          (-0.0001, 1.0001, -0.0001, -0.0001, 1.0001, 1.0001, Radius_geom),
          (1.0001, 1.0001, -0.0001, 1.0001, 1.0001, 1.0001, Radius_geom),
          (-0.0001, -0.0001, -0.0001, 1.0001, -0.0001, -0.0001, Radius_geom),
          (-0.0001, -0.0001, -0.0001, -0.0001, 1.0001, -0.0001, Radius_geom),
          (1.0001, 1.0001, -0.0001, -0.0001, 1.0001, -0.0001, Radius_geom),
          (1.0001, 1.0001, -0.0001, 1.0001, -0.0001, -0.0001, Radius_geom),
          (-0.0001, -0.0001, 1.0001, 1.0001, -0.0001, 1.0001, Radius_geom),
          (-0.0001, -0.0001, 1.0001, -0.0001, 1.0001, 1.0001, Radius_geom),
          (1.0001, 1.0001, 1.0001, -0.0001, 1.0001, 1.0001, Radius_geom),
          (1.0001, 1.0001, 1.0001, 1.0001, -0.0001, 1.0001, Radius_geom)]
    OctahedronZ = [(0.5, -0.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, -0.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 1.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 0.5, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, -0.0001, 0.5, 0.5, 1.0001, Radius_geom)]
    OctahedronZcross = [(0.5, -0.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, -0.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 1.0001, 0.5, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, 1.0001, Radius_geom),
        (1.0001, 0.5, 0.5, 0.5, 0.5, -0.0001, Radius_geom),
        (0.5, 0.5, -0.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, -0.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, 1.0001, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 0.5, 1.0001, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, -0.0001, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (-0.0001, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (1.0001, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
        (0.5, 1.0001, 0.5, 0.5, 0.5, 0.5, Radius_geom)]
    Kelvin = [(0.5, 0.25, -0.0001, 0.25, 0.5, -0.0001, Radius_geom),
              (0.5, 0.25, -0.0001, 0.75, 0.5, -0.0001, Radius_geom),
              (0.5, 0.75, -0.0001, 0.25, 0.5, -0.0001, Radius_geom),
              (0.5, 0.75, -0.0001, 0.75, 0.5, -0.0001, Radius_geom),
              (0.5, 0.25, 1.0001, 0.25, 0.5, 1.0001, Radius_geom),
              (0.5, 0.25, 1.0001, 0.75, 0.5, 1.0001, Radius_geom),
              (0.5, 0.75, 1.0001, 0.25, 0.5, 1.0001, Radius_geom),
              (0.5, 0.75, 1.0001, 0.75, 0.5, 1.0001, Radius_geom),
              (0.5, -0.0001, 0.25, 0.25, -0.0001, 0.5, Radius_geom),
              (0.5, -0.0001, 0.25, 0.75, -0.0001, 0.5, Radius_geom),
              (0.5, -0.0001, 0.75, 0.25, -0.0001, 0.5, Radius_geom),
              (0.5, -0.0001, 0.75, 0.75, -0.0001, 0.5, Radius_geom),
              (0.5, 1.0001, 0.25, 0.25, 1.0001, 0.5, Radius_geom),
              (0.5, 1.0001, 0.25, 0.75, 1.0001, 0.5, Radius_geom),
              (0.5, 1.0001, 0.75, 0.25, 1.0001, 0.5, Radius_geom),
              (0.5, 1.0001, 0.75, 0.75, 1.0001, 0.5, Radius_geom),
              (-0.0001, 0.5, 0.25, -0.0001, 0.25, 0.5, Radius_geom),
              (-0.0001, 0.5, 0.25, -0.0001, 0.75, 0.5, Radius_geom),
              (-0.0001, 0.5, 0.75, -0.0001, 0.25, 0.5, Radius_geom),
              (-0.0001, 0.5, 0.75, -0.0001, 0.75, 0.5, Radius_geom),
              (1, 0.5, 0.25, 1.0001, 0.25, 0.5, Radius_geom),
              (1, 0.5, 0.25, 1.0001, 0.75, 0.5, Radius_geom),
              (1, 0.5, 0.75, 1.0001, 0.25, 0.5, Radius_geom),
              (1, 0.5, 0.75, 1.0001, 0.75, 0.5, Radius_geom),
              (0.5, 0.25, -0.0001, 0.5, -0.0001, 0.25, Radius_geom),
              (0.25, 0.5, -0.0001, -0.0001, 0.5, 0.25, Radius_geom),
              (0.75, 0.5, -0.0001, 1.0001, 0.5, 0.25, Radius_geom),
              (0.5, 0.75, -0.0001, 0.5, 1.0001, 0.25, Radius_geom),
              (0.25, -0.0001, 0.5, -0.0001, 0.25, 0.5, Radius_geom),
              (0.75, -0.0001, 0.5, 1.0001, 0.25, 0.5, Radius_geom),
              (0.75, 1.0001, 0.5, 1.0001, 0.75, 0.5, Radius_geom),
              (0.25, 1.0001, 0.5, -0.0001, 0.75, 0.5, Radius_geom),
              (0.5, -0.0001, 0.75, 0.5, 0.25, 1.0001, Radius_geom),
              (-0.0001, 0.5, 0.75, 0.25, 0.5, 1.0001, Radius_geom),
              (0.5, 1.0001, 0.75, 0.5, 0.75, 1.0001, Radius_geom),
              (1, 0.5, 0.75, 0.75, 0.5, 1.0001, Radius_geom)]
    CubicV2 = [(0.5, -0.0001, 0.5, 0.5, 0.5, 0.5, Radius_geom), #Tricks -0.0001 and 1.0001 for trigger cutting part outside cell
               (-0.0001, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 1.0001, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (1.0001, 0.5, 0.5, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, -0.0001, 0.5, 0.5, 0.5, Radius_geom),
               (0.5, 0.5, 1.0001, 0.5, 0.5, 0.5, Radius_geom)]
    if (Lattice == 0):
        return BCC
    if (Lattice == 1):
        return Octet
    if (Lattice == 2):
        return OctetExt
    if (Lattice == 3):
        return OctetInt
    if (Lattice == 4):
        return BCCZ
    if (Lattice == 5):
        return Cubic
    if (Lattice == 6):
        return OctahedronZ
    if (Lattice == 7):
        return OctahedronZcross
    if (Lattice == 8):
        return Kelvin
    if (Lattice == 9):
        return CubicV2

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
        Type = 'OctahedronZ'
    if (Lattice == 7):
        Type = 'OctahedronZcross'
    if (Lattice == 8):
        Type = 'Kelvin'
    if (Lattice == 9):
        Type = 'CubicV2'
    if (Lattice == 10):
        Type = 'CubicV3'
    if (Lattice == 11):
        Type = 'CubicV4'
    return Type

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************

name_model = 'BCC_test'
name_Job = 'Job_1'
name_Part = 'Lattice_Part'
name_Assembly = 'Lattice_assembly'
VectorOrientation = [0,0,-1]
Radius = 0.5
cell_size = 7.5
cell_size_X = cell_size
cell_size_Y = cell_size
cell_size_Z = cell_size
number_cell = 5
number_cell_X = number_cell
number_cell_Y = number_cell
number_cell_Z = 8

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

# Gradient on cell dimensions
GradDimRule = 'constant'
GradDimDirection = [1,0,0]
GradDimParameters = [1.5,0.0,2.0] #Float
# Gradient on radius of beams
GradRadRule = 'constant'
GradRadDirection = [0,0,1]
GradRadParameters = [1.0,0.0,2.0]
#Gradient Rule
# - constant
# - linear
# - parabolic
# - sinusoide
# - exponential

Multimat = 0
# 0 -> materiaux
# 1 -> multimat par couche
GradMaterialDirection = 1 # 1:X / 2:Y / 3:Z

AnalysisType = 1
# 0 Modelisation lattice only
# 1 Compression Z static
# 2 Compression Z Implicit
# 3 Compression BCC Solid Plastic Static
compressionPourcent = 20

MethodSim = 1
# 0 No modification
# 1 Node Modification

# delete_all_models()
#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Main

#*******************************************************************************************************************
#*******************************************************************************************************************
# Gradient properties
gradDim = gradSettings(GradDimRule,GradDimDirection,GradDimParameters,number_cell_X,number_cell_Y, number_cell_Z) # sortir des floats
gradRadius = gradSettings(GradRadRule,GradRadDirection,GradRadParameters,number_cell_X,number_cell_Y, number_cell_Z) # sortir des floats
gradMat = gradMaterialSetting(Multimat,GradMaterialDirection,number_cell_X,number_cell_Y, number_cell_Z)


#Generate data from lattice
lattice = Lattice(cell_size_X,cell_size_Y,cell_size_Z, number_cell_X,number_cell_Y,number_cell_Z,Lattice_Type, Radius,gradRadius,gradDim,gradMat,MethodSim)

if AnalysisType != 3:
    #Generate lattice on abaqus
    if Multimat == 0:
        constructLatticeAbaqus(name_model,name_Part,name_Assembly,Lattice_Type,MethodSim)
    # else:
        # constructLatticeAbaqusMultiMat(name_model,name_Part,name_Assembly,lattice.nodes,lattice.beams,Radius,Lattice_Type)


if AnalysisType == 1:
    # Create assembly for lattice 
    Assembly_beam(name_model,name_Part,name_Assembly)
    # Create set for fixed nodes and loaded nodes
    CreateSet(name_model,lattice,'Fixed_nodes',"z-",lattice.nodes,name_Assembly)
    CreateSet(name_model,lattice,'loaded_nodes',"z+",lattice.nodes,name_Assembly)
    # Create loading condition and submit job
    name_step = 'Step-1'
    Create_Step(name_model,name_step,'Initial',AnalysisType)
    Create_BC_Fixed(name_model,'Initial','Fixed_nodes')
    displacementCompression = (lattice.zMax-lattice.zMin)*compressionPourcent/100
    Create_Loads(name_model,name_step,'loaded_nodes','Load_1',[0,0,-displacementCompression])
    Create_GetDisplacement(name_model,name_step,'loaded_nodes')
    Create_GetReactionForce(name_model,name_step,'Fixed_nodes')
    Delete_output_default(name_model)
    Submit_Job(name_Job,name_model)
    visualizationSimulation(name_Job)
    dataRF3, dataU3, dataTime = get_result(name_Job,name_step)
    save_result(Lattice_Type, number_cell, AnalysisType, MethodSim, dataRF3, dataU3, dataTime)
elif AnalysisType == 2:
    # Create superior surface to create compression on lattice
    name_set_RF = 'Set-RF'
    nameSurfacePart = CreateSurface(name_model,name_set_RF)
    mesh_Part(name_model,nameSurfacePart,0.1)
    name_Part_Assembly, name_surface_Assembly = Assembly_beam_Surface(name_model,name_Part,nameSurfacePart)
    # Fix bottom nodes 
    CreateSet(name_model,lattice,'Fixed_nodes',"z-",lattice.nodes,name_Part_Assembly)
    Create_BC_Fixed(name_model,'Initial','Fixed_nodes')
    createContactProperty(name_model,'contact_prop','contact_surface_beam','Initial')
    Create_Step(name_model,'Step-1','Initial',AnalysisType)
    displacementCompression = (lattice.zMax-lattice.zMin)*compressionPourcent/100
    Create_Loads_Surface(name_model,'Step-1',name_surface_Assembly,'Displacement_surface',[0,0,-displacementCompression])
    Submit_Job(name_Job,name_model)
    visualizationSimulation()
elif AnalysisType == 3:
    Lattice_geom = Lattice_geometry_solid(Lattice_Type,Radius)
    Maillage_size = 0.5
    init_name_model = 'BCC_05_'
    init_name_Job = 'Job_BCC_05_'
    for i in range(6):
        name_model = init_name_model+str(i)
        name_Job = init_name_Job+str(i)
        number_cell = i+1
        number_cell_X = number_cell
        number_cell_Y = number_cell
        number_cell_Z = number_cell
        generatedModel = Kernel.Library.Lattice.generateGeneralLattice(name_model,(cell_size_X,cell_size_Y,cell_size_Z),Lattice_geom,(number_cell_X, number_cell_Y, number_cell_Z),Maillage_size)
        # Change material to Johnson Cook
        mdb.models[name_model].materials['Ti-6Al-4V'].plastic.setValues(
            hardening=JOHNSON_COOK, scaleStress=None, table=((460.0, 1450.0, 1.31, 0.85, 1500.0, 1000.0), ))
        # Create Step
        mdb.models[name_model].StaticStep(name='StaticPlastic', previous='Initial')
        #Create Set for BC
        a = mdb.models[name_model].rootAssembly
        instances = a.instances
        for instance_name, instance in instances.items():
            instance_name_save = instance_name
        f1 = a.instances[str(instance_name_save)].faces
        delta = 0.1
        faces = f1.getByBoundingBox(xMin = -delta,yMin =-delta, zMin=cell_size_Z * number_cell_Z -delta, xMax=cell_size_X * number_cell_X +delta, yMax=cell_size_Y * number_cell_Y +delta, zMax=cell_size_Z * number_cell_Z +delta)
        region = a.Set(faces=faces, name='Top_Surface')
        displacement = (number_cell*cell_size_Z)/2
        mdb.models[name_model].DisplacementBC(name='Displacement', 
            createStepName='StaticPlastic', region=region, u1=0.0, u2=0.0, u3=-displacement, 
            ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
            distributionType=UNIFORM, fieldName='', localCsys=None)
        faces = f1.getByBoundingBox(xMin = -delta,yMin =-delta, zMin=-delta, xMax=cell_size_X * number_cell_X +delta, yMax=cell_size_Y * number_cell_Y +delta, zMax=+delta)
        region = a.Set(faces=faces, name='Bottom_region')
        mdb.models[name_model].EncastreBC(name='Fixed', createStepName='StaticPlastic', region=region, localCsys=None)
        # Create Output Request
        regionDef=mdb.models[name_model].rootAssembly.sets['Top_Surface']
        mdb.models[name_model].HistoryOutputRequest(name='Displacement', 
            createStepName='StaticPlastic', variables=('U3', ), region=regionDef, 
            sectionPoints=DEFAULT, rebar=EXCLUDE)
        regionDef=mdb.models[name_model].rootAssembly.sets['Bottom_region']
        mdb.models[name_model].HistoryOutputRequest(name='ReactionForce', 
            createStepName='StaticPlastic', variables=('RF3', ), 
            region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
        mdb.models[name_model].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'U'))
        del mdb.models[name_model].historyOutputRequests['H-Output-1']
        # Create Job
        mdb.Job(name=name_Job, model=name_model, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
            multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    # Submit_Job(name_Job,name_model)
    


# Affichage final
if AnalysisType == 0:
    #Visualization on abaqus
    SetAbaqusWindows(name_model,name_Part)

