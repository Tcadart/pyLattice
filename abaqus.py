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
import time
import numpy as np
from Lattice import *


def CreateNodes(name_model, name_Part, node_data):
    p = mdb.models[name_model].parts[name_Part]
    for i in range(len(node_data)):
        p.DatumPointByCoordinate(coords=(node_data[i][1], node_data[i][2], node_data[i][3]))


def CreateBeams(name_model, name_Part, Beam_data):
    p = mdb.models[name_model].parts[name_Part]
    d2 = p.datums
    for i in range(len(Beam_data)):
        p.WirePolyLine(points=((d2[int(Beam_data[i][1]) + 2], d2[int(Beam_data[i][2]) + 2]),), mergeType=IMPRINT,
                       meshable=ON)


lattice = Lattice(1, 1, 1, 2, 2, 2)

# lattice.generate_custom_lattice(7)
cell = lattice.generate_random_lattice(20, 40, 0.01)[0]

# node_data = affichage_points_console(lattice)
node_data = lattice.affichage_points_console()
print("node_data", node_data)

# Beam_data = affichage_beams_console(lattice)
Beam_data = lattice.affichage_beams_console()
print("Beam_data", Beam_data)
# angles = lattice.Getangle()
angles = lattice.angles
print("les angles", angles)
name_model = 'Lattice_cube'
name_Part = 'Lattice_Part'
name_Job = 'Job_1'


def CreateModel(name_model):
    # Create Model
    mdb.Model(name=name_model, modelType=STANDARD_EXPLICIT)


def CreatePart(name_model, name_Part):
    # Create Part
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
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[name_model].sketches['__profile__']
    # Delete Part
    del p.features['Wire-1']


time.sleep(2)
# print("node_data", node_data)
# affichecell()
CreateModel(name_model)
CreatePart(name_model, name_Part)
CreateNodes(name_model, name_Part, node_data)
CreateBeams(name_model, name_Part, Beam_data)
