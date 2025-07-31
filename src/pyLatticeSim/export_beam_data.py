import math
import os
import numpy as np

from dolfinx import io, fem
import gmsh
import ufl

from .utils import clear_directory

class exportBeamData:
    """
    Beam data exportation to Paraview

    Parameters:
    ------------
    name_file: string
        Name of the file to export data in
    """

    def __init__(self, name_file:str, simulation_model):
        self.simulation_model = simulation_model
        self.name_file = name_file
        self._vtk = io.VTKFile(self.simulation_model.domain.comm, self.name_file, "w")

        self.result_to_export = []
        self.M_function = None
        self.element_mesh = {}


    def define_function_space(self):
        """
        Define the function space for forces and moments
        """
        if self.M_function is None:
            element = ufl.VectorElement("DG", self.simulation_model.domain.ufl_cell(), 0, dim=3)
            self.M_function = fem.functionspace(self.simulation_model.domain, element)


    def export_displacement_rotation(self, FE_result, case:int=0):
        """
        Export results on dispacement and rotation

        Parameters:
        --------------
        FE_result: Finite element solution
        case: load case
        """
        disp = FE_result.sub(0).collapse()
        if case != 0:
            disp.name = "Displacement_" + str(case)
        else:
            disp.name = "Displacement"
        self.result_to_export.append(disp)
        theta = FE_result.sub(1).collapse()
        if case != 0:
            theta.name = "Rotation_" + str(case)
        else:
            theta.name = "Rotation"
        self.result_to_export.append(theta)


    def export_moments(self, FE_result, case:int=0):
        """
        Export results on moments

        Parameters:
        --------------
        u: Finite element solution
        case: load case
        """
        # Calculate moments in simulation model
        self.simulation_model.calculate_moments(FE_result)

        self.define_function_space()
        Moment = fem.Function(self.M_function)
        Moment_data = fem.Expression(self.simulation_model.moments,
                                     self.M_function.element.interpolation_points())
        Moment.interpolate(Moment_data)
        if case != 0:
            Moment.name = "Moment_" + str(case)
        else:
            Moment.name = "Moment"
        self.result_to_export.append(Moment)


    def export_macro_strain(self, case:int):
        """
        Export Macro Strain with defined case

        Parameters:
        -----------
        case : integer
            Strain loading case
        """
        self.define_function_space()
        Macro = fem.Function(self.M_function)
        Moment_data = fem.Expression(self.simulation_model.findImposedStrain(case), self.M_function.element.interpolation_points())
        Macro.interpolate(Moment_data)
        Macro.name = "Macro"
        self.result_to_export.append(Macro)

    def export_local_coordinates_system(self):
        """
        Export local coordinates axis system for all element in the domain
        """
        self.define_function_space()

        a1_ = fem.Function(self.M_function)
        a1__data = fem.Expression(self.simulation_model.BeamModel.a1, self.M_function.element.interpolation_points())
        a1_.interpolate(a1__data)
        a1_.name = "a1"
        self.result_to_export.append(a1_)

        a2_ = fem.Function(self.M_function)
        a2__data = fem.Expression(self.simulation_model.BeamModel.a2, self.M_function.element.interpolation_points())
        a2_.interpolate(a2__data)
        a2_.name = "a2"
        self.result_to_export.append(a2_)

        t_ = fem.Function(self.M_function)
        t__data = fem.Expression(self.simulation_model.BeamModel.t, self.M_function.element.interpolation_points())
        t_.interpolate(t__data)
        t_.name = "t"
        self.result_to_export.append(t_)

    def export_finalize(self, time:float=0.0):
        """
        Finalize the export process and close file
        """
        self._vtk.write_function(self.result_to_export, t=time)
        self._vtk.close()

    def write_function(self, time:float=0.0):
        """
        Write Function
        """
        print("Saving simulation results to file : " + self.name_file)
        self._vtk.write_function(self.result_to_export, t=time)

    def close_file(self):
        """
        Close the VTK file
        """
        self._vtk.close()


    def full_export(self, case:int):
        """
        Export all data for case i

        Parameters:
        -----------
        case : integer
            Strain loading case
        """
        self.export_displacement_rotation(self.simulation_model.u)
        self.export_moments(self.simulation_model.u)
        self.export_macro_strain(case)
        self.export_local_coordinates_system()
        self.export_finalize()

    def export_data_homogenization(self):
        """
        Export all data from 6 loading cases in homogenization
        """
        for case, simu_result in enumerate(self.simulation_model.saveDataToExport):
            self.export_displacement_rotation(simu_result)
            self.export_moments(simu_result)
            self.export_internal_force(simu_result)
            self.export_local_coordinates_system()
            self._vtk.write_function(self.result_to_export, t=case)
            self.result_to_export = []
        self._vtk.close()

    def export_internal_force(self, u, case:int=0):
        """
        Export results on Force

        Parameters:
        --------------
        u: Finite element solution
        case: load case
        """
        # Calculate moments in simulation model
        self.simulation_model.calculate_forces(u)

        self.define_function_space()
        Force = fem.Function(self.M_function)
        Force_data = fem.Expression(self.simulation_model.forces,
                                    self.M_function.element.interpolation_points())
        Force.interpolate(Force_data)
        if case != 0:
            Force.name = "Forces_" + str(case)
        else:
            Force.name = "Forces"
        self.result_to_export.append(Force)

    def export_vizualisation_3D(self, save_directory:str, number_point_ext: int = 8, mesh_radius_int: int = 3):
        """
        Export beam data on cylinder 3D
        """
        #TODO : a voir quoi faire ici
        def find_vector_director(dofmap, geometry):
            """Find the vector director of the beam element"""
            return [geometry.x[dofmap[1]][i] - geometry.x[dofmap[0]][i] for i in range(3)]

        # Clear old files
        clear_directory(save_directory)

        # Create mesh of beam element
        nameMesh = self.generate_beam_element(number_point_ext, mesh_radius_int)

        radius = self.simulation_model.BeamModel.radius.x.array
        nodePos = self.simulation_model.domain.geometry.x
        for i, dofmap in enumerate(self.simulation_model.domain.geometry.dofmap):
            if radius[i]>0.0:
                # Open mesh
                domainElement, markers, facets = io.gmshio.read_from_msh(nameMesh, self.simulation_model.domain.comm)
                # Apply modifications
                self.change_radius(domainElement, radius[i])
                self.resize_element(domainElement, dofmap, nodePos)
                self.rotate_element(domainElement, find_vector_director(dofmap, self.simulation_model.domain.geometry))
                self.move_element(domainElement, nodePos[dofmap[0]])
                # Interpolate data on mesh
                # To Do
                # Save mesh and data
                self.save_VTK_data(save_directory + str(i), domainElement)

        # Convert to open in Paraview in one step
        self.create_PVD_file(save_directory)

    def change_radius(self, domainElement, radius:float):
        """
        Adapt radius of the element

        Parameters:
        ------------
        domainElement: mesh
            Mesh of the domain to change
        radius: float
            Radius wanted for the element
        """
        domainElement.geometry.x[:, 0] *= radius
        domainElement.geometry.x[:, 1] *= radius

    def move_element(self, domainElement, nodePosition):
        """
       Move element

        Parameters:
        ------------
        domainElement: mesh
            Mesh of the domain to change
        nodePosition: array of dim 3
            Initial position wanted for the element
        """
        domainElement.geometry.x[:] += nodePosition

    def resize_element(self, domainElement, dofmap, nodePos):
        """
        Resize lenght of the element

        Parameters:
        ------------
        domainElement: mesh
            Mesh of the domain to change
        dofmap: Degree of freedom map of the element
        nodePos: domain.geometry.x
        """
        elementLenght = np.sqrt((nodePos[dofmap[1]][0]-nodePos[dofmap[0]][0])**2+(nodePos[dofmap[1]][1]-nodePos[dofmap[0]][1])**2+(nodePos[dofmap[1]][2]-nodePos[dofmap[0]][2])**2)
        domainElement.geometry.x[:, 2] *=  elementLenght

    def rotate_element(self, domainElement, newDirection):
        """
        Rotate element

        Parameters:
        ------------
        domainElement: mesh
            Mesh of the domain to change
        newDirection: array of dim 3
            Vector director wanted for the element
        """
        currentDirection = np.array([0, 0, 1])
        newDirection = np.array(newDirection)
        newDirection = newDirection / np.linalg.norm(newDirection)  # Normalize

        # Axis of rotation (cross product of current and new direction)
        axis = np.cross(currentDirection, newDirection)
        axis_length = np.linalg.norm(axis)
        if axis_length == 0:
            return  # The directions are already aligned

        axis = axis / axis_length  # Normalize the axis
        angle = np.arccos(np.dot(currentDirection, newDirection))  # Angle to rotate

        # Rodrigues' rotation formula components
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        I = np.eye(3)
        R = I + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

        # Apply rotation to all points
        domainElement.geometry.x[:, :3] = np.dot(domainElement.geometry.x[:, :3],R.T)

    def save_VTK_data(self, saveDirectory:str, domainElement):
        vtk = io.VTKFile(domainElement.comm, saveDirectory, "w")
        element = ufl.VectorElement("DG", domainElement.ufl_cell(), 0, dim=3)
        function = fem.functionspace(domainElement, element)
        A = fem.Function(function)
        Adata = fem.Constant(domainElement, (0.5, 0.5, 0.5))
        Force_data = fem.Expression(Adata,function.element.interpolation_points())
        A.interpolate(Force_data)
        A.name = "A"
        vtk.write_function(A)
        vtk.close()

    def find_beam_element_mesh_name(self, number_point_ext: int, mesh_radius_int: int):
        """
        Determine normalized mesh name and if it already exist

        Returns:
        ---------
        Path to mesh: String
        Mesh exist: boolean
        """
        meshName = f"BeamElement_{number_point_ext}{mesh_radius_int}.msh"
        script_dir = os.path.dirname(os.path.abspath(__file__))
        directory = os.path.join(script_dir, "Mesh", "Beam_Element")

        # Check if directory exists, if not, create it
        if not os.path.exists(directory):
            os.makedirs(directory)

        filePath = os.path.join(directory, meshName)
        return filePath, os.path.exists(filePath)

    def generate_beam_element(self, number_point_ext:int, mesh_radius_int:int):
        """
        Generate mesh of a beam element

        Parameters:
        ------------
        numberPointExt: integer
            Number of point on the mesh on the exterior of the beam
            Need to be modulo 4
        meshRadiusInt: integer
            Number of point on the mesh on the radius of the beam
        """
        radius = 1

        def modify_number_point_ext_beam(numberPointExt:int):
            """
            If necessary modify the number of point on the mesh of the exterior of the beam
            to work with the methodology
            """
            while numberPointExt % 4 != 0:
                numberPointExt += 1
            return numberPointExt

        nameMesh, meshExist = self.find_beam_element_mesh_name(number_point_ext, mesh_radius_int)
        if meshExist == False:
            # Initialize Gmsh and add a model
            gmsh.initialize()
            gmsh.model.add("BeamMesh")

            # Parameters
            number_point_ext = modify_number_point_ext_beam(number_point_ext)

            squareLength = radius/3
            point1 = gmsh.model.occ.addPoint(-squareLength,squareLength,0)
            point2 = gmsh.model.occ.addPoint(squareLength,squareLength,0)
            point3 = gmsh.model.occ.addPoint(squareLength,-squareLength,0)
            point4 = gmsh.model.occ.addPoint(-squareLength,-squareLength,0)

            line1 = gmsh.model.occ.addLine(point1,point2)
            line2 = gmsh.model.occ.addLine(point2,point3)
            line3 = gmsh.model.occ.addLine(point3,point4)
            line4 = gmsh.model.occ.addLine(point4,point1)

            pointCenter = gmsh.model.occ.addPoint(0,0,0)
            angleDegrees = math.radians(45)
            point1C = gmsh.model.occ.addPoint(-radius * math.cos(angleDegrees), radius * math.sin(angleDegrees), 0)
            point2C = gmsh.model.occ.addPoint(radius * math.cos(angleDegrees), radius * math.sin(angleDegrees), 0)
            point3C = gmsh.model.occ.addPoint(radius * math.cos(-angleDegrees), radius * math.sin(-angleDegrees), 0)
            point4C = gmsh.model.occ.addPoint(-radius * math.cos(-angleDegrees), radius * math.sin(-angleDegrees), 0)

            arc_circle1 = gmsh.model.occ.addCircleArc(point1C, pointCenter, point2C)
            arc_circle2 = gmsh.model.occ.addCircleArc(point2C, pointCenter, point3C)
            arc_circle3 = gmsh.model.occ.addCircleArc(point3C, pointCenter, point4C)
            arc_circle4 = gmsh.model.occ.addCircleArc(point4C, pointCenter, point1C)

            lineInt1 = gmsh.model.occ.addLine(point1C, point1)
            lineInt2 = gmsh.model.occ.addLine(point2C, point2)
            lineInt3 = gmsh.model.occ.addLine(point3C, point3)
            lineInt4 = gmsh.model.occ.addLine(point4C, point4)


            gmsh.model.occ.synchronize()
            meshCircleExt = int(number_point_ext / 4 + 1)
            gmsh.model.mesh.setTransfiniteCurve(arc_circle1, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(arc_circle2, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(arc_circle3, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(arc_circle4, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(lineInt1, mesh_radius_int)
            gmsh.model.mesh.setTransfiniteCurve(lineInt2, mesh_radius_int)
            gmsh.model.mesh.setTransfiniteCurve(lineInt3, mesh_radius_int)
            gmsh.model.mesh.setTransfiniteCurve(lineInt4, mesh_radius_int)

            gmsh.model.mesh.setTransfiniteCurve(line1, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(line2, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(line3, meshCircleExt)
            gmsh.model.mesh.setTransfiniteCurve(line4, meshCircleExt)

            loop = gmsh.model.occ.addCurveLoop([line1,line2,line3,line4])
            plane_surface = gmsh.model.occ.addPlaneSurface([loop])

            loopCircle1 = gmsh.model.occ.addCurveLoop([lineInt1,arc_circle1,lineInt2,line1])
            loopCircle2 = gmsh.model.occ.addCurveLoop([lineInt2,arc_circle2,lineInt3,line2])
            loopCircle3 = gmsh.model.occ.addCurveLoop([lineInt3,arc_circle3,lineInt4,line3])
            loopCircle4 = gmsh.model.occ.addCurveLoop([lineInt4,arc_circle4,lineInt1,line4])
            plane_surface_arc1 = gmsh.model.occ.addPlaneSurface([loopCircle1])
            plane_surface_arc2 = gmsh.model.occ.addPlaneSurface([loopCircle2])
            plane_surface_arc3 = gmsh.model.occ.addPlaneSurface([loopCircle3])
            plane_surface_arc4 = gmsh.model.occ.addPlaneSurface([loopCircle4])

            gmsh.model.occ.synchronize()
            gmsh.model.mesh.setTransfiniteSurface(plane_surface)
            gmsh.model.mesh.setTransfiniteSurface(plane_surface_arc1)
            gmsh.model.mesh.setTransfiniteSurface(plane_surface_arc2)
            gmsh.model.mesh.setTransfiniteSurface(plane_surface_arc3)
            gmsh.model.mesh.setTransfiniteSurface(plane_surface_arc4)
            dim = 2
            gmsh.model.mesh.setRecombine(dim, plane_surface)
            gmsh.model.mesh.setRecombine(dim, plane_surface_arc1)
            gmsh.model.mesh.setRecombine(dim, plane_surface_arc2)
            gmsh.model.mesh.setRecombine(dim, plane_surface_arc3)
            gmsh.model.mesh.setRecombine(dim, plane_surface_arc4)

            extrudeSurface = []
            for i, surface in enumerate([plane_surface,plane_surface_arc1,plane_surface_arc2,plane_surface_arc3,plane_surface_arc4]):
                extrudeSurface.append((dim,surface))

            extrudeVolume = gmsh.model.occ.extrude(extrudeSurface, 0, 0, 1, recombine=True, numElements=[1])

            gmsh.model.occ.synchronize()

            dim = 3
            volumeTags = [vol[1] for vol in gmsh.model.occ.getEntities(dim)]

            # Créer des groupes physiques pour une meilleure gestion dans les simulations
            gmsh.model.addPhysicalGroup(dim, volumeTags, 1)
            gmsh.model.setPhysicalName(dim, 1, "VolumeTag")

            # Générer le maillage
            gmsh.model.mesh.generate(dim)

            # Save and finalize
            gmsh.write(nameMesh)
            gmsh.finalize()
        return nameMesh

    def convert_mesh_for_paraview(self, nameMesh:str):
        domain, a, b = io.gmshio.read_from_msh(nameMesh, self.simulation_model.domain.comm)
        vtk = io.VTKFile(self.simulation_model.domain.comm, "Result/Try", "w")
        element = ufl.VectorElement("DG", domain.ufl_cell(), 0, dim=3)
        function = fem.functionspace(domain, element)
        A = fem.Function(function)
        Adata = fem.Constant(domain, (0.5, 0.5, 0.5))
        Force_data = fem.Expression(Adata,function.element.interpolation_points())
        A.interpolate(Force_data)
        A.name = "A"
        vtk.write_function(A)
        vtk.close()

    def create_PVD_file(self, vtk_directory, output_filename="#0_AllElements.pvd"):
        """
        Generate a PVD file to open all beam results in one time

        Parameters:
        ------------
        vtkDirectory: String
            Path to files
        outputFilename: String
            Name of output file
        """
        if not output_filename.endswith('.pvd'):
            output_filename += '.pvd'

        pvdContent = '<?xml version="1.0"?>\n'
        pvdContent += '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
        pvdContent += '    <Collection>\n'

        print(vtk_directory)
        vtkFiles = [f for f in os.listdir(vtk_directory) if f.endswith('.vtu')]
        vtkFiles.sort()

        for vtkFile in vtkFiles:
            pvdContent += f'        <DataSet timestep="0" part="0" file="{vtkFile}"/>\n'

        pvdContent += '    </Collection>\n'
        pvdContent += '</VTKFile>\n'

        # Écriture du fichier PVD
        with open(os.path.join(vtk_directory, output_filename), 'w') as pvdFile:
            pvdFile.write(pvdContent)
        print("Fichier à ouvrir sur Paraview : ")
        print(f'Fichier PVD généré : {output_filename}')


    def export_reaction_force(self, lattice_data):
        """
        Export reaction force
        """
        alreadyDone = []
        self.define_function_space()
        functionReaction = fem.FunctionSpace(self.simulation_model.domain, ("CG", 1))
        ReactionForce = fem.Function(functionReaction)
        dof_coordinates = functionReaction.tabulate_dof_coordinates()
        nodePosition = np.round(dof_coordinates, 3)
        triplet_tuples = [tuple(row) for row in nodePosition]
        dictNode = {triplet: idx for idx, triplet in enumerate(triplet_tuples)}
        print(dictNode)
        for cell in lattice_data.cells:
            for beam in cell.beams:
                for node in [beam.point1, beam.point2]:
                    if node.indexBoundary is not None and node.indexBoundary not in alreadyDone:
                        nodeIndices = dictNode.get(tuple(np.array([node.x,node.y,node.z])), None)
                        print(nodeIndices)
                        print(ReactionForce.vector[nodeIndices], node.reactionForceValue[:3])
                        ReactionForce.vector[nodeIndices] = node.reactionForceValue[:3]

        ReactionForce.interpolate(ReactionForce)
        ReactionForce.name = "ReactionForce"
        self.result_to_export.append(ReactionForce)