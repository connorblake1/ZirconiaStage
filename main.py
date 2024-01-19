import cadquery as cq
import numpy as np
import gmsh
import meshio
import time
# import sfepy is done implicitly when it's called in poisson_neumann.py
import matplotlib.pyplot as plt
start_time = time.time()
doRemesh = True
doFEM = True
# Key Dimensions
def generateSTL(input_array): # TODO make as actual inputs
    outer_layer_diam = 34.9# mm
    inner_layer_diam = 10# mm
    outer_cut= .9*.5*outer_layer_diam
    inner_cut = .55*inner_layer_diam
    layer_thickness = .6# mm
    pillar_rad = .45# mm
    pillar_height = .7# mm
    hex_center2center = 1.3# mm
    hex_apothem = .435
    M = 10 #M = 23 #pillar/layer
    N = 5 # N = 13 #layers

    top_height = layer_thickness+(N-1)*(layer_thickness+pillar_height)
    print("Top:",top_height)
    A = (.001*.001)*np.pi*(outer_layer_diam/2)**2 # m^2
    dT = 100 # k
    l = .001*top_height # m
    Factor = l/dT/A
    print("Critical Factor:",Factor)

    scalef = 7
    basis_vectors = np.array([[0, 1], [np.sqrt(3) / 2, 0.5]])
    hex_vectors = np.array([[1/2,np.sqrt(3)/2],[-1/2,np.sqrt(3)/2]])
    pscale = hex_center2center*scalef
    cella = scalef*hex_center2center / np.sqrt(3)
    a_basis = np.array([0,1/np.sqrt(3)])
    b_basis = np.array([0,1/np.sqrt(3)])
    nums = int(outer_layer_diam/cella/1.3)

    # punch hex
    def hexagonal_lattice_points(N, R_o,R_i):
        lattice_points = []
        for i in range(-N,N):
            for j in range(-N,N):
                    lattice_points.append(basis_vectors[0] * cella * i + basis_vectors[1]*cella*j)
        lattice_points = np.array(lattice_points)
        x = lattice_points[:,0]
        y = lattice_points[:,1]
        mask = (np.sqrt(x ** 2 + y ** 2) <= R_o) & (np.sqrt(x ** 2 + y ** 2) > R_i)
        coordinates = np.column_stack([x[mask], y[mask]])
        return coordinates
    def hexagonalprism(center,apothem,extrusionheight):
        hexagon = (
            cq.Workplane("XY")
                .center(center[0],center[1])
                .polygon(6,2/np.sqrt(3)*apothem) #TODO rotate
                .extrude(extrusionheight)
        )
        return hexagon
    centers = hexagonal_lattice_points(nums,outer_cut,inner_cut)
    stage = cq.Workplane('XY').circle(outer_layer_diam/2).extrude(layer_thickness)
    stage = stage.cut(cq.Workplane('XY').circle(inner_layer_diam/2).extrude(layer_thickness))
    for center in centers:
        stage = stage.cut(hexagonalprism(center,scalef*hex_apothem,extrusionheight=layer_thickness))
    # stack layers
    structure = stage.translate((0,0,0)) # TODO speed up by storing this and then adding and saving intenrally after hexes punched
    for i in range(N):
        stage_up = stage.translate((0,0,i*(layer_thickness+pillar_height)))
        structure = structure.union(stage_up)
    # place pillars
    def get_xy(tup):
        n,m,p = tup
        return pscale*np.array(n*hex_vectors[0] + m*hex_vectors[1] + a_basis + p*(b_basis-a_basis))+scalef*(basis_vectors[0]+basis_vectors[1])
    def is_valid(x,R_i,R_o):
        return np.linalg.norm(x) <= R_o and np.linalg.norm(x) >= R_i
    ind_range = 3
    for layer in range(N-1):
        done = []
        for _ in range(M):
            while True:
                pillar = [np.random.randint(-ind_range, ind_range), np.random.randint(-ind_range,ind_range),int(.5+np.random.rand())]
                if pillar in done:
                    continue
                x,y = get_xy(pillar)
                val = is_valid(np.array([x,y]),inner_cut,outer_cut)
                if val:
                    done.append(pillar)
                    pill = cq.Workplane("XY").circle(pillar_rad).extrude(pillar_height).translate((x,y,layer_thickness+layer*(pillar_height+layer_thickness)))
                    structure = structure.union(pill)
                    break
    cq.exporters.export(structure,'stage.stl') # TODO make return readable file handles
    cq.exporters.export(structure, 'stage.step')
    print("Generation Done: ",str(np.round(time.time()-start_time,2)))
def generateMesh(step_file):
    gmsh.initialize() # starts all gmsh code
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.open(step_file)
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0),.4) # get entities gets all the points which are what matter, size of mesh element
    gmsh.option.setNumber("Mesh.MeshSizeMin", .4)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 1)
    gmsh.model.mesh.generate(3)
    # gmsh.fltk.run()
    gmsh.write('stage.msh')
    # meshio is just for reconverting the file
    mesh = meshio.read('stage.msh')
    mesh2 = meshio.Mesh(mesh.points,mesh.cells)
    mesh2.write("stage.vtk") # not bdf, mesh, inp, node
    print("Meshing Done: ",str(np.round(time.time()-start_time,2)))
def computeFEM(vtk_file):
    import os
    folder_path = r"C:\Users\theco\PycharmProjects\ZirconiaStage"
    script_path = folder_path + r"\poisson_neumann.py"
    sfepy_path = r"C:\ProgramData\anaconda3\Scripts\sfepy-run.exe"  # Replace with the actual path
    os.system(f'cd /D "{folder_path}" && "{sfepy_path}" "{script_path}"') # TODO stage.vtk as input
    print("Total Elapsed:",str(np.round(time.time()-start_time,2)))
generateSTL(None)
if doRemesh:
    generateMesh("stage.step")
if doRemesh and doFEM:
    computeFEM('stage.vtk')


