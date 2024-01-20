import cadquery as cq
import numpy as np
import gmsh
import meshio
import time
import re
# import sfepy is done implicitly when it's called in poisson_neumann.py
import matplotlib.pyplot as plt
start_time = time.time()
doRemesh = True
doFEM = True
# Key Dimensions
def generateSTL(input_array,dummyCylinder=False,outerDiam=-1,extrusion=-1,fname='cylinder'):
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
    A = np.pi*(outer_layer_diam/2)**2 # m^2
    dT = 100  # k # TODO input
    l = top_height  # m
    Factor = l/dT/A
    print("Critical Factor:",Factor)
    # Replace BCs
    if dummyCylinder:
        new_number = str(extrusion-.001)
    else:
        new_number = str(l-.001)
    with open('poisson_neumann.py', 'r') as file:
        lines = file.readlines()
    for i, line in enumerate(lines):
        if 'Gamma_L' in line and 'vertices in (z >' in line:
            match = re.search(r'z > (\d+(\.\d+)?)', line)
            if match:
                old_number = match.group(1)
                lines[i] = line.replace(old_number, str(new_number))
    with open('poisson_neumann.py', 'w') as file:
        file.writelines(lines)
    if dummyCylinder:
        cyl = cq.Workplane('XY').circle(outerDiam/2).extrude(extrusion)
        A2 = np.pi*(outerDiam/2)**2
        cq.exporters.export(cyl, fname+"_OD"+str(np.round(outerDiam,2))+"_L"+str(np.round(extrusion,2))+'.stl')
        cq.exporters.export(cyl, fname+"_OD"+str(np.round(outerDiam,2))+"_L"+str(np.round(extrusion,2))+'.step')
        return extrusion/dT/A2
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
                .polygon(6,2/np.sqrt(3)*apothem)
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
    cq.exporters.export(structure,'stage.stl')
    cq.exporters.export(structure,'stage.step')
    print("Generation Done: ",str(np.round(time.time()-start_time,2)))
    return Factor
def generateMesh(step_file,meshSize=.2,meshMin=.1,meshMax=.5):
    gmsh.initialize() # starts all gmsh code
    fname = step_file.rsplit('.', 1)[0]
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.open(step_file)
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0),meshSize) # get entities gets all the points which are what matter, size of mesh element
    gmsh.option.setNumber("Mesh.MeshSizeMin", meshMin)
    gmsh.option.setNumber("Mesh.MeshSizeMax", meshMax)
    gmsh.model.mesh.generate(3)
    # gmsh.fltk.run()
    gmsh.write(fname+'.msh')
    # meshio is just for reconverting the file
    mesh = meshio.read(fname+'.msh')
    mesh2 = meshio.Mesh(mesh.points,mesh.cells)
    mesh2.write(fname+".vtk") # not bdf, mesh, inp, node
    print("Meshing Done: ",str(np.round(time.time()-start_time,2)))
def computeFEMFlux(vtk_file):
    import os
    logfile = "FEM_output"
    folder_path = r"C:\Users\theco\PycharmProjects\ZirconiaStage"
    script_path = folder_path + r"\poisson_neumann.py"
    sfepy_path = r"C:\ProgramData\anaconda3\Scripts\sfepy-run.exe"  # Replace with the actual path
    # Replace filename
    file_path = 'poisson_neumann.py'
    with open(file_path, 'r') as file:
        content = file.read()
    content = re.sub(r'filename_mesh\s*=\s*".*\.vtk"', f'filename_mesh = "{vtk_file}"', content)
    with open(file_path, 'w') as file:
        file.write(content)
    # Execute
    os.system(f'cd /D "{folder_path}" && "{sfepy_path}" "{script_path}" --log "{logfile}"') # TODO stage.vtk as input
    with open(logfile, 'r') as file:
        lines_with_sfepy_q = [abs(float(line.strip()[9:])) for line in file if line.startswith("sfepy: Q=")]
    print(lines_with_sfepy_q)
    print("Total Elapsed:",str(np.round(time.time()-start_time,2)))
    return sum(lines_with_sfepy_q)/2

doBenchmarking = False
if doBenchmarking:
    ODn = 10
    Ln = 10
    ODs = np.linspace(.2,10,ODn)
    Ls = np.linspace(.4,10,Ln)
    kappas = np.zeros((Ln,ODn))
    for i,L in enumerate(Ls):
        for j,OD in enumerate(ODs):
            prefix = 'cylinder'
            loTA = generateSTL(None,True,outerDiam=OD,extrusion=L)
            prefix += "_OD"+str(np.round(OD,2))+"_L"+str(np.round(L,2))
            print("Computing " + prefix)
            if doRemesh:
                sname = prefix+".step"
                generateMesh(sname)
            if doRemesh and doFEM:
                vtkname = prefix+".vtk"
                Q = computeFEMFlux(vtkname)
                kappa = Q*loTA
                kappas[i,j] = kappa
                print("Computed Kappa_eff_z: "+str(np.round(kappa,2)))
    for i in range(len(kappas)):
        print(kappas[i])
    for i in range(Ln):
        plt.plot(ODs,kappas[i, :], label="L = " + str(np.round(Ls[i],2)) )
        plt.xlabel('Outer Diameters')
        plt.ylabel('Kappa_eff_z')
        plt.ylim([2,3])
        plt.legend()
        plt.savefig('L'+str(np.round(Ls[i],2))+'.png')
        plt.clf()

    for j in range(ODn):
        plt.plot(Ls,kappas[:, j], label="OD = " + str(np.round(ODs[j],2)))
        plt.xlabel('Lengths')
        plt.ylabel('Kappa_eff_z')
        plt.ylim([2, 3])
        plt.legend()
        plt.savefig('OD'+str(np.round(ODs[j],2))+'.png')
        plt.clf()

doMeshConvergence = True
if doMeshConvergence:
    meshSplits = np.flip(np.linspace(.05,.5,10))
    for meshS in meshSplits:
        prefix = 'stage'
        loTA = generateSTL(None)
        print("Computing " + prefix)
        if doRemesh:
            sname = prefix+".step"
            generateMesh(sname,meshSize=meshS,meshMin=meshS/2,meshMax=2*meshS)
        if doRemesh and doFEM:
            vtkname = prefix+".vtk"
            Q = computeFEMFlux(vtkname)
            kappa = Q*loTA
            print("Mesh: " + str(np.round(meshS,3)))
            print("Computed Kappa_eff_z: "+str(np.round(kappa,4)))
