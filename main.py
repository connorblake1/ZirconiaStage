import cadquery as cq
from cadquery import exporters
import numpy as np
import datetime
from tqdm import tqdm
import gmsh
import meshio
import time
import json
import re
import os
# import sfepy is done implicitly when it's called in poisson_neumann.py
import matplotlib.pyplot as plt
start_time = time.time()
doRemesh = True
doFEM = True
# Key Dimensions
N_super = 2  # N = 13 #layers # TODO
M_super = 20 # M = 23 #pillar/layer#
outer_layer_diam_super = 34.9  # mm
inner_layer_diam_super = 10  # mm
outer_cut_super = .9 * .5 * outer_layer_diam_super
inner_cut_super = .6 * inner_layer_diam_super
layer_thickness_super = .6  # mm            (used in a_r)
pillar_rad_super = .45  # mm                w_p
pillar_height_super = .7  # mm              d_p
hex_center2center_super = 1.3  # mm         d_r (rung width is roughly half this)
hex_apothem_super = .435
dT_super = 100  # k # TODO input
top_height_super = layer_thickness_super + (N_super - 1) * (layer_thickness_super + pillar_height_super)
A_super = np.pi * (outer_layer_diam_super / 2) ** 2  # m^2
l_super = top_height_super  # mm
Factor_super = l_super / dT_super / A_super
metadata = {
    "N_Layers":N_super,
    "M_PillPerLayer":M_super,
    "outer_layer_diam_mm": outer_layer_diam_super,
    "inner_layer_diam_mm": inner_layer_diam_super,
    "outer_pillar_cutoff_mm":outer_cut_super,
    "inner_pillar_cutoff_mm":inner_cut_super,
    "layer_thickness_mm":layer_thickness_super,
    "pillar_radius_mm":pillar_rad_super,
    "pillar_height_mm":pillar_height_super,
    "hex_center2center_mm":hex_center2center_super,
    "hex_apothem_mm":hex_apothem_super,
    "top_face_z":top_height_super,
    "delta_T":dT_super,
    "solid_area_mm2":A_super,
    "l_over_deltaT_A":Factor_super,
    "pillar_mode":None,
    "scale_factor":3
}
# FOR LOADING IN METADATA FROM FILES
# with open("dummy_path", "r") as file:
#     metadata = json.load(file)
#     print("METADATA LOADED")
#     print(metadata)

def generateSTL(input_array,dummyCylinder=False,outerDiam=-1,extrusion=-1,fname='cylinder',selectionType="random",metadata_in=metadata,foldername=None,intermediateName=None):
    (N_Layers, M_PillPerLayer, outer_layer_diam_mm, inner_layer_diam_mm, outer_pillar_cutoff_mm,
     inner_pillar_cutoff_mm, layer_thickness_mm, pillar_radius_mm, pillar_height_mm,
     hex_center2center_mm, hex_apothem_mm, top_face_z, delta_T, solid_area_mm2,
     l_over_deltaT_A, pillar_mode,scalef) = metadata_in.values()
    metadata["pillar_mode"] = selectionType
    l = top_face_z
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
        cq.exporters.export(cyl, fname+'.stl')
        cq.exporters.export(cyl, fname+'.step')
        return extrusion/delta_T/A2
    print("Critical Factor:",l_over_deltaT_A)

    metadata["scale_factor"] = scalef
    basis_vectors = np.array([[0, 1], [np.sqrt(3) / 2, 0.5]])
    hex_vectors = np.array([[1/2,np.sqrt(3)/2],[-1/2,np.sqrt(3)/2]])
    pscale = hex_center2center_mm*scalef
    cella = scalef*hex_center2center_mm / np.sqrt(3)
    a_basis = np.array([0,1/np.sqrt(3)])
    b_basis = np.array([0,1/np.sqrt(3)])
    nums = int(outer_layer_diam_mm/cella/1.3)

    if intermediateName is None:
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
        centers = hexagonal_lattice_points(nums,outer_pillar_cutoff_mm,inner_pillar_cutoff_mm)
        stage = cq.Workplane('XY').circle(outer_layer_diam_mm/2).extrude(layer_thickness_mm)
        stage = stage.cut(cq.Workplane('XY').circle(inner_layer_diam_mm/2).extrude(layer_thickness_mm))
        print("Patterning Hexagons...")
        for i in tqdm(range(len(centers))):
            center = centers[i]
            stage = stage.cut(hexagonalprism(center,scalef*hex_apothem_mm,extrusionheight=layer_thickness_mm))
        stage_metadata = {
            "outer_layer_diam_mm": outer_layer_diam_mm,
            "inner_layer_diam_mm": inner_layer_diam_mm,
            "outer_pillar_cutoff_mm": outer_pillar_cutoff_mm,
            "inner_pillar_cutoff_mm": inner_pillar_cutoff_mm,
            "layer_thickness_mm": layer_thickness_mm,
            "hex_center2center_mm":hex_center2center_mm,
            "hex_apothem_mm":hex_apothem_mm,
            "solid_area_mm2":solid_area_mm2,
            "scale_factor":scalef
        }
        intName = foldername+"intermediatestage_"+current_datetime.strftime("%Y_%m_%d_%H_%M")
        with open(intName+".json", "w") as file:
            json.dump(stage_metadata, file, indent=4)
        cq.exporters.export(stage, intName+".step")
    else:
        with open(intermediateName + '.json', 'r') as file:
            loadedStageDict = json.load(file)
            allgood = True
            for key in loadedStageDict:
                if loadedStageDict[key] != metadata_in[key]:
                    allgood = False
                    break
            if not allgood:
                print("INTERMEDIATE STAGE NOT COMPATIBLE. EXITING...")
                exit()
            else:
                stage = cq.importers.importStep(intermediateName+".step")


    # stack layers
    print("Stacking...")
    structure = stage.translate((0,0,0))
    for i in tqdm(range(N_Layers)):
        stage_up = stage.translate((0,0,i*(layer_thickness_mm+pillar_height_mm)))
        structure = structure.union(stage_up)
    # place pillars
    def get_xy(tup):
        n,m,p = tup
        return pscale*np.array(n*hex_vectors[0] + m*hex_vectors[1] + a_basis + p*(b_basis-a_basis))+scalef*(basis_vectors[0]+basis_vectors[1])
    def is_valid(x,R_i,R_o):
        return np.linalg.norm(x) <= R_o and np.linalg.norm(x) >= R_i
    ind_range = 10 # TODO
    pillar_list = []
    metadataout = metadata_in.copy()
    if selectionType == "innerouter":
        checks = []
        close = []
        for xi in range(-ind_range,ind_range+1):
            for yi in range(-ind_range,ind_range+1):
                for bi in range(2):
                    pillar = [xi,yi,bi]
                    x, y = get_xy(pillar)
                    d2 = x ** 2 + y ** 2
                    xii, yii, d2i = np.round(x, 2), np.round(y, 2), np.round(d2, 2)
                    pxy = (xii, yii, d2i)
                    if is_valid(np.array([x,y]),inner_pillar_cutoff_mm,outer_pillar_cutoff_mm) and pxy not in checks:
                        checks.append(pxy)
                        lpxy = list(pxy)
                        lpxy.append(tuple(pillar))
                        close.append(lpxy)
        close_list = sorted(close, key=lambda x: x[2])
        for layer in tqdm(range(N_Layers-1)):
            layer_list = []
            for ii in range(M_PillPerLayer):
                if layer % 2 == 0:
                    x,y = close_list[ii][0],close_list[ii][1]
                    layer_list.append(close_list[ii][3])
                else:
                    x,y = close_list[-1-ii][0], close_list[-1-ii][1]
                    layer_list.append(close_list[-1-ii][3])
                pill = cq.Workplane("XY").circle(pillar_radius_mm).extrude(pillar_height_mm).translate(
                    (x, y, layer_thickness_mm + layer * (pillar_height_mm + layer_thickness_mm)))
                structure = structure.union(pill)
            pillar_list.append(layer_list)
    elif selectionType == "random": # TODO duplicates
        for layer in tqdm(range(N_Layers - 1)):
            done = []
            layer_list = []
            for _ in range(M_PillPerLayer):
                while True:
                    pillar = (np.random.randint(-ind_range, ind_range), np.random.randint(-ind_range, ind_range),
                              int(.5 + np.random.rand()))
                    x, y = get_xy(pillar)
                    val = is_valid(np.array([x, y]), inner_pillar_cutoff_mm, outer_pillar_cutoff_mm)
                    xi, yi = np.round(x, 2), np.round(y, 2)
                    pxy = (xi, yi)
                    if pxy in done:
                        continue
                    if val:
                        done.append(pxy)
                        layer_list.append(pillar)
                        pill = cq.Workplane("XY").circle(pillar_radius_mm).extrude(pillar_height_mm).translate(
                            (x, y, layer_thickness_mm + layer * (pillar_height_mm + layer_thickness_mm)))
                        structure = structure.union(pill)
                        break
            pillar_list.append(layer_list)
    elif selectionType == "zigzag":
        l1 = []
        l1c = []
        l2 = []
        l2c = []
        for xi in range(-ind_range, ind_range + 1,2):
            for yi in range(-ind_range, ind_range + 1,2):
                for bi in range(2):
                    pillar = [xi, yi, bi]
                    x, y = get_xy(pillar)
                    d2 = x ** 2 + y ** 2
                    xii, yii, d2i = np.round(x, 2), np.round(y, 2), np.round(d2, 2)
                    pxy = (xii, yii, d2i)
                    if is_valid(np.array([x, y]), inner_pillar_cutoff_mm, outer_pillar_cutoff_mm*1.1) and pxy not in l1c:
                        l1c.append(pxy)
                        lpxy = list(pxy)
                        lpxy.append(tuple(pillar))
                        l1.append(lpxy)
        for xi in range(-ind_range-1, ind_range + 1, 2):
            for yi in range(-ind_range-1, ind_range + 1, 2):
                for bi in range(2):
                    pillar = [xi, yi, bi]
                    x, y = get_xy(pillar)
                    d2 = x ** 2 + y ** 2
                    xii, yii, d2i = np.round(x, 2), np.round(y, 2), np.round(d2, 2)
                    pxy = (xii, yii, d2i)
                    if is_valid(np.array([x, y]), inner_pillar_cutoff_mm, outer_pillar_cutoff_mm*1.1) and pxy not in l2c:
                        l2c.append(pxy)
                        lpxy = list(pxy)
                        lpxy.append(tuple(pillar))
                        l2.append(lpxy)
        for layer in tqdm(range(N_Layers - 1)):
            layer_list = []
            lengthOfLists = min([M_PillPerLayer,len(l1),len(l2)])
            metadataout["M_PillPerLayer"] = lengthOfLists
            print(lengthOfLists)
            for ii in range(lengthOfLists):
                if layer % 2 == 0:
                    x, y = l1[ii][0], l1[ii][1]
                    layer_list.append(l1[ii][3])
                else:
                    x, y = l2[-1 - ii][0], l2[-1 - ii][1]
                    layer_list.append(l2[ii][3])
                pill = cq.Workplane("XY").circle(pillar_radius_mm).extrude(pillar_height_mm).translate(
                    (x, y, layer_thickness_mm + layer * (pillar_height_mm + layer_thickness_mm)))
                structure = structure.union(pill)
            pillar_list.append(layer_list)
    elif selectionType == 'straight':
        done = []
        layer_list = []
        for _ in tqdm(range(M_PillPerLayer)):
            while True:
                pillar = (np.random.randint(-ind_range, ind_range), np.random.randint(-ind_range, ind_range),
                          int(.5 + np.random.rand()))
                x, y = get_xy(pillar)
                val = is_valid(np.array([x, y]), inner_pillar_cutoff_mm, outer_pillar_cutoff_mm)
                xi,yi = np.round(x,2),np.round(y,2)
                pxy = (xi,yi)
                if pxy in done:
                    continue
                if val:
                    done.append(pxy)
                    layer_list.append(pillar)
                    for layer in range(N_Layers - 1):
                        pill = cq.Workplane("XY").circle(pillar_radius_mm).extrude(pillar_height_mm).translate(
                            (x, y, layer_thickness_mm + layer * (pillar_height_mm + layer_thickness_mm)))
                        structure = structure.union(pill)
                    break
        for layer in tqdm(range(N_Layers - 1)):
            pillar_list.append(layer_list)
    elif selectionType == "random_nostack":
        exclude_list = []
        for layer in tqdm(range(N_Layers - 1)):
            done = []
            layer_list = []
            for _ in range(M_PillPerLayer):
                while True:
                    pillar = (np.random.randint(-ind_range, ind_range), np.random.randint(-ind_range, ind_range),
                              int(.5 + np.random.rand()))
                    x, y = get_xy(pillar)
                    val = is_valid(np.array([x, y]), inner_pillar_cutoff_mm, outer_pillar_cutoff_mm)
                    xi,yi = np.round(x,2),np.round(y,2)
                    pxy = (xi,yi)
                    if pxy in done or pxy in exclude_list:
                        continue
                    if val:
                        done.append(pxy)
                        layer_list.append(pillar)
                        pill = cq.Workplane("XY").circle(pillar_radius_mm).extrude(pillar_height_mm).translate(
                            (x, y, layer_thickness_mm + layer * (pillar_height_mm + layer_thickness_mm)))
                        structure = structure.union(pill)
                        break
            exclude_list = done.copy()
            pillar_list.append(layer_list)

    metadataout["pillar_list"] = pillar_list
    with open(fname+".json", "w") as file:
        json.dump(metadataout, file, indent=4)
    cq.exporters.export(structure,fname+'.stl')
    cq.exporters.export(structure,fname+'.step')
    print("Generation Done: ",str(np.round(time.time()-start_time,2)))
    return l_over_deltaT_A
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
def computeFEMFlux(vtk_file,foldernamein): # TODO maybe broken?
    import os
    logfile = "FEM_output"
    folder_path = fr"C:\Users\theco\PycharmProjects\ZirconiaStage\{foldernamein}"
    script_path = folder_path + r"\poisson_neumann.py"
    sfepy_path = r"C:\ProgramData\anaconda3\Scripts\sfepy-run.exe"  # Replace with the actual path
    # Replace filename
    file_path = 'poisson_neumann.py'
    with open(file_path, 'r') as file:
        content = file.read()
    newpath = foldernamein+"\\\\"+vtk_file
    content = re.sub(r'filename_mesh\s*=\s*".*\.vtk"', fr'filename_mesh = "{newpath}"', content)
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
    ODn = 5
    Ln = 5
    ODs = np.linspace(.2,10,ODn)
    Ls = np.linspace(.4,10,Ln)
    Mn = 3
    meshes = np.flip(np.linspace(.15,.5,Mn))
    kappas = np.zeros((Ln,ODn,Mn))
    for i,L in enumerate(Ls):
        for j,OD in enumerate(ODs):
            for k,meshSize in enumerate(meshes):
                prefix = 'cylinder'
                meshSize = min(OD,meshSize)
                prefix += "_OD"+str(np.round(OD,2))+"_L"+str(np.round(L,2))+"_M" + str(np.round(meshSize,2))
                loTA = generateSTL(None,True,outerDiam=OD,extrusion=L,fname=prefix)
                print("Computing " + prefix)
                if doRemesh:
                    sname = prefix+".step"
                    generateMesh(sname,meshSize=meshSize,meshMin=meshSize/2,meshMax=2*meshSize)
                if doRemesh and doFEM:
                    vtkname = prefix+".vtk"
                    Q = computeFEMFlux(vtkname,"_") # TODO foldername
                    kappa = Q*loTA
                    kappas[i,j,k] = kappa
                    print("Computed Kappa_eff_z: "+prefix+ "    " + str(np.round(kappa,2)))
    for i in range(len(kappas)):
        print(kappas[i])
    for i in range(Ln):
        for j in range(Mn):
            plt.plot(ODs,kappas[i, :, j], label="L (mm) = " + str(np.round(Ls[i],2)) + " Mesh (mm) = " + str(np.round(meshes[j],2)))
        plt.xlabel('Outer Diameter (mm)')
        plt.ylabel('Kappa_eff_z (W/mK)')
        plt.ylim([2,3])
        plt.legend()
        plt.savefig('L'+str(np.round(Ls[i],2))+'.png')
        plt.clf()
    exit()
    for j in range(ODn):

        plt.plot(Ls,kappas[:, j], label="OD (mm)= " + str(np.round(ODs[j],2)))
        plt.xlabel('Length (mm)')
        plt.ylabel('Kappa_eff_z (W/mK)')
        plt.ylim([2, 3])
        plt.legend()
        plt.savefig('OD'+str(np.round(ODs[j],2))+'.png')
        plt.clf()

doMeshConvergence = False
if doMeshConvergence:
    meshSplits = np.flip(np.linspace(.15,.5,5))
    for meshS in meshSplits:
        prefix = 'stage'+"_M" +str(np.round(meshS,3)) + "_io"
        loTA = generateSTL(None,fname=prefix,selectionType="innerouter")
        print("Computing " + prefix)
        if doRemesh:
            sname = prefix+".step"
            generateMesh(sname,meshSize=meshS,meshMin=meshS/2,meshMax=2*meshS)
        if doRemesh and doFEM:
            vtkname = prefix+".vtk"
            Q = computeFEMFlux(vtkname,"")#TODO fix foldername
            kappa = Q*loTA
            print("Mesh: " + str(np.round(meshS,3)))
            print("Computed Kappa_eff_z: "+str(np.round(kappa,4)))

doGenerateNostack = True
if doGenerateNostack:
    Nstages = 1
    current_datetime = datetime.datetime.now()
    fileoutname = "SampleStages" + current_datetime.strftime("%Y_%m_%d_%H_%M")
    if not os.path.exists(fileoutname):
        os.makedirs(fileoutname)

    offset = 0
    for i in range(Nstages):
        prefix = "ZigZag_" + str(i+offset)
        pathname = fileoutname + "\\\\" + prefix
        print("Computing " + pathname)
        loTA = generateSTL(None,fname=pathname,selectionType="zigzag",metadata_in=metadata,foldername="",intermediateName="intermediatestage_2024_02_19_10_23")
        # if doRemesh:
        #     sname = pathname + ".step"
        #     meshS = .15
        #     generateMesh(sname, meshSize=meshS, meshMin=meshS / 2, meshMax=2 * meshS)
        # if doRemesh and doFEM:
        #     vtkname = prefix + ".vtk"
        #     Q = computeFEMFlux(vtkname,fileoutname)
        #     kappa = Q * loTA
        #     print("Mesh: " + str(np.round(meshS, 3)))
        #     print("Computed Kappa_eff_z: " + str(np.round(kappa, 4)))

        # generate column labels
        # generate file
        # save file, save pkl with stagelable
        # do some dummy computations
