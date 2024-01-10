import cadquery as cq
import numpy as np
import gmsh
import meshio
doRemesh = True
doFEM = True
# Key Dimensions
pillar_rad = .45# mm
outer_layer_diam = 34.9# mm
inner_layer_diam = 10# mm
layer_thickness = .6# mm
pillar_height = .7# mm
hex_center2center = 1.3# mm
hex_apothem = .435
outer_cut,inner_cut = .9*.5*outer_layer_diam,.55*inner_layer_diam
M = 10 #M = 23 #pillar/layer
N = 5 # N = 13 #layers

scalef =7
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
structure = stage.translate((0,0,0))
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

# export
cq.exporters.export(structure,'stage.stl')

if doRemesh:
    # MESH IT
    cq.exporters.export(structure, 'stage.step')
    gmsh.initialize()
    gmsh.open('stage.step')
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    # gmsh.fltk.run()
    gmsh.write('stage.msh')
    mesh = meshio.read('stage.msh')
    mesh2 = meshio.Mesh(mesh.points,mesh.cells)
    mesh2.write("stage.vtk") # not bdf, mesh, inp, node
if doRemesh and doFEM:
    # RUN FEM
    import os
    folder_path = r"C:\Users\theco\PycharmProjects\ZirconiaStage"
    script_path = folder_path + r"\poisson_neumann.py"
    sfepy_path = r"C:\ProgramData\anaconda3\Scripts\sfepy-run.exe"  # Replace with the actual path
    os.system(f'cd /D "{folder_path}" && "{sfepy_path}" "{script_path}"')


