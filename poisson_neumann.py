"""
Run with the default settings (no refinement)::

  sfepy-run sfepy/examples/diffusion/poisson_neumann.py

Refine the mesh twice::

  sfepy-run sfepy/examples/diffusion/poisson_neumann.py -O "'refinement_level' : 2"
 sfepy-run --list=terms    prints out all possible flags

"""
from __future__ import absolute_import
import numpy as nm
from sfepy.base.base import output, Struct
from sfepy import data_dir
def post_process(out, pb, state, extend=False):
    """
    Calculate :math:`\nabla t` and compute boundary fluxes.
    """
    totals = nm.zeros(1)
    for gamma in ['Gamma_L', 'Gamma_R']:
        flux = pb.evaluate('ev_surface_flux.i.%s(m.K, t)' % gamma,verbose=False)
        output("Q="+str(nm.round(flux,3)))
    return out

# filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
filename_mesh = "SampleStages2024_02_14_08_12\Stage_1.vtk"
# filename_mesh = 'cylinder.vtk'
materials = {
    'm' : ({'K' : 3.0 * nm.eye(3)},),
}
regions = {
    'Omega' : 'all',
    'Gamma_L' : ('vertices in (z > 5.798999999999999)', 'facet'),
    'Gamma_R' : ('vertices in (z < 0.001)', 'facet'),
}
fields = {
    'temperature' : ('real', 1, 'Omega', 1),
}
variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}
ebcs = {
    't1' : ('Gamma_L', {'t.0' : 100}),
    't2' : ('Gamma_R',{'t.0' : 0})
}
integrals = {
    'i' : 2
}
equations = {
    'Temperature' : """
           dw_diffusion.i.Omega(m.K, s, t)
         = 0"""
   # """
}

solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-12,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    # 'refinement_level' : 2,
    'post_process_hook' : 'post_process',
}
