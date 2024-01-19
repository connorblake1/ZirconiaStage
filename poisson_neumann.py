"""
Run with the default settings (no refinement)::

  sfepy-run sfepy/examples/diffusion/poisson_neumann.py

Refine the mesh twice::

  sfepy-run sfepy/examples/diffusion/poisson_neumann.py -O "'refinement_level' : 2"
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import output, Struct
from sfepy import data_dir

def post_process(out, pb, state, extend=False):
    """
    Calculate :math:`\nabla t` and compute boundary fluxes.
    """
    dv = pb.evaluate('ev_diffusion_velocity.i.Omega(m.K, t)', mode='el_avg',
                     verbose=False)
    out['dv'] = Struct(name='output_data', mode='cell',
                       data=dv, dofs=None)

    totals = nm.zeros(1)
    for gamma in ['Gamma_L', 'Gamma_R']:

        flux = pb.evaluate('ev_surface_flux.i.%s(m.K, t)' % gamma,
                           verbose=False)
        flux_data = (gamma, flux)
        totals += flux_data[1:]

        output('%8s flux: % 8.3f'
               % flux_data)

    # totals[2] = totals[0] / totals[1]
    # output('   total flux: % 8.3f length: % 8.3f flux/length: % 8.3f'
    #        % tuple(totals))

    return out

# filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'
filename_mesh = "stage.vtk"

materials = {
    'm' : ({'K' : 3.0 * nm.eye(3)},),
}

regions = {
    'Omega' : 'all',
    'Gamma_L' : ('vertices in (z > 5.799)', 'facet'),
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
        'i_max' : 3,
        'eps_a' : 1e-12,
    }),
}

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    # 'refinement_level' : 2,
    'post_process_hook' : 'post_process',
}
