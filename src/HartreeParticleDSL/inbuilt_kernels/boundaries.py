import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.c_types import c_double
from HartreeParticleDSL.HartreeParticleDSL import part, config

@kernels.perpart_interaction
def _periodic_boundaries(part1: part, config: config):
    create_variable(c_double, x_size, config.space.box_dims.x_max - config.space.box_dims.x_min)
    create_variable(c_double, y_size, config.space.box_dims.y_max - config.space.box_dims.y_min)
    create_variable(c_double, z_size, config.space.box_dims.z_max - config.space.box_dims.z_min)
    if part1.core_part.position.x > config.space.box_dims.x_max:
        part1.core_part.position.x = part1.core_part.position.x - x_size
    if part1.core_part.position.y > config.space.box_dims.y_max:
        part1.core_part.position.y = part1.core_part.position.y - y_size
    if part1.core_part.position.z > config.space.box_dims.z_max:
        part1.core_part.position.z = part1.core_part.position.z - z_size

    if part1.core_part.position.x < config.space.box_dims.x_min:
        part1.core_part.position.x = part1.core_part.position.x + x_size
    if part1.core_part.position.y < config.space.box_dims.y_min:
        part1.core_part.position.y = part1.core_part.position.y + y_size
    if part1.core_part.position.z < config.space.box_dims.z_min:
        part1.core_part.position.z = part1.core_part.position.z + z_size

periodic_boundaries = _periodic_boundaries
'''
Simple 3 dimensional periodic boundary condition function. Compares particle
position with config.space.box_dims and ensures the range is
dim_min <= position < dim_max.
'''

periodic_boundaries.__doc__ = '''
Simple 3 dimensional periodic boundary condition function. Compares particle
position with config.space.box_dims and ensures the range is
dim_min <= position < dim_max.
'''


periodic_boundaries_2 = _periodic_boundaries
'''Test'''
