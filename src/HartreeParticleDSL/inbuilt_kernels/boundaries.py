import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.c_types import c_double
from HartreeParticleDSL.HartreeParticleDSL import part, config

@kernels.perpart_interaction
def periodic_boundaries(part1: part, config: config):
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

periodic_boundaries = periodic_boundaries
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


@kernels.sink_boundary
def x_max_sink_boundary(part1: part, config: config):
    if part1.core_part.position.x > config.space.box_dims.x_max:
        part1.neighbour_part_deletion_flag = 1
    else:
        part1.neighbour_part_deletion_flag = 0

x_max_sink_boundary = x_max_sink_boundary
'''
Simple sink boundary that flags all particles above the x_max value
for deletion.
'''

x_max_sink_boundary.__doc__ = '''
Simple sink boundary that flags all particles above the x_max value
for deletion.
'''

@kernels.source_boundary
def x_min_neutral_source_boundary(part1: part, config: config):
    create_variable(c_double, v)
    create_variable(c_double, r_v1)
    create_variable(c_double, r_v2)
    create_variable(c_double, r_v3)
    create_variable(c_double, r_x)
    v = 2 * config.vth / sqrt(pi())
    r_v1 = random_number()
    r_v2 = random_number()
    r_v3 = random_number()
    r_x = random_number()

    part1.core_part.velocity[0] = v * myErfInf2(sqrt(r_v1))
    part1.core_part.velocity[0] = v * myErfInf2(sqrt(r_v2))
    part1.core_part.velocity[0] = v * myErfInf2(sqrt(r_v3))
    part1.core_part.position.x = -config.dt * part1.core_part.velocity.x * r_x
    part1.weight = config.n0 * config.vth * config.dt / 2.0 / sqrt(pi()) / 1000000.0
    part1.charge = 0.0

x_min_neutral_source_boundary.set_source_count(1000000)

x_min_neutral_source_boundary = x_min_neutral_source_boundary
'''
    Neutral source boundary that creates a number of neutral particles.
    Requires support for random numbers, MPI and an implementation of
    inverse error function.
'''

x_min_neutral_source_boundary.__doc__ = '''
    Neutral source boundary that creates a number of neutral particles.
    Requires support for random numbers, MPI and an implementation of
    inverse error function.
'''
