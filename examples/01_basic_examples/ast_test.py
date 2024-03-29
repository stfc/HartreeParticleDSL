import ast as ast
import inspect as inspect

import HartreeParticleDSL.kernel_types.kernels as kernels
import HartreeParticleDSL.IO_modules.random_IO.random_IO as io_modules
import HartreeParticleDSL.backends.C_AOS.C_AOS as C_AOS
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.c_types import c_int, c_double


@kernels.pairwise_interaction
def foo(part1, part2, r2, config):
    create_variable(c_int, "a", 0)
    create_variable(c_double, "b", 2.63)
    a = a + 2
    for i in range(3):
        part1.force[i] = b
        part1.two_dim_array[0][i] = b
    

@kernels.pairwise_interaction
def random_velocity(part1, part2, r2, config):
    for i in range(3):
        create_variable(c_double, "r", random_double())
        r = random_double() + random_double()
        part1.core_part.velocity[i] = part1.core_part.velocity[i] + r
        part2.core_part.velocity[i] = part2.core_part.velocity[i] - r

@kernels.perpart_interaction
def move_part(part1, config):
    for i in range(3):
        part1.core_part.position[i] = part1.core_part.position[i] + part1.core_part.velocity[i] * config.dt


# Test includes
config = HartreeParticleDSL.Config()
config.add_element( "part_mass", "double" )
config.add_element( "time", "double")
config.add_element( "dt", "double" )

part = HartreeParticleDSL.Particle()
part.add_element( "force", "double[3]")
part.add_element( "two_dim_array", "double[4][3]")

io_module = io_modules.Random_Particles()

HartreeParticleDSL.set_backend(C_AOS.C_AOS())
HartreeParticleDSL.set_particle_type(part)
HartreeParticleDSL.set_config_type(config)
HartreeParticleDSL.set_io_modules(io_module, io_module)
print("hello", type(HartreeParticleDSL.get_backend()))
HartreeParticleDSL.gen_code()

@kernels.main_declaration
def main():
    initialise(particle_count=1000, filename="abc.def")
    config.time = 0.0
    config.dt = 0.1
    create_variable(c_int, "z", 2)
    set_cutoff(0.5, C_AOS.CONSTANT)
    c_function()
    while config.time < 1.0:
        invoke(random_velocity)
        invoke(move_part)
        config.time = config.time + config.dt
        println("%f %f", config.time, config.dt)
    cleanup()
