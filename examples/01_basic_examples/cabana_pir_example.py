import ast as ast
import inspect as inspect

import HartreeParticleDSL.kernel_types.kernels as kernels
import HartreeParticleDSL.IO_modules.random_IO.random_IO as io_modules
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import config, part
from HartreeParticleDSL.c_types import c_int, c_double


@kernels.perpart_interaction
def move_part(part1: part, config: config):
        part1.core_part.position[0] = part1.core_part.position[0] + part1.velocity[0] * config.dt
        part1.core_part.position[1] = get_pointer(config.part_mass)
        create_variable(c_double, delta_x, 0.1)
        create_variable(c_double, part_vy, 0.1)
        create_variable(c_double, part_vz, 0.1)
        
# Test includes
#config = HartreeParticleDSL.Config()
config.add_element( "part_mass", "double" )
config.add_element( "time", "double")
config.add_element( "dt", "double" )
config.add_element( "cutoff", "double")

#part = HartreeParticleDSL.Particle()
part.add_element( "force", "double[3]")
part.add_element( "two_dim_array", "double[4][3]")
part.add_element( "velocity", "double[3]")

io_module = io_modules.Random_Particles()

cpir = Cabana_PIR()
HartreeParticleDSL.set_backend(cpir)
HartreeParticleDSL.set_particle_type(part)
HartreeParticleDSL.set_config_type(config)
HartreeParticleDSL.set_io_modules(io_module, io_module)
HartreeParticleDSL.gen_code()

@kernels.main_declaration
def main():
    initialise(particle_count=1000, filename="abc.def")
    config.time = 0.0
    config.dt = 0.1
    config.cutoff = 0.5
    create_variable(c_int, z, 2)
    while config.time < 1.0:
        invoke(move_part)
        config.time = config.time + config.dt
        println("", config.time,"\" \"", config.dt)
    cleanup()
