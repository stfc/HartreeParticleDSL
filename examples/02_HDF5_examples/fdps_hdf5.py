import ast as ast
import inspect as inspect

import HartreeParticleDSL.kernel_types.kernels as kernels
import HartreeParticleDSL.IO_modules.random_IO.random_IO as io_modules
import HartreeParticleDSL.IO_modules.HDF5_IO.hdf5_IO as hdf5_IO
from HartreeParticleDSL.backends.FDPS_backend.FDPS import FDPS
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.c_types import c_int, c_double


@kernels.perpart_interaction
def move_part(part1, config):
        part1.core_part.position.x = part1.core_part.position.x + part1.core_part.velocity[0] * config.dt
        

# Test includes
config = HartreeParticleDSL.Config()
config.add_element( "part_mass", "double" )
config.add_element( "time", "double")
config.add_element( "dt", "double" )
config.add_element( "cutoff", "double")                                                                                                                                                                                                                                                                                   
part = HartreeParticleDSL.Particle()
part.add_element( "force", "double[3]")
part.add_element( "two_dim_array", "double[4][3]")

io_module = io_modules.Random_Particles()
out_module = hdf5_IO.HDF5_IO()
out_module.add_output("x_pos", "core_part.position.x")

fdps = FDPS()
HartreeParticleDSL.set_backend(fdps)
HartreeParticleDSL.set_particle_type(part)
HartreeParticleDSL.set_config_type(config)
HartreeParticleDSL.set_io_modules(io_module, out_module)
HartreeParticleDSL.gen_code()

@kernels.main_declaration
def main():
    initialise(1000, "abc.def")
    config.time = 0.0
    config.dt = 0.1
    config.cutoff = 0.5
    while config.time < 1.0:
        invoke(move_part)
        config.time = config.time + config.dt
        println("", config.time,"\" \"", config.dt)
    write_output("File_out.hdf5")
    cleanup()
