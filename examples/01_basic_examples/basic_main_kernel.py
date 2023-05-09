import HartreeParticleDSL.kernel_types.kernels as kernels
import HartreeParticleDSL.IO_modules.random_IO.random_IO as io_modules
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import config, part

# Add some time properties to keep track of in our simulation.
config.add_element("time", "double")
config.add_element("dt", "double")
config.add_element("step", "int")

# Use a random IO system for now (the particles don't do anything in this
# simulation anyway).
io_module = io_modules.Random_Particles()

# Load a backend for codegen.
cpir = Cabana_PIR()
# Setup the DSL for codegen
HartreeParticleDSL.set_backend(cpir)
HartreeParticleDSL.set_particle_type(part)
HartreeParticleDSL.set_config_type(config)
HartreeParticleDSL.set_io_modules(io_module, io_module)
# Call gen_code before declaring the main kernel.
HartreeParticleDSL.gen_code()


# Now create a simple main kernel.
# We will setup the simulation, set a dt and end time in the config,
# run steps until we reach the end time, and then print out the number
# of steps we ran for.
@kernels.main_declaration
def main():
    initialise(particle_count=0, filename="")
    config.time = 0.0
    config.dt = 0.1
    config.step = 0
    create_variable(double, end_time, 13.0)
    while config.time < end_time
         config.step = config.step + 1
         config.time = config.time + config.dt
    println("Ran for ", config.step, " steps.")
    cleanup()
