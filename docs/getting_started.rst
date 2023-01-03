Getting started with HartreeParticleDSL
=======================================

Installing HartreeParticleDSL
-----------------------------
HartreeParticleDSL should be downloaded from `https://github.com/stfc/HartreeParticleDSL`.

HartreeParticleDSL can be installed using pip, simply use `pip install .`
from the root directory to install HartreeParticleDSL.

It is not currently available outside of github.


Writing a HartreeParticleDSL script
------------------------------------

HartreeParticleDSL scripts are basic python scripts that describe the set of
operations to be executed on particles. Each independent set of operations is
combined into a kernel, e.g.

    @HartreeParticleDSl.kernels.perpart_interaction
    def my_kernel(part1: part, config: config):
        part1.core_part.position[0] += part1.core_part.velocity[0] * config.dt

Kernels can be executed in the main function using an `invoke(my_kernel)` call.
This invoke call tells the DSL to execute the operations across all particles in
the system at that time.

By default, the particle object only contains a `core_part` substructure, which contains
3D position and 3D velocity members. This position member is used by the DSL to
perform certain optimisations and to guarantee correctness of generated code,
so should be used by users to store their particles positions.

To add more members to the particle object, the user can call the `add_element` call
from the outermost-level of their python script, e.g.

    from HartreeParticleDSL.HartreeParticleDSL import part
    part.add_element("mass", "double")

would add a double precision member that could be accessed as `part1.mass` from
inside a kernel.

Similarly, by default the config object contains only a `space` member, which
itself is a structure containing a `box_dims` structure, which contains the 3D
min/max coordinates as `x_min`, `x_max`, etc. for the global box dimensions.

Other elements can be added to the space object in a similar way to the particle:

    from HartreeParticleDSL.HartreeParticleDSL import config
    config.add_element("dt", "double")

This would add a double precision member called `dt` to the config object. Unlike
particle members, config members can also be accessed from inside the main
definition.


Selecting options and backends
------------------------------
HartreeParticleDSL is designed to support a wide variety of backends, though the
current choices are limited.

Backends should be created by users in their script, and then provided to the DSL,
for example:

    from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
    import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
    cabana = Cabana_PIR()
    HartreeParticleDSL.set_backend(cabana)

Other options, such as `set_mpi` or `set_cuda` are available as options to the DSL
as well.

Each backend may have its own subset of features supported, so check the
documentation for a specific backend to check if a feature is supported.

Writing a main function
-----------------------

Main functions are used by most backends to determine the code that will be run
in a simulation by the generated code. Some backends in the future may choose
not to support main functions, as their purpose will be to generate kernel
code to fit into existing larger codebases.

The main function declaration usually contains:
1. Calls to `initialise` and `cleanup`. These generate any backend-specific code
   required for a backend to startup, setup a testcase and close down any
   libraries used during its runtime. Note that the `initialise` call takes
   `particle_count` and `filename` keyword arguments, and depending on the IO
    method chosen, one (or both) may be unused during the program runtime.
2. A timestepping loop (or similar). This usually will call some kernel
   invocations, do IO, and/or make calls to coupled systems.

For example, you might do:

    @kernels.main_declaration
    def main():
        initialise(particle_count=1, filename="myfile.hdf5")
        # NB This example is using HDF5 to read the file, so the
        # particle_count value is ignored.
        config.dt = 0.1
        config.time = 0.0
        while config.time < 1.0:
            invoke(mykernel)
            config.time = config.time + config.dt
        write_output("end.hdf5")
        cleanup()

for a very basic main function.

Running your script and using the output
----------------------------------------

Running a HartreeParticleDSL script is similar to any other python script,
e.g. `python myscript.py`. However, different backends may output the resulting
files to different places.

Old backends (that don't make use of the Particle intermediate representation)
will create a `part.h` file in the chosen location, but will output the
main code to standard out.

Newer backends will instead create a `part.h` and `code.cpp`, plus copy any other
required files into the relevant directory.

Currently, no backends create a Makefile or CMakeLists.txt file to compile the
output, so manual compilation is required. However, the new-style backends
will be required to provide build files soon to minimise the challenge of
compilation. These will also be added into the relevant directory.

Users can also specify the directory for output to be placed from their script:

    HartreeParticleDSL.set_output_dir("./output_dir")

This call would mean all outputs would be placed in `./output_dir`, and that
directory would be created if possible if it doesn't exist.
