HartreeParticleDSL.IO\_modules.PHDF5\_IO package
================================================

The PHDF5_IO module is the only module that currently supports PIR
backends (such as the preferred Cabana_PIR_backend) and MPI communication.

This module uses Parallel HDF5 to perform file operations, and allows users
to specify inputs and outputs with `add_input` and `add_output` respectively.

The module doesn't yet allow users to read values into the config from file
attributes, if this is a useful feature please open an issue in the github.

Submodules
----------

HartreeParticleDSL.IO\_modules.PHDF5\_IO.PHDF5\_IO module
---------------------------------------------------------

.. automodule:: HartreeParticleDSL.IO_modules.PHDF5_IO.PHDF5_IO
    :members: PHDF5_IO, add_input, add_output
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: HartreeParticleDSL.IO_modules.PHDF5_IO
    :members:
    :undoc-members:
    :show-inheritance:
