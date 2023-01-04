User documentation for HartreeParticleDSL
=========================================

This is the main page for the pydoc generated documentation. Most functions
that users of HartreeParticleDSL should need are documented here.


HartreeParticleDSL.HartreeParticleDSL module
--------------------------------------------

This module contains the main DSL entry points and documentation. The 
``_HartreeParticleDSL`` class should not be used directly, and the visible
functions handle accessing that class.

.. automodule:: HartreeParticleDSL.HartreeParticleDSL
   :members: set_particle_type, set_config_type, set_backend, set_io_modules,
    set_output_dir, set_mpi, get_mpi, set_cuda, get_cuda, Particle,
    Config, part, config
   :undoc-members:
   :show-inheritance:

HartreeParticleDSLExceptions
------------------------------------------------------

This module contains various Exceptions that may be thrown in
HartreeParticleDSL programs.

.. automodule:: HartreeParticleDSL.HartreeParticleDSLExceptions
   :members:
   :undoc-members:
   :show-inheritance:


Inbuilt Kernels
---------------
The inbuild kernels package contains a number of kernels available for
use in HartreeParticleDSL code. These can just be imported and invoked
without the need to specify any extra functionality. If they require the
particle or config to contain specific fields, those will be documented
in their documentation


Boundary Conditions
^^^^^^^^^^^^^^^^^^^
.. automodule:: HartreeParticleDSL.inbuilt_kernels
    :members: boundaries
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::
   :maxdepth: 1

   HartreeParticleDSL.IO_modules
   HartreeParticleDSL.backends
   HartreeParticleDSL.kernel_types
