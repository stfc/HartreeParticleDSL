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


Defining Kernels
----------------

The kernel types are defined in the ``HartreeParticleDSL.kernel_types.kernels``
module. There are a variety of types supported, however not all backends will
support all types of kernel, so check backend documentation for more detail.

.. automodule:: HartreeParticleDSL.kernel_types.kernels
    :members: perpart_interaction, main_declaration, pairwise_interaction
    

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
The inbuilt kernels package contains a number of kernels available for
use in HartreeParticleDSL code. These can just be imported and invoked
without the need to specify any extra functionality. If they require the
particle or config to contain specific fields, those will be documented
in their documentation


Boundary Conditions
^^^^^^^^^^^^^^^^^^^

The boundaries module contains the inbuilt boundary conditions supported
for particles in HartreeParticleDSL code. These will be expanded in the future,
including optimisations for 1D and 2D variants (note that the 3D version will
still work correctly for 1D/2D variants).
 .. automodule:: HartreeParticleDSL.inbuilt_kernels.boundaries
    :members:
    :undoc-members:
    :private-members:
    :special-members:
    :show-inheritance:
    :no-value:

Subpackages
-----------

.. toctree::
   :maxdepth: 1

   HartreeParticleDSL.IO_modules
   HartreeParticleDSL.backends
   HartreeParticleDSL.kernel_types
