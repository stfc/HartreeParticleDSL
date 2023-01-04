Pydoc documentation
==========================

This is the main page for the pydoc generated documentation. Most functions
should be documented here or in subpages with their usage and descriptions.


Subpackages
-----------

.. toctree::
   :maxdepth: 0

   HartreeParticleDSL.IO_modules
   HartreeParticleDSL.backends
   HartreeParticleDSL.kernel_types

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

HartreeParticleDSL.HartreeParticleDSLExceptions module
------------------------------------------------------

This module contains various Exceptions that may be thrown in
HartreeParticleDSL programs.

.. automodule:: HartreeParticleDSL.HartreeParticleDSLExceptions
   :members:
   :undoc-members:
   :show-inheritance:

HartreeParticleDSL.c\_types module
----------------------------------

Deprecated.

.. automodule:: HartreeParticleDSL.c_types
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: HartreeParticleDSL
   :members:
   :undoc-members:
   :show-inheritance:
