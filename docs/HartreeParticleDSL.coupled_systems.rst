HartreeParticleDSL.coupled\_systems package
===========================================

Couple systems provide non-particle specific functionality to
HartreeParticleDSL programs.

Currently all implemented coupled systems are dedicated to solving
and interpolating to FDTD grid solvers. These grid solvers are currently
standalone implementations, however future implementations here do not
have to be and could interact with other existing codebases.

The FDTD_MPI_Kokkos backend is the most up-to-date, support Particle_IR
and the Cabana_PIR setups, while the other backends are legacy and support
older variants with less functionality.


Subpackages
-----------

.. toctree::

    HartreeParticleDSL.coupled_systems.FDTD_MPI_Kokkos
    HartreeParticleDSL.coupled_systems.FDTD
    HartreeParticleDSL.coupled_systems.FDTD_Kokkos

Module contents
---------------

.. automodule:: HartreeParticleDSL.coupled_systems
    :members:
    :undoc-members:
    :show-inheritance:
