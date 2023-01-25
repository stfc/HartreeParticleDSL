HartreeParticleDSL.coupled\_systems.FDTD\_MPI\_Kokkos package
=============================================================

A Finite Difference Time Domain solver with (optional) MPI using Kokkos. Supports
the Cabana_PIR backend.

Currently supports 1D grids with tophat interpolation.

This system contains a variety of .cpp and .hpp files, which are copied into the
chosen directory as part of the script running.

Submodules
----------

HartreeParticleDSL.coupled\_systems.FDTD\_MPI\_Kokkos.FDTD\_MPI\_Kokkos module
------------------------------------------------------------------------------

.. automodule:: HartreeParticleDSL.coupled_systems.FDTD_MPI_Kokkos.FDTD_MPI_Kokkos
    :members: set_interpolator, call_init_grid, setup_testcase, call_cleanup_grid,
     call_finish_initialisation_grid, call_eb_fields_first_halfstep, 
     call_eb_fields_final_halfstep, call_reset_current, call_finish_current,
     output_grid, call_interpolate_to_particles, gather_forces_to_grid
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: HartreeParticleDSL.coupled_systems.FDTD_MPI_Kokkos
    :members:
    :undoc-members:
    :show-inheritance:
