#ifndef HPDSL_KOKKOS_FDTD_MPI_BOUNDARIES_H
#define HPDSL_KOKKOS_FDTD_MPI_BOUNDARIES_H
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "FDTD_MPI_field.hpp"

void field_bc_mpi(field_type field, int nx, int ng);
void efield_bcs(field_type ex, field_type ey, field_type ez,
                int nx, int ng);
void bfield_bcs(field_type bx, field_type by, field_type bz,
                int nx, int ng, bool mpi_only);
void processor_summation_boundaries_mpi(field_type field, int nx, int ng);
void current_bcs(field_type jx, field_type jy, field_type jz,
                 int nx, int ng);
void current_finish(field_type jx, field_type jy, field_type jz,
                    int nx, int ng);
void bfield_final_bcs(field_type bx, field_type by, field_type bz, int nx, int ng);

#endif
