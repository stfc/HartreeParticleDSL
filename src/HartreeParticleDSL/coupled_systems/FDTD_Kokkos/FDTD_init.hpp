#ifndef HPDSL_KOKKOS_FDTD_INIT_H
#define HPDSL_KOKKOS_FDTD_INIT_H

#include "part.h"
#include "FDTD_field.hpp"

void kokkos_fdtd_initialize_1D(FDTD_field &field, int nx, int ng);

void kokkos_fdtd_cleanup_1D(FDTD_field &field);

void kokkos_fdtd_initialize_example_1D(FDTD_field &field, int nx, int ng);
#endif
