#ifndef HPDSL_KOKKOS_FDTD_FIELD_H
#define HPDSL_KOKKOS_FDTD_FIELD_H
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

using MemorySpace = Kokkos::HostSpace;
using field_type = Kokkos::View<double* , MemorySpace>;

struct field{
    double hdt;
    double hdtx;
    double cnx;
    double fac;
    int field_order;
    double fng;
    double cfl;
    double x_grid_min_local;
    double x_grid_max_local;
};

struct FDTD_field{
    struct field field;
    field_type ex;
    field_type ey;
    field_type ez;
    field_type bx;
    field_type by;
    field_type bz;
    field_type jx;
    field_type jy;
    field_type jz;
    int nx;
    int ng;
    int jng;
};
#endif
