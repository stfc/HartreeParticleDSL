#ifndef HPDSL_KOKKOS_FDTD_MPI_FIELD_H
#define HPDSL_KOKKOS_FDTD_MPI_FIELD_H
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_ScatterView.hpp>

using MemorySpace = Kokkos::HostSpace;
using field_type = Kokkos::View<double* , MemorySpace>;
using scatter_field_type = Kokkos::Experimental::ScatterView<double*>;

struct field_struc{
    double hdt;
    double hdtx;
    double cnx;
    double fac;
    int field_order;
    double fng;
    double cfl;
    double x_grid_min_local;
    double x_grid_max_local;
    double x_min_local;
    double x_max_local;
};

struct FDTD_field{
    struct field_struc field;
    field_type ex;
    field_type ey;
    field_type ez;
    field_type bx;
    field_type by;
    field_type bz;
    field_type jx;
    field_type jy;
    field_type jz;
    scatter_field_type scatter_jx;
    scatter_field_type scatter_jy;
    scatter_field_type scatter_jz;
    int nxglobal;
    int nx;
    int min_local_cell;
    int max_local_cell;
    int ng;
    int jng;
    double dx;
};
#endif
