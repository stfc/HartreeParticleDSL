#ifndef HPDSL_KOKKOS_FDTD_BOUNDARIES_H
#define HPDSL_KOKKOS_FDTD_BOUNDARIES_H
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "FDTD_field.hpp"

//Periodic only
struct field_bc{

    field_type _field;
    int _ng;
    int _nx;

    KOKKOS_INLINE_FUNCTION
    field_bc(field_type field, int nx, int ng): _field(field),
        _ng(ng), _nx(nx){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        _field(_ng+_nx+ix) = _field(_ng+ix);
        _field(ix) = _field(_nx+ix);
    }

};

//Periodic
struct processor_summation_boundaries{

    field_type _field;
    int _ng;
    int _nx;

    KOKKOS_INLINE_FUNCTION
    processor_summation_boundaries(field_type field, int nx, int ng): _field(field),
        _ng(ng), _nx(nx){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        _field(_ng+ix) += _field(_ng+_nx+ix);
        _field(_nx+_ng-ix) += _field(_ng-ix);
    }

};

void efield_bcs(field_type ex, field_type ey, field_type ez,
                int nx, int ng);

void bfield_bcs(field_type bx, field_type by, field_type bz,
                int nx, int ng, bool mpi_only);

void current_bcs(field_type jx, field_type jy, field_type jz,
                 int nx, int ng);

void current_finish(field_type jx, field_type jy, field_type jz,
                 int nx, int ng);

void bfield_final_bcs(field_type bx, field_type by, field_type bz, int nx, int ng);
#endif
