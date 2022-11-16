#include "part.h"
#include "FDTD_MPI_step.hpp"
#include "FDTD_MPI_boundaries.hpp"

void update_eb_fields_half_1D(field_type &ex, field_type &ey, field_type &ez, int nx,
                              field_type &jx, field_type &jy, field_type &jz,
                              field_type &bx, field_type &by, field_type &bz,
                              double dt, double dx, field_struct_host host_field,
                              field_struct_type field, int ng,
                              update_e_field_functor &update_e_field,
                              update_b_field_functor &update_b_field,
                              Kokkos::RangePolicy<> rangepolicy){

    host_field(0).hdt = 0.5 * dt;
    host_field(0).hdtx = host_field(0).hdt / dx;
    host_field(0).cnx = host_field(0).hdtx * c*c;
    host_field(0).fac = host_field(0).hdt / epsilon0;
    Kokkos::deep_copy(field, host_field);

    Kokkos::parallel_for("update_e_field", rangepolicy, update_e_field);
    Kokkos::fence();

    efield_bcs(ex, ey, ez, nx, ng);

    Kokkos::parallel_for("update_b_field", rangepolicy, update_b_field);
    Kokkos::fence();

    bfield_bcs(bx, by, bz, nx, ng, true);
}

void update_eb_fields_final_1D(field_type &ex, field_type &ey, field_type &ez, int nx,
                               field_type &jx, field_type &jy, field_type &jz,
                               field_type &bx, field_type &by, field_type &bz,
                               double dt, double dx, field_struct_host host_field, field_struct_type field, int ng,
                               update_e_field_functor &update_e_field,
                               update_b_field_functor &update_b_field,
                               Kokkos::RangePolicy<> rangepolicy){

    host_field(0).hdt = 0.5 * dt;
    host_field(0).hdtx = host_field(0).hdt / dx;
    host_field(0).cnx = host_field(0).hdtx * c*c;
    host_field(0).fac = host_field(0).hdt / epsilon0;

    Kokkos::deep_copy(field, host_field);
    Kokkos::parallel_for("update_b_field", rangepolicy, update_b_field);
    Kokkos::fence();

    bfield_final_bcs(bx, by, bz, nx, ng);

    Kokkos::parallel_for("update_e_field", rangepolicy, update_e_field);
    Kokkos::fence();

    efield_bcs(ex, ey, ez, nx, ng);
}

struct current_start_functor{

    field_type _jx;
    field_type _jy;
    field_type _jz;

    KOKKOS_INLINE_FUNCTION
    current_start_functor(field_type jx, field_type jy, field_type jz) :
        _jx(jx), _jy(jy), _jz(jz) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        _jx(ix) = 0.0;
        _jy(ix) = 0.0;
        _jz(ix) = 0.0;
    }
};

void current_start(field_type jx, field_type jy, field_type jz, int nxlocal, int jng){
    current_start_functor current_start(jx, jy, jz);
    Kokkos::parallel_for("Current start", nxlocal + 2*jng, current_start);
    Kokkos::fence();
}
