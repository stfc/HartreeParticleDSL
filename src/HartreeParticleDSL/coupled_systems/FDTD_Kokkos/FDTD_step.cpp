#include "FDTD_step.hpp"
#include "FDTD_boundaries.hpp"


//const double c = 299792458.00000000;
//const double epsilon0 = 8.854187817620389850536563031710750e-12;

void update_eb_fields_half_1D(FDTD_field &field, const int nx, const int ng, const double dt, const double dx, 
                              update_e_field_functor &update_e_field,
                              update_b_field_functor &update_b_field,
                              Kokkos::RangePolicy<> rangepolicy){
    auto host_field = Kokkos::create_mirror_view(field.field);
    host_field(0).hdt = 0.5 * dt;
    host_field(0).hdtx = host_field(0).hdt / dx;
    host_field(0).cnx = host_field(0).hdtx * c*c;
    host_field(0).fac = host_field(0).hdt / epsilon0;
    Kokkos::deep_copy(field.field, host_field);

    Kokkos::parallel_for("update_e_field", rangepolicy, update_e_field);
    Kokkos::fence();

    //TODO
    efield_bcs(field.ex, field.ey, field.ez, nx, ng);

    Kokkos::parallel_for("update_b_field", rangepolicy, update_b_field);
    Kokkos::fence();

    //TODO
    bfield_bcs(field.bx, field.by, field.bz, nx, ng, true);
}

void update_eb_fields_final_1D(FDTD_field &field, const int nx, const int ng, const double dt, const double dx,
                              update_e_field_functor &update_e_field,
                              update_b_field_functor &update_b_field,
                              Kokkos::RangePolicy<> rangepolicy){

    auto host_field = Kokkos::create_mirror_view(field.field);
    host_field(0).hdt = 0.5 * dt;
    host_field(0).hdtx = host_field(0).hdt / dx;
    host_field(0).cnx = host_field(0).hdtx * c*c;
    host_field(0).fac = host_field(0).hdt / epsilon0;
    Kokkos::deep_copy(field.field, host_field);

    Kokkos::parallel_for("update_b_field", rangepolicy, update_b_field);
    Kokkos::fence();

    //TODO
    bfield_final_bcs(field.bx, field.by, field.bz, nx, ng);

    Kokkos::parallel_for("update_e_field", rangepolicy, update_e_field);
    Kokkos::fence();

    //TODO
    efield_bcs(field.ex, field.ey, field.ez, nx, ng);
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


void current_start(FDTD_field &field, const int nx, const int jng){
        current_start_functor current_start(field.jx, field.jy, field.jz);
        Kokkos::parallel_for("Current start", nx + 2*jng,
                current_start);
        Kokkos::fence();
}
