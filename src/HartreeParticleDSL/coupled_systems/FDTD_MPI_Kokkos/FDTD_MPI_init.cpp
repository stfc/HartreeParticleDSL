#include "FDTD_MPI_init.hpp"
#include "part.h"

void kokkos_fdtd_initialize_1D(FDTD_field &field, int nx, int ng){
    // Assuming kokkos is initialized
    field.ex = field_type("ex", nx + 2*ng);
    field.ey = field_type("ey", nx + 2*ng);
    field.ez = field_type("ez", nx + 2*ng);
    field.bx = field_type("bx", nx + 2*ng);
    field.by = field_type("by", nx + 2*ng);
    field.bz = field_type("bz", nx + 2*ng);
    field.jx = field_type("jx", nx + 2*ng);
    field.jy = field_type("jy", nx + 2*ng);
    field.jz = field_type("jz", nx + 2*ng);

    field.scatter_jx = scatter_field_type(config.jx);
    field.scatter_jy = scatter_field_type(config.jy);
    field.scatter_jz = scatter_field_type(config.jz);
}

void kokkos_fdtd_cleanup_1D(FDTD_field &field){

}


using host_mirror_type = field_type::HostMirror;


struct init_eb_arrays_functor{
    host_mirror_type _ex, _ey, _ez;
    host_mirror_type _bx, _by, _bz;

    KOKKOS_INLINE_FUNCTION
    init_eb_arrays_functor(host_mirror_type ex, host_mirror_type ey, host_mirror_type ez,
                           host_mirror_type bx, host_mirror_type by, host_mirror_type bz) :
                           _ex(ex), _ey(ey), _ez(ez), _bx(bx),
                            _by(by), _bz(bz) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        _ex(ix) = 0.0;
        _ey(ix) = 0.0;
        _ez(ix) = 0.0;
        _bx(ix) = 0.0;
        _by(ix) = 0.0;
        _bz(ix) = 2.1;
    }
};

struct init_j_arrays_functor{                                                                                                                                                  host_mirror_type _jx, _jy, _jz;

    KOKKOS_INLINE_FUNCTION
    init_j_arrays_functor(host_mirror_type jx, host_mirror_type jy, host_mirror_type jz) :
                      _jx(jx), _jy(jy), _jz(jz) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        _jx(ix) = 0.0;
        _jy(ix) = 0.0;
        _jz(ix) = 0.0;
    }
};

void kokkos_fdtd_initialize_example_1D(FDTD_field &field, int nx, int ng){
    
        // Create host mirrors for initialisation
        auto ex = Kokkos::create_mirror_view(field.ex);
        auto ey = Kokkos::create_mirror_view(field.ey);
        auto ez = Kokkos::create_mirror_view(field.ez);
        auto bx = Kokkos::create_mirror_view(field.bx);
        auto by = Kokkos::create_mirror_view(field.by);
        auto bz = Kokkos::create_mirror_view(field.bz);
        auto jx = Kokkos::create_mirror_view(field.jx);
        auto jy = Kokkos::create_mirror_view(field.jy);
        auto jz = Kokkos::create_mirror_view(field.jz);

        init_eb_arrays_functor init_eb(ex, ey, ez, bx, by, bz);
        init_j_arrays_functor init_j(jx, jy, jz);

        Kokkos::parallel_for("Init e&b arrays", Kokkos::RangePolicy<Kokkos::OpenMP>(0,nx+2*ng),
                             init_eb);
        Kokkos::parallel_for("Init j arrays", Kokkos::RangePolicy<Kokkos::OpenMP>(0,nx + 2*ng),
                             init_j);
        field.field.field_order = 2;
        field.field.fng = (double)(field.field.field_order) / 2.0;
        if (field.field.field_order == 2){
            field.field.cfl = 1.0;
        }else if(field.field.field_order == 4){
            field.field.cfl = 6.0/7.0;
        }else{
            field.field.cfl = 120.0/149.0;
        }
        Kokkos::fence();
        Kokkos::deep_copy(field.ex, ex);
        Kokkos::deep_copy(field.ey, ey);
        Kokkos::deep_copy(field.ez, ez);
        Kokkos::deep_copy(field.bx, bx);
        Kokkos::deep_copy(field.by, by);
        Kokkos::deep_copy(field.bz, bz);
        Kokkos::deep_copy(field.jx, jx);
        Kokkos::deep_copy(field.jy, jy);
        Kokkos::deep_copy(field.jz, jz);

}
