#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "FDTD_boundaries.hpp"

void efield_bcs(field_type ex, field_type ey, field_type ez,
                int nx, int ng){
    field_bc ex_fbc(ex, nx, ng);
    field_bc ey_fbc(ey, nx, ng);
    field_bc ez_fbc(ez, nx, ng);

    auto rp = Kokkos::RangePolicy<>(0, ng);
    Kokkos::parallel_for("ex_bcs", rp, ex_fbc);
    Kokkos::parallel_for("ey_bcs", rp, ey_fbc);
    Kokkos::parallel_for("ez_bcs", rp, ez_fbc);
    Kokkos::fence();
}

void bfield_bcs(field_type bx, field_type by, field_type bz,
                int nx, int ng, bool mpi_only){

    field_bc bx_fbc(bx, nx, ng);
    field_bc by_fbc(by, nx, ng);
    field_bc bz_fbc(bz, nx, ng);
    auto rp = Kokkos::RangePolicy<>(0, ng);
    Kokkos::parallel_for("bx_bcs", rp, bx_fbc);
    Kokkos::parallel_for("by_bcs", rp, by_fbc);
    Kokkos::parallel_for("bz_bcs", rp, bz_fbc);
    Kokkos::fence();
}

void current_bcs(field_type jx, field_type jy, field_type jz,
                 int nx, int ng){

    processor_summation_boundaries jx_sum(jx, nx, ng);
    processor_summation_boundaries jy_sum(jy, nx, ng);
    processor_summation_boundaries jz_sum(jz, nx, ng);
    auto rp = Kokkos::RangePolicy<>(0, ng);
    Kokkos::parallel_for("jx_summation_boundaries", rp, jx_sum);
    Kokkos::parallel_for("jy_summation_boundaries", rp, jy_sum);
    Kokkos::parallel_for("jz_summation_boundaries", rp, jz_sum);
    Kokkos::fence();
}

void current_finish(field_type jx, field_type jy, field_type jz,
                 int nx, int ng){

    current_bcs(jx, jy, jz, nx, ng);

    field_bc jx_fbc(jx, nx, ng);
    field_bc jy_fbc(jy, nx, ng);
    field_bc jz_fbc(jz, nx, ng);
    auto rp = Kokkos::RangePolicy<>(0, ng);
    Kokkos::parallel_for("jx_bcs", rp, jx_fbc);
    Kokkos::parallel_for("jy_bcs", rp, jy_fbc);
    Kokkos::parallel_for("jz_bcs", rp, jz_fbc);
    Kokkos::fence();
}

void bfield_final_bcs(field_type bx, field_type by, field_type bz, int nx, int ng){
    //Ignore update_laser_omegs
    bfield_bcs(bx, by, bz, nx, ng, true);
    //x_min_boundary ignore
    //x_max_boundary ignore                                                                                                                                                    //
    bfield_bcs(bx, by, bz, nx, ng, false);
}
