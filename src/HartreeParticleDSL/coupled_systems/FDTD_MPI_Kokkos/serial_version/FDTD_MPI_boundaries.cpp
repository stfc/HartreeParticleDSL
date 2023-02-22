#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "FDTD_MPI_boundaries.hpp"

void field_bc_mpi(field_type field, int nx, int ng){
    //Copy ends of real domain into ghost cells.
    for(int i = 0; i < ng; i++){
        field(ng+nx+i) = field(ng+i);
        field(i) = field(nx+i);
    }
}

void efield_bcs(field_type ex, field_type ey, field_type ez,
                int nx, int ng){
 field_bc_mpi(ex, nx, ng);
 field_bc_mpi(ey, nx, ng);
 field_bc_mpi(ez, nx, ng);
}

void bfield_bcs(field_type bx, field_type by, field_type bz,
                int nx, int ng, bool mpi_only){
 field_bc_mpi(bx, nx, ng);
 field_bc_mpi(by, nx, ng);
 field_bc_mpi(bz, nx, ng);
}

void processor_summation_boundaries_mpi(field_type field, int nx, int ng){

    // Summation boundary
    for(int i = 0; i < ng; i++){
        field(nx+i) += field(i);
        field(ng+i) += field(nx+ng+i);
    }
}

void current_bcs(field_type jx, field_type jy, field_type jz,
                 int nx, int ng){
    processor_summation_boundaries_mpi(jx, nx, ng);
    processor_summation_boundaries_mpi(jy, nx, ng);
    processor_summation_boundaries_mpi(jz, nx, ng);
}

void current_finish(field_type jx, field_type jy, field_type jz,
                    int nx, int ng){
  current_bcs(jx, jy, jz, nx, ng);
  field_bc_mpi(jx, nx, ng);
  field_bc_mpi(jy, nx, ng);
  field_bc_mpi(jz, nx, ng);
}

void bfield_final_bcs(field_type bx, field_type by, field_type bz, int nx, int ng){
    bfield_bcs(bx, by, bz, nx, ng, true);
    bfield_bcs(bx, by, bz, nx, ng, false);
}
