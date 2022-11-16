#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "mpi.h"
#include "FDTD_MPI_boundaries.hpp"

void field_bc_mpi(field_type field, int nx, int ng){

    int myrank = 0; MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    int nranks = 1; MPI_Comm_size( MPI_COMM_WORLD, &nranks );
    int previous_rank = ( myrank == 0 ) ? nranks - 1 : myrank - 1;
    int next_rank = ( myrank == nranks - 1 ) ? 0 : myrank + 1;
    double* bottom = (double*) malloc(sizeof(double) * ng);
    double* top = (double*) malloc(sizeof(double) * ng);

    MPI_Request send_bottom_request, recv_bottom_request;
    MPI_Request send_top_request, recv_top_request;
    MPI_Isend( &(field.data()[ng]), ng, MPI_DOUBLE, previous_rank, 1, MPI_COMM_WORLD, &send_bottom_request); //Send bottom
    MPI_Isend( &(field.data()[nx]), ng, MPI_DOUBLE, next_rank, 2, MPI_COMM_WORLD, &send_top_request); //Send top
    MPI_Irecv( top, ng, MPI_DOUBLE, next_rank, 1, MPI_COMM_WORLD, &recv_top_request); //Recv top
    MPI_Irecv( bottom, ng, MPI_DOUBLE, previous_rank, 2, MPI_COMM_WORLD, &recv_bottom_request); //Recv bottom

    // Do the communication
    MPI_Request arrayreq[4];
    arrayreq[0] = send_bottom_request;
    arrayreq[1] = recv_bottom_request;
    arrayreq[2] = send_top_request;
    arrayreq[3] = recv_top_request;
    MPI_Waitall(4, arrayreq, MPI_STATUSES_IGNORE);

    for(int i = 0; i < ng; i++){
        field(ng+nx+i) = top[i];
        field(i) = bottom[i];
    }

    free(bottom);
    free(top);
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

    int myrank = 0; MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    int nranks = 1; MPI_Comm_size( MPI_COMM_WORLD, &nranks );
    int previous_rank = ( myrank == 0 ) ? nranks - 1 : myrank - 1;
    int next_rank = ( myrank == nranks - 1 ) ? 0 : myrank + 1;
    double* bottom = (double*) malloc(sizeof(double) * ng);
    double* top = (double*) malloc(sizeof(double) * ng);

    MPI_Request send_bottom_request, recv_bottom_request;
    MPI_Request send_top_request, recv_top_request;
    MPI_Isend( &(field.data()[0]), ng, MPI_DOUBLE, previous_rank, 1, MPI_COMM_WORLD, &send_bottom_request); //Send bottom
    MPI_Isend( &(field.data()[nx+ng]), ng, MPI_DOUBLE, next_rank, 2, MPI_COMM_WORLD, &send_top_request); //Send top
    MPI_Irecv( top, ng, MPI_DOUBLE, next_rank, 1, MPI_COMM_WORLD, &recv_top_request); //Recv top
    MPI_Irecv( bottom, ng, MPI_DOUBLE, previous_rank, 2, MPI_COMM_WORLD, &recv_bottom_request); //Recv bottom

    // Do the communication
    MPI_Request arrayreq[4];
    arrayreq[0] = send_bottom_request;
    arrayreq[1] = recv_bottom_request;
    arrayreq[2] = send_top_request;
    arrayreq[3] = recv_top_request;
    MPI_Waitall(4, arrayreq, MPI_STATUSES_IGNORE);

    for(int i = 0; i < ng; i++){
        field(nx+i) += top[i];
        field(ng+i) += bottom[i];
    }
    free(bottom);
    free(top);

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
