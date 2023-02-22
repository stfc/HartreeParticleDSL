#ifndef _FDTD_MPI_IO_HDF5_HPP
#define _FDTD_MPI_IO_HDF5_HPP
#include "FDTD_MPI_field.hpp"

void store_domain_decomposition(struct FDTD_field &field, boundary &box);
void load_grid_hdf5(struct FDTD_field &field, char* filename,
                    int myrank, int nranks, boundary &box);

void grid_hdf5_output(struct FDTD_field field, char* filename, int myrank, int nranks);
#endif
