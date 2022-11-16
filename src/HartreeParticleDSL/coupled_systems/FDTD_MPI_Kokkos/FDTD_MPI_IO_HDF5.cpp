#include <Kokkos_Core.hpp>
#include <Kokkos_Scatterview.hpp>
#include "mpi.h"
#include "FDTD_MPI_field.hpp"
#include "hdf5.h"
#include "unistd.h"

void store_domain_decomposition(struct FDTD_field &field, boundary &box){
    box.local_x_min = field.field.x_min_local;
    box.local_x_max = field.field.x_max_local;
    box.local_y_min = 0.0;
    box.local_y_max = 0.0;
    box.local_z_min = 0.0;
    box.local_z_max = 0.0;
}

void load_grid_hdf5(struct FDTD_field &field, char* filename,
         int myrank, int nranks, boundary &box){
    // Open HDF5 file
    hid_t file_id = H5Fopen("input.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Failed to open input.hdf5\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    hid_t temp_space;
    hsize_t dims[1];
// TODO

    hid_t f_ex = H5Dopen2(file_id, "Electric_Field_Ex", H5P_DEFAULT);
    if(f_ex < 0){
        printf("Failed to find Electric_Field_Ex in the input file, aborting.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // Load the grid size
    temp_space = H5Dget_space(f_ex);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    int size_grid = dims[0];

    field.nxglobal = size_grid;
    field.dx = (box.x_max - box.x_min) / ((double)field.nxglobal);
    field.nx = 0;

    if(field.nxglobal % nranks == 0){
        field.nx = field.nxglobal / nranks;
    }else{
        field.nx = field.nxglobal / nranks;
        if(field.nxglobal % nranks > myrank){
            (field.nx)++;
        }
    }
    int min_local_cell = (myrank * (field.nxglobal / nranks));
    if(field.nxglobal % nranks != 0){
        if( myrank > field.nxglobal % nranks ){
            field.min_local_cell += (field.nxglobal % nranks);
        }else{
            field.min_local_cell += (myrank);
        }
    }
    //Exclusive
    field.max_local_cell = field.min_local_cell + field.nxlocal;

    //Now we got the local grid info, load the local grid.
    size_grid = field.nxlocal;
    field.ex = field_type("ex", size_grid + 2*field.ng);
    field.ey = field_type("ey", size_grid + 2*field.ng);
    field.ez = field_type("ez", size_grid + 2*field.ng);
    field.bx = field_type("bx", size_grid + 2*field.ng);
    field.by = field_type("by", size_grid + 2*field.ng);
    field.bz = field_type("bz", size_grid + 2*field.ng);
    field.jx = field_type("jx", size_grid + 2*field.jng);
    field.jy = field_type("jy", size_grid + 2*field.jng);
    field.jz = field_type("jz", size_grid + 2*field.jng);

    field.scatter_jx = scatter_field_type(field.jx);
    field.scatter_jy = scatter_field_type(field.jy);
    field.scatter_jz = scatter_field_type(field.jz);

    auto host_ex = Kokkos::create_mirror_view(field.ex);
    auto host_ey = Kokkos::create_mirror_view(field.ey);
    auto host_ez = Kokkos::create_mirror_view(field.ez);
    auto host_bx = Kokkos::create_mirror_view(field.bx);
    auto host_by = Kokkos::create_mirror_view(field.by);
    auto host_bz = Kokkos::create_mirror_view(field.bz);
    auto host_jx = Kokkos::create_mirror_view(field.jx);
    auto host_jy = Kokkos::create_mirror_view(field.jy);
    auto host_jz = Kokkos::create_mirror_view(field.jz);

    double* ex_temp_array = (double*) malloc(sizeof(double) * field.nxglobal);
    //already done hid_t f_ex = H5Dopen2(file_id, "Electric_Field_Ex", H5P_DEFAULT);
    H5Dread(f_ex, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ex_temp_array);
    for(int i = 0; i < size_grid; i++){
        host_ex(i+field.ng) = ex_temp_array[min_local_cell + i];
    }
    free(ex_temp_array);
    H5Dclose(f_ex);

    hid_t f_ey = H5Dopen2(file_id, "Electric_Field_Ey", H5P_DEFAULT);
    if(f_ey >= 0){
        double* ey_temp_array = (double*) malloc(sizeof(double) * field.nxglobal);
        H5Dread(f_ey, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ey_temp_array);
        for(int i = 0; i < size_grid; i++){
            host_ey(i+field.ng) = ey_temp_array[min_local_cell + i];
        }
        free(ey_temp_array);
        H5Dclose(f_ey);
    }else{
        // If we don't find a field we assume its 0
        for(int i = 0; i < size_grid; i++){
            host_ey(i+field.ng) = 0.0;
        }
    }

    hid_t f_ez = H5Dopen2(file_id, "Electric_Field_Ez", H5P_DEFAULT);
    if(f_ez >= 0){
        double* ez_temp_array = (double*) malloc(sizeof(double) * field.nxglobal);
        H5Dread(f_ez, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ez_temp_array);
        for(int i = 0; i < size_grid; i++){
            host_ez(i+field.ng) = ez_temp_array[min_local_cell + i];
        }
        free(ez_temp_array);
        H5Dclose(f_ez);
    }else{
        for(int i = 0; i < size_grid; i++){
            host_ez(i+field.ng) = 0.0;
        }
    }

    hid_t f_bx = H5Dopen2(file_id, "Magnetic_Field_Bx", H5P_DEFAULT);
    if(f_bx >= 0){
        double* bx_temp_array = (double*) malloc(sizeof(double) * *nxglobal);
        H5Dread(f_bx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bx_temp_array);
        for(int i = 0; i < size_grid; i++){
            host_bx(i+field.ng) = bx_temp_array[min_local_cell + i];
        }
        free(bx_temp_array);
        H5Dclose(f_bx);
    }else{
        for(int i = 0; i < size_grid; i++){
            host_bx(i+field.ng) = 0.0;
        }
    }

    hid_t f_by = H5Dopen2(file_id, "Magnetic_Field_By", H5P_DEFAULT);
    if(f_by >= 0){
        double* by_temp_array = (double*) malloc(sizeof(double) * *nxglobal);
        H5Dread(f_by, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, by_temp_array);
        for(int i = 0; i < size_grid; i++){
            host_by(i+field.ng) = by_temp_array[min_local_cell + i];
        }
        free(by_temp_array);
        H5Dclose(f_by);
    }else{
        for(int i = 0; i < size_grid; i++){
            host_by(i+field.ng) = 0.0;
        }
    }

    hid_t f_bz = H5Dopen2(file_id, "Magnetic_Field_Bz", H5P_DEFAULT);
    if(f_bz >= 0){
        double* bz_temp_array = (double*) malloc(sizeof(double) * *nxglobal);
        H5Dread(f_bz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bz_temp_array);
        printf("Read %i elements into host_bz, size of kokkos ds is %i\n", size_grid, host_bz.size());
        for(int i = 0; i < size_grid; i++){
            host_bz(i+field.ng) = bz_temp_array[min_local_cell + i];
        }
        free(bz_temp_array);
        H5Dclose(f_bz);
    }else{
        for(int i = 0; i < size_grid; i++){
            host_bz(i+field.ng) = 0.0;
        }
    }

    for(int i = 0; i < field.ng; i++){
        host_ex(i) = 0.0;
        host_ey(i) = 0.0;
        host_ez(i) = 0.0;
        host_bx(i) = 0.0;
        host_by(i) = 0.0;
        host_bz(i) = 0.0;
        host_ex(field.ng+size_grid+i) = 0.0;
        host_ey(field.ng+size_grid+i) = 0.0;
        host_ez(field.ng+size_grid+i) = 0.0;
        host_bx(field.ng+size_grid+i) = 0.0;
        host_by(field.ng+size_grid+i) = 0.0;
        host_bz(field.ng+size_grid+i) = 0.0;
    }
    for(int i = 0; i < size_grid+2*field.jng; i++){
        host_jx(i) = 0.0;
        host_jy(i) = 0.0;
        host_jz(i) = 0.0;
    }

    field.field.x_min_local = field.min_local_cell * field.dx;
    field.field.x_max_local = field.max_local_cell * field.dx;
    field.field.x_grid_min_local = field.field.x_min_local + field.dx/2.0;
    field.field.x_grid_max_local = field.field.x_max_local - field.dx/2.0;
}


void grid_hdf5_output(struct FDTD_field field, char* filename, int myrank, int nranks){
// TODO
    auto ex = Kokkos::create_mirror_view(field.ex);
    auto ey = Kokkos::create_mirror_view(field.ey);
    auto ez = Kokkos::create_mirror_view(field.ez);
    auto bx = Kokkos::create_mirror_view(field.bx);
    auto by = Kokkos::create_mirror_view(field.by);
    auto bz = Kokkos::create_mirror_view(field.bz);
    auto jx = Kokkos::create_mirror_view(field.jx);
    auto jy = Kokkos::create_mirror_view(field.jy);
    auto jz = Kokkos::create_mirror_view(field.jz);

    Kokkos::deep_copy(ex, field.ex);
    Kokkos::deep_copy(ey, field.ey);
    Kokkos::deep_copy(ez, field.ez);
    Kokkos::deep_copy(bx, field.bx);
    Kokkos::deep_copy(by, field.by);
    Kokkos::deep_copy(bz, field.bz);
    Kokkos::deep_copy(jx, field.jx);
    Kokkos::deep_copy(jy, field.jy);
    Kokkos::deep_copy(jz, field.jz);

    // Setup parallel access
    hid_t acc_template = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info; MPI_Info_create(&info);
    H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
    H5Pset_coll_metadata_write(acc_template, 1);
    hid_t file_id;
    if(access(filename, F_OK) == 0){
        //Exists
        file_id = H5Fopen(filename, H5F_ACC_RDWR, acc_template);
    }else{
        //Doesn't exist
        file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    }

    hid_t xf_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_COLLECTIVE);

    int grid_size = field.nx;
    double* temp = (double*)malloc(sizeof(double) * grid_size);
    int ng = field.ng;
    int grid_size = field.nx;
    hsize_t h_dims[1];
    h_dims[0] = grid_size;
    hsize_t nxg = field.nxglobal;
    global_dim = H5Screate_simple(1, &nxg, NULL);

    int my_offset = 0;
    if(myrank == 0 && nranks > 1){
        int nxlocal = grid_size;
        MPI_Send(&nxlocal, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }else if (myrank == nranks-1){
        MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }else{
        MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int nxlocal = my_offset + grid_size;
        MPI_Send(&nxlocal, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD);
    }

    hsize_t offset[2];
    hsize_t count[2];
    offset[0] = my_offset;
    offset[1] = 0;
    count[0] = h_dims[0];
    count[1] = 0;
    H5Sselect_hyperslab(global_dim, H5S_SELECT_SET, offset, NULL, count, NULL);

    h_dims[0] = grid_size;

    memspace = H5Screate_simple(1, h_dims, NULL);

    //Offset in memspace is 0
    offset[0] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    //Write the grid
    hid_t Electric_Field_Ex = H5Dcreate2(file_id, "Electric_Field_Ex", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ex(i+ng);
    }
    H5Dwrite(Electric_Field_Ex, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Electric_Field_Ex);

    hid_t Electric_Field_Ey = H5Dcreate2(file_id, "Electric_Field_Ey", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ey(i+ng);
    }
    H5Dwrite(Electric_Field_Ey, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Electric_Field_Ey);

    hid_t Electric_Field_Ez = H5Dcreate2(file_id, "Electric_Field_Ez", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ez(i+ng);
    }
    H5Dwrite(Electric_Field_Ez, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Electric_Field_Ez);

    hid_t Magnetic_Field_Bx = H5Dcreate2(file_id, "Magnetic_Field_Bx", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = bx(i+ng);
    }
    H5Dwrite(Magnetic_Field_Bx, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Magnetic_Field_Bx);

    hid_t Magnetic_Field_By = H5Dcreate2(file_id, "Magnetic_Field_By", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = by(i+ng);
    }
    H5Dwrite(Magnetic_Field_By, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Magnetic_Field_By);

    hid_t Magnetic_Field_Bz = H5Dcreate2(file_id, "Magnetic_Field_Bz", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = bz(i+ng);
    }
    H5Dwrite(Magnetic_Field_Bz, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Magnetic_Field_Bz);

    hid_t Current_Jx = H5Dcreate2(file_id, "Current_Jx", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jx(i+ng);
    }
    H5Dwrite(Current_Jx, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Current_Jx);

    hid_t Current_Jy = H5Dcreate2(file_id, "Current_Jy", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jy(i+ng);
    }
    H5Dwrite(Current_Jy, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Current_Jy);

    hid_t Current_Jz = H5Dcreate2(file_id, "Current_Jz", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jz(i+ng);
    }
    H5Dwrite(Current_Jz, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, temp);
    H5Dclose(Current_Jz);

    free(temp);
    H5Fclose(file_id);
}
