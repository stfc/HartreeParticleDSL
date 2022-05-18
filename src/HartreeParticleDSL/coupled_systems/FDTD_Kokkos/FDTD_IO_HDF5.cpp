#include <Kokkos_Core.hpp>
#include "FDTD_field.hpp"
#include "hdf5.h"
#include "unistd.h"

void grid_hdf5_output(struct FDTD_field field, char* filename){

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

    hid_t file_id;
    if(access(filename, F_OK) == 0){
        //Exists
        file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }else{
        //Doesn't exist
        file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    double* temp = (double*)malloc(sizeof(double) * field.nx);
    int ng = field.ng;
    int grid_size = field.nx;
    hsize_t h_dims[1];
    h_dims[0] = field.nx;
    hid_t dim = H5Screate_simple(1, h_dims, NULL);
    //Write the grid
    hid_t Electric_Field_Ex = H5Dcreate2(file_id, "Electric_Field_Ex", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ex(i+ng);
    }
    H5Dwrite(Electric_Field_Ex, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Electric_Field_Ex);

    hid_t Electric_Field_Ey = H5Dcreate2(file_id, "Electric_Field_Ey", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ey(i+ng);
    }
    H5Dwrite(Electric_Field_Ey, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Electric_Field_Ey);

    hid_t Electric_Field_Ez = H5Dcreate2(file_id, "Electric_Field_Ez", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = ez(i+ng);
    }
    H5Dwrite(Electric_Field_Ez, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Electric_Field_Ez);

    hid_t Magnetic_Field_Bx = H5Dcreate2(file_id, "Magnetic_Field_Bx", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = bx(i+ng);
    }
    H5Dwrite(Magnetic_Field_Bx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Magnetic_Field_Bx);

    hid_t Magnetic_Field_By = H5Dcreate2(file_id, "Magnetic_Field_By", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = by(i+ng);
    }
    H5Dwrite(Magnetic_Field_By, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Magnetic_Field_By);

    hid_t Magnetic_Field_Bz = H5Dcreate2(file_id, "Magnetic_Field_Bz", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = bz(i+ng);
    }
    H5Dwrite(Magnetic_Field_Bz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Magnetic_Field_Bz);

    hid_t Current_Jx = H5Dcreate2(file_id, "Current_Jx", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jx(i+ng);
    }
    H5Dwrite(Current_Jx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Current_Jx);

    hid_t Current_Jy = H5Dcreate2(file_id, "Current_Jy", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jy(i+ng);
    }
    H5Dwrite(Current_Jy, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Current_Jy);

    hid_t Current_Jz = H5Dcreate2(file_id, "Current_Jz", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i = 0; i < grid_size; i++){
        temp[i] = jz(i+ng);
    }
    H5Dwrite(Current_Jz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    H5Dclose(Current_Jz);

    free(temp);
    H5Fclose(file_id);

}
