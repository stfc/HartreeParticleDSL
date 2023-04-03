import HartreeParticleDSL.IO_modules.PHDF5_IO.PHDF5_IO as PHDF5_IO
import pytest
from HartreeParticleDSL.HartreeParticleDSL import Particle, set_mpi, reset_for_tests
from HartreeParticleDSL.HartreeParticleDSLExceptions import UnsupportedCodeError

def test_phdf5_init():
    a = PHDF5_IO.PHDF5_IO(indent=8)
    assert a._inputs == {}
    assert a._outputs == {}
    assert a._indent == 8
    assert a._current_indent == 0
    a = PHDF5_IO.PHDF5_IO()
    assert a._indent == 4
    assert a._current_indent == 0

def test_phdf5_add_input():
    a = PHDF5_IO.PHDF5_IO()
    z = "h_in"
    y = "part_in"
    a.add_input(z, y)
    assert a._inputs[z] == y

def test_phdf5_add_output():
    a = PHDF5_IO.PHDF5_IO()
    z = "h_out"
    y = "part_out"
    a.add_output(z, y)
    assert a._outputs[z] == y

def test_phdf5_indent():
    a = PHDF5_IO.PHDF5_IO(indent=1)
    assert a.indent() == ""

def test_phdf5_increment_indent():
    a = PHDF5_IO.PHDF5_IO(indent=1)
    a._increment_indent()
    assert a.indent() == " "

def test_phdf5_decrement_indent():
    a = PHDF5_IO.PHDF5_IO(indent=1)
    a._increment_indent()
    assert a.indent() == " "
    a._decrement_indent()
    assert a.indent() == ""

def test_phdf5_get_header_includes_cabana_pir():
    set_mpi(True)
    a = PHDF5_IO.PHDF5_IO()
    ins = a.get_header_includes_cabana_pir()
    assert "\"hdf5.h\"" in ins
    assert "\"mpi.h\"" in ins
    set_mpi(False)
    ins = a.get_header_includes_cabana_pir()
    assert "\"mpi.h\"" not in ins
    reset_for_tests()

def test_phdf5_get_includes_cabana_pir():
    set_mpi(True)
    a = PHDF5_IO.PHDF5_IO()
    ins = a.get_includes_cabana_pir()
    assert "\"hdf5.h\"" in ins
    assert "\"mpi.h\"" in ins
    set_mpi(False)
    ins = a.get_header_includes_cabana_pir()
    assert "\"mpi.h\"" not in ins
    reset_for_tests()

def test_phdf5_gen_code_cabana_pir():
    a = PHDF5_IO.PHDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_cabana_pir(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")

    with pytest.raises(UnsupportedCodeError) as excinfo:
        code = a.gen_code_cabana_pir(part)
    assert("thing element not supported in Parallel HDF5 IO." in str(excinfo.value))

    set_mpi(True)
    a = PHDF5_IO.PHDF5_IO()
    part = Particle()
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    out = a.gen_code_cabana_pir(part)
    correct = '''void get_box_size(boundary &box, const char* filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        std::cout << "Failed to open file " << filename << "\\n";
        exit(1);
    }
    hid_t temp_space;
    hsize_t dims[1];

    hid_t boxsize = H5Dopen2(file_id, "Box_Size", H5P_DEFAULT);
    temp_space = H5Dget_space(boxsize);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    int size_box = dims[0];
    double* box_temp = (double*) malloc(sizeof(double) * size_box);
    H5Dread(boxsize, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_temp);
    box.x_min = box_temp[0];
    box.x_max = box_temp[1];
    if(size_box > 2){
        box.y_min = box_temp[2];
        box.y_max = box_temp[3];
    }else{
        box.y_min = 0.0;
        box.y_max = 0.0;
    }
    if(size_box > 4){
        box.z_min = box_temp[4];
        box.z_max = box_temp[5];
    }else{
        box.z_min = 0.0;
        box.z_max = 0.0;
    }

    free(box_temp);
    H5Dclose(boxsize);
    H5Fclose(file_id);
}

template <class aosoa_class, class aosoa_host_class> void hdf5_input(aosoa_class &non_host_aosoa,
    aosoa_host_class &particle_aosoa, config_type &config, boundary &box, const char* filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        std::cout << "Failed to open file " << filename << "\\n";
        exit(1);
    }
    hid_t temp_space;
    hsize_t dims[1];

    hid_t pos_x_test_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    if( pos_x_test_var < 0){
        printf("Failed to find dataset pos_x\\n");
        exit(1);
    }
    H5Dclose(pos_x_test_var);
    hid_t pos_y_test_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    if( pos_y_test_var < 0){
        printf("Failed to find dataset pos_y\\n");
        exit(1);
    }
    H5Dclose(pos_y_test_var);
    hid_t pos_z_test_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    if( pos_z_test_var < 0){
        printf("Failed to find dataset pos_z\\n");
        exit(1);
    }
    H5Dclose(pos_z_test_var);
    hid_t vel_x_test_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    if( vel_x_test_var < 0){
        printf("Failed to find dataset vel_x\\n");
        exit(1);
    }
    H5Dclose(vel_x_test_var);
    hid_t vel_y_test_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    if( vel_y_test_var < 0){
        printf("Failed to find dataset vel_y\\n");
        exit(1);
    }
    H5Dclose(vel_y_test_var);
    hid_t vel_z_test_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    if( vel_z_test_var < 0){
        printf("Failed to find dataset vel_z\\n");
        exit(1);
    }
    H5Dclose(vel_z_test_var);
    hid_t part_count_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    hid_t space = H5Dget_space(part_count_var);
    int ndims = H5Sget_simple_extent_ndims(space);
    hid_t pos_x_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_x_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_x_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_x_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_x_temp_array);
    H5Dclose(pos_x_var);
    hid_t pos_y_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_y_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_y_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_y_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_y_temp_array);
    H5Dclose(pos_y_var);
    hid_t pos_z_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_z_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_z_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_z_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_z_temp_array);
    H5Dclose(pos_z_var);
    int global_parts = dims[0];
    int num_parts = 0;
    for(int i = 0; i < global_parts; i++){
        if((pos_x_temp_array[i] < box.local_x_max && pos_x_temp_array[i] >= box.local_x_min) && (pos_y_temp_array[i] < box.local_y_max && pos_y_temp_array[i] >= box.local_y_min) && (pos_z_temp_array[i] < box.local_z_max && pos_z_temp_array[i] >= box.local_z_min)){
            num_parts++;
        }
    }
    config.config_host(0).space.nparts = global_parts;
    int new_size = static_cast<int>(num_parts);
    particle_aosoa.resize(new_size);
    non_host_aosoa.resize(new_size);
    int counter = 0;
    auto pos_slice = Cabana::slice<core_part_position>(particle_aosoa);
    for(int i = 0; i < global_parts; i++){
        if((pos_x_temp_array[i] < box.local_x_max && pos_x_temp_array[i] >= box.local_x_min) && (pos_y_temp_array[i] < box.local_y_max && pos_y_temp_array[i] >= box.local_y_min) && (pos_z_temp_array[i] < box.local_z_max && pos_z_temp_array[i] >= box.local_z_min)){
            pos_slice(counter, 0) = pos_x_temp_array[i];
            pos_slice(counter, 1) = pos_y_temp_array[i];
            pos_slice(counter, 2) = pos_z_temp_array[i];
            counter++;
        }
    }
    hid_t vel_x_read_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    space = H5Dget_space(vel_x_read_var);
    auto vel_x_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_x_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_x_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_x_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        if((pos_x_temp_array[i] < box.local_x_max && pos_x_temp_array[i] >= box.local_x_min) && (pos_y_temp_array[i] < box.local_y_max && pos_y_temp_array[i] >= box.local_y_min) && (pos_z_temp_array[i] < box.local_z_max && pos_z_temp_array[i] >= box.local_z_min)){
            vel_x_slice(counter, 0) = vel_x_temp_array[i];
            counter++;
        }
    }
    free(vel_x_temp_array);
    H5Dclose(vel_x_read_var);
    hid_t vel_y_read_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    space = H5Dget_space(vel_y_read_var);
    auto vel_y_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_y_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_y_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_y_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        if((pos_x_temp_array[i] < box.local_x_max && pos_x_temp_array[i] >= box.local_x_min) && (pos_y_temp_array[i] < box.local_y_max && pos_y_temp_array[i] >= box.local_y_min) && (pos_z_temp_array[i] < box.local_z_max && pos_z_temp_array[i] >= box.local_z_min)){
            vel_y_slice(counter, 1) = vel_y_temp_array[i];
            counter++;
        }
    }
    free(vel_y_temp_array);
    H5Dclose(vel_y_read_var);
    hid_t vel_z_read_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    space = H5Dget_space(vel_z_read_var);
    auto vel_z_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_z_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_z_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_z_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        if((pos_x_temp_array[i] < box.local_x_max && pos_x_temp_array[i] >= box.local_x_min) && (pos_y_temp_array[i] < box.local_y_max && pos_y_temp_array[i] >= box.local_y_min) && (pos_z_temp_array[i] < box.local_z_max && pos_z_temp_array[i] >= box.local_z_min)){
            vel_z_slice(counter, 2) = vel_z_temp_array[i];
            counter++;
        }
    }
    free(vel_z_temp_array);
    H5Dclose(vel_z_read_var);
    free(pos_x_temp_array);
    free(pos_y_temp_array);
    free(pos_z_temp_array);
    H5Fclose(file_id);
    Cabana::deep_copy(non_host_aosoa, particle_aosoa);
}'''
    assert correct in out
    set_mpi(False)

    out = a.gen_code_cabana_pir(part)
    correct = '''void get_box_size(boundary &box, const char* filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        std::cout << "Failed to open file " << filename << "\\n";
        exit(1);
    }
    hid_t temp_space;
    hsize_t dims[1];

    hid_t boxsize = H5Dopen2(file_id, "Box_Size", H5P_DEFAULT);
    temp_space = H5Dget_space(boxsize);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    int size_box = dims[0];
    double* box_temp = (double*) malloc(sizeof(double) * size_box);
    H5Dread(boxsize, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_temp);
    box.x_min = box_temp[0];
    box.x_max = box_temp[1];
    if(size_box > 2){
        box.y_min = box_temp[2];
        box.y_max = box_temp[3];
    }else{
        box.y_min = 0.0;
        box.y_max = 0.0;
    }
    if(size_box > 4){
        box.z_min = box_temp[4];
        box.z_max = box_temp[5];
    }else{
        box.z_min = 0.0;
        box.z_max = 0.0;
    }

    free(box_temp);
    H5Dclose(boxsize);
    H5Fclose(file_id);
}

template <class aosoa_class, class aosoa_host_class> void hdf5_input(aosoa_class &non_host_aosoa,
    aosoa_host_class &particle_aosoa, config_type &config, boundary &box, const char* filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        std::cout << "Failed to open file " << filename << "\\n";
        exit(1);
    }
    hid_t temp_space;
    hsize_t dims[1];

    hid_t pos_x_test_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    if( pos_x_test_var < 0){
        printf("Failed to find dataset pos_x\\n");
        exit(1);
    }
    H5Dclose(pos_x_test_var);
    hid_t pos_y_test_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    if( pos_y_test_var < 0){
        printf("Failed to find dataset pos_y\\n");
        exit(1);
    }
    H5Dclose(pos_y_test_var);
    hid_t pos_z_test_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    if( pos_z_test_var < 0){
        printf("Failed to find dataset pos_z\\n");
        exit(1);
    }
    H5Dclose(pos_z_test_var);
    hid_t vel_x_test_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    if( vel_x_test_var < 0){
        printf("Failed to find dataset vel_x\\n");
        exit(1);
    }
    H5Dclose(vel_x_test_var);
    hid_t vel_y_test_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    if( vel_y_test_var < 0){
        printf("Failed to find dataset vel_y\\n");
        exit(1);
    }
    H5Dclose(vel_y_test_var);
    hid_t vel_z_test_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    if( vel_z_test_var < 0){
        printf("Failed to find dataset vel_z\\n");
        exit(1);
    }
    H5Dclose(vel_z_test_var);
    hid_t part_count_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    hid_t space = H5Dget_space(part_count_var);
    int ndims = H5Sget_simple_extent_ndims(space);
    hid_t pos_x_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_x_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_x_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_x_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_x_temp_array);
    H5Dclose(pos_x_var);
    hid_t pos_y_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_y_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_y_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_y_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_y_temp_array);
    H5Dclose(pos_y_var);
    hid_t pos_z_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    temp_space = H5Dget_space(pos_z_var);
    H5Sget_simple_extent_dims(temp_space, dims, NULL);
    double* pos_z_temp_array = (double*) malloc(sizeof(double) * dims[0]);
    H5Dread(pos_z_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_z_temp_array);
    H5Dclose(pos_z_var);
    int global_parts = dims[0];
    int num_parts = 0;
    num_parts = global_parts;
    config.config_host(0).space.nparts = global_parts;
    int new_size = static_cast<int>(num_parts);
    particle_aosoa.resize(new_size);
    non_host_aosoa.resize(new_size);
    int counter = 0;
    auto pos_slice = Cabana::slice<core_part_position>(particle_aosoa);
    for(int i = 0; i < global_parts; i++){
        pos_slice(counter, 0) = pos_x_temp_array[i];
        pos_slice(counter, 1) = pos_y_temp_array[i];
        pos_slice(counter, 2) = pos_z_temp_array[i];
        counter++;
    }
    hid_t vel_x_read_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    space = H5Dget_space(vel_x_read_var);
    auto vel_x_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_x_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_x_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_x_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        vel_x_slice(counter, 0) = vel_x_temp_array[i];
        counter++;
    }
    free(vel_x_temp_array);
    H5Dclose(vel_x_read_var);
    hid_t vel_y_read_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    space = H5Dget_space(vel_y_read_var);
    auto vel_y_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_y_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_y_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_y_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        vel_y_slice(counter, 1) = vel_y_temp_array[i];
        counter++;
    }
    free(vel_y_temp_array);
    H5Dclose(vel_y_read_var);
    hid_t vel_z_read_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    space = H5Dget_space(vel_z_read_var);
    auto vel_z_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    double* vel_z_temp_array = (double*) malloc(sizeof(double) * global_parts);
    H5Dread(vel_z_read_var, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, vel_z_temp_array);
    counter = 0;
    for(int i = 0; i < global_parts; i++){
        vel_z_slice(counter, 2) = vel_z_temp_array[i];
        counter++;
    }
    free(vel_z_temp_array);
    H5Dclose(vel_z_read_var);
    free(pos_x_temp_array);
    free(pos_y_temp_array);
    free(pos_z_temp_array);
    H5Fclose(file_id);
    Cabana::deep_copy(non_host_aosoa, particle_aosoa);
}'''

    assert correct in out

    a = PHDF5_IO.PHDF5_IO()
    part = Particle()
    a.add_output("pos_x", "core_part.position[0]")
    a.add_output("pos_y", "core_part.position.y")
    a.add_output("pos_z", "core_part.position.z")
    a.add_output("vel_x", "core_part.velocity[0]")
    a.add_output("vel_y", "core_part.velocity[1]")
    a.add_output("vel_z", "core_part.velocity[2]")
    set_mpi(True)
    out = a.gen_code_cabana_pir(part)
    correct = '''template <class aosoa_host_class> void hdf5_output(aosoa_host_class &particle_aosoa,
    const char* filename, boundary &box, config_type &config, int myrank, int nranks){
    hid_t acc_template = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info; MPI_Info_create(&info);
    H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
    H5Pset_coll_metadata_write(acc_template, 1); //metadata writes are collective
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    if(file_id < 0){
        std::cout << "[" << myrank << "] failed to open " << filename << "\\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    hsize_t h_dims[1];
    h_dims[0] = particle_aosoa.size();
    int my_offset = 0;
    if(myrank == 0 && nranks > 1){
        int npart = particle_aosoa.size();        MPI_Send(&npart, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }else if (myrank == nranks-1){
        MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }else{
        MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int npart = my_offset + particle_aosoa.size();
        MPI_Send(&npart, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD);
    }
    hid_t xf_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_COLLECTIVE);
    int global_part_count = 0;
    MPI_Allreduce( &h_dims[0], &global_part_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    hsize_t gpc = (hsize_t) global_part_count;
    hid_t global_dim = H5Screate_simple(1, &gpc, NULL);
    hsize_t offset[2];
    hsize_t count[2];
    offset[0] = my_offset;
    offset[1] = 0;
    count[0] = h_dims[0];
    count[1] = 0;
    H5Sselect_hyperslab(global_dim, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t memspace = H5Screate_simple(1, h_dims, NULL);
    //Offset in memspace is 0
    offset[0] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t pos_x_output_field = H5Dcreate2(file_id, "pos_x", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_x_output_array[i] = Cabana::get<core_part_position>(part, 0);
    }

    H5Dwrite(pos_x_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_x_output_array);
    free(pos_x_output_array);
    H5Dclose(pos_x_output_field);
    hid_t pos_y_output_field = H5Dcreate2(file_id, "pos_y", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_y_output_array[i] = Cabana::get<core_part_position>(part, 1);
    }

    H5Dwrite(pos_y_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_y_output_array);
    free(pos_y_output_array);
    H5Dclose(pos_y_output_field);
    hid_t pos_z_output_field = H5Dcreate2(file_id, "pos_z", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_z_output_array[i] = Cabana::get<core_part_position>(part, 2);
    }

    H5Dwrite(pos_z_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_z_output_array);
    free(pos_z_output_array);
    H5Dclose(pos_z_output_field);
    hid_t vel_x_output_field = H5Dcreate2(file_id, "vel_x", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_x_output_array[i] = Cabana::get<core_part_velocity>(part, 0);
    }

    H5Dwrite(vel_x_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_x_output_array);
    free(vel_x_output_array);
    H5Dclose(vel_x_output_field);
    hid_t vel_y_output_field = H5Dcreate2(file_id, "vel_y", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_y_output_array[i] = Cabana::get<core_part_velocity>(part, 1);
    }

    H5Dwrite(vel_y_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_y_output_array);
    free(vel_y_output_array);
    H5Dclose(vel_y_output_field);
    hid_t vel_z_output_field = H5Dcreate2(file_id, "vel_z", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_z_output_array[i] = Cabana::get<core_part_velocity>(part, 2);
    }

    H5Dwrite(vel_z_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_z_output_array);
    free(vel_z_output_array);
    H5Dclose(vel_z_output_field);


    H5Fclose(file_id);
}'''
    assert correct in out
    set_mpi(False)
    out = a.gen_code_cabana_pir(part)
    correct = '''template <class aosoa_host_class> void hdf5_output(aosoa_host_class &particle_aosoa,
    const char* filename, boundary &box, config_type &config, int myrank, int nranks){
    hid_t acc_template = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    if(file_id < 0){
        std::cout << "[" << myrank << "] failed to open " << filename << "\\n";
        abort();
    }

    hsize_t h_dims[1];
    h_dims[0] = particle_aosoa.size();
    int my_offset = 0;
    hid_t xf_id = H5P_DEFAULT;
    int global_part_count = 0;
    global_part_count = h_dims[0];
    hsize_t gpc = (hsize_t) global_part_count;
    hid_t global_dim = H5Screate_simple(1, &gpc, NULL);
    hsize_t offset[2];
    hsize_t count[2];
    offset[0] = my_offset;
    offset[1] = 0;
    count[0] = h_dims[0];
    count[1] = 0;
    H5Sselect_hyperslab(global_dim, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t memspace = H5Screate_simple(1, h_dims, NULL);
    //Offset in memspace is 0
    offset[0] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t pos_x_output_field = H5Dcreate2(file_id, "pos_x", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_x_output_array[i] = Cabana::get<core_part_position>(part, 0);
    }

    H5Dwrite(pos_x_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_x_output_array);
    free(pos_x_output_array);
    H5Dclose(pos_x_output_field);
    hid_t pos_y_output_field = H5Dcreate2(file_id, "pos_y", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_y_output_array[i] = Cabana::get<core_part_position>(part, 1);
    }

    H5Dwrite(pos_y_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_y_output_array);
    free(pos_y_output_array);
    H5Dclose(pos_y_output_field);
    hid_t pos_z_output_field = H5Dcreate2(file_id, "pos_z", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_z_output_array[i] = Cabana::get<core_part_position>(part, 2);
    }

    H5Dwrite(pos_z_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, pos_z_output_array);
    free(pos_z_output_array);
    H5Dclose(pos_z_output_field);
    hid_t vel_x_output_field = H5Dcreate2(file_id, "vel_x", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_x_output_array[i] = Cabana::get<core_part_velocity>(part, 0);
    }

    H5Dwrite(vel_x_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_x_output_array);
    free(vel_x_output_array);
    H5Dclose(vel_x_output_field);
    hid_t vel_y_output_field = H5Dcreate2(file_id, "vel_y", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_y_output_array[i] = Cabana::get<core_part_velocity>(part, 1);
    }

    H5Dwrite(vel_y_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_y_output_array);
    free(vel_y_output_array);
    H5Dclose(vel_y_output_field);
    hid_t vel_z_output_field = H5Dcreate2(file_id, "vel_z", H5T_NATIVE_DOUBLE, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_z_output_array[i] = Cabana::get<core_part_velocity>(part, 2);
    }

    H5Dwrite(vel_z_output_field, H5T_NATIVE_DOUBLE, memspace, global_dim, xf_id, vel_z_output_array);
    free(vel_z_output_array);
    H5Dclose(vel_z_output_field);


    H5Fclose(file_id);
}''' 
    assert correct in out

    part.add_element("thing", "c_int")
    a.add_output("thing", "thing")
    with pytest.raises(UnsupportedCodeError) as excinfo:
        out = a.gen_code_cabana_pir(part)
    assert("thing element not supported in HDF5 IO." in str(excinfo.value))


def test_phdf5_get_box_size_cabana_pir():
    a = PHDF5_IO.PHDF5_IO()
    out = a.call_get_box_size_pir(123, "a")
    assert "    get_box_size(config.config_host(0).space.box_dims, a);\n" == out

def test_phdf5_call_input_cabana_pir():
    a = PHDF5_IO.PHDF5_IO()
    out = a.call_input_cabana_pir(123, "a")
    correct = '''/* Create structures of size 1 to initialise, the HDF5 function will resize them */
    Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 1);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 1);
    hdf5_input<decltype(particle_aosoa), decltype(particle_aosoa_host)>(particle_aosoa, particle_aosoa_host, config, config.config_host(0).space.box_dims, a);
'''
    assert correct in out

def test_phdf5_call_output_cabana_pir():
    set_mpi(True)
    a = PHDF5_IO.PHDF5_IO()
    out = a.call_output_cabana_pir(123, "\"a\"", "varname")
    correct = '''{
        char filename[300];
        sprintf(filename, "a%.4d.hdf5", varname);
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &nranks );
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out
    set_mpi(False)
    out = a.call_output_cabana_pir(123, "\"a\"", "varname")
    correct = '''{
        char filename[300];
        sprintf(filename, "a%.4d.hdf5", varname);
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        myrank = 0;
        nranks = 1;
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out
    set_mpi(True)
    out = a.call_output_cabana_pir(123, "a", "varname")
    correct = '''{
        char filename[300];
        sprintf(filename, "a%.4d.hdf5", varname);
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &nranks );
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out
    set_mpi(False)
    out = a.call_output_cabana_pir(123, "a", "varname")
    correct = '''{
        char filename[300];
        sprintf(filename, "a%.4d.hdf5", varname);
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        myrank = 0;
        nranks = 1;
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out

    set_mpi(True)
    out = a.call_output_cabana_pir(123, "a")
    correct = '''{
        char filename[300] = "a";
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &nranks );
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out

    set_mpi(False)
    out = a.call_output_cabana_pir(123, "a")
    correct = '''{
        char filename[300] = "a";
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        myrank = 0;
        nranks = 1;
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out

    set_mpi(True)
    out = a.call_output_cabana_pir(123, "\"a\"")
    correct = '''{
        char filename[300] = "a";
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &nranks );
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out

    set_mpi(False)
    out = a.call_output_cabana_pir(123, "\"a\"")
    correct = '''{
        char filename[300] = "a";
        Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
        int myrank, nranks;
        myrank = 0;
        nranks = 1;
        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);
    }'''
    assert correct in out
    reset_for_tests()

def test_phdf5_get_linked_libraries():
    a = PHDF5_IO.PHDF5_IO()
    assert a.get_linked_libraries() == ["${HDF5_C_LIBRARIES}"]
