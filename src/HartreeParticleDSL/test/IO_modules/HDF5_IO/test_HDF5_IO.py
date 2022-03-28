import HartreeParticleDSL.IO_modules.HDF5_IO.hdf5_IO as hdf5_IO
import pytest
from HartreeParticleDSL.HartreeParticleDSL import Particle

def test_hdf5_init():
    a = hdf5_IO.HDF5_IO(indent=8)
    assert a._inputs == {}
    assert a._outputs == {}
    assert a._indent == 8
    assert a._current_indent == 0
    a = hdf5_IO.HDF5_IO()
    assert a._indent == 4
    assert a._current_indent == 0

def test_hdf5_add_input():
    a = hdf5_IO.HDF5_IO()
    z = "h_in"
    y = "part_in"
    a.add_input(z, y)
    assert a._inputs[z] == y

def test_hdf5_add_output():
    a = hdf5_IO.HDF5_IO()
    z = "h_out"
    y = "part_out"
    a.add_output(z, y)
    assert a._outputs[z] == y

def test_hdf5_indent():
    a = hdf5_IO.HDF5_IO(indent=1)
    assert a.indent() == ""

def test_hdf5_increment_indent():
    a = hdf5_IO.HDF5_IO(indent=1)
    a.increment_indent()
    assert a.indent() == " "

def test_hdf5_decrement_indent():
    a = hdf5_IO.HDF5_IO(indent=1)
    a.increment_indent()
    assert a.indent() == " "
    a.decrement_indent()
    assert a.indent() == ""

def test_get_includes_c():
    a = hdf5_IO.HDF5_IO()
    x = a.get_includes_c()
    assert "<stdlib.h>" in x
    assert "\"hdf5.h\"" in x

def test_gen_code_c():
    a = hdf5_IO.HDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_c(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    correct = '''struct part* hdf5_input(const char* filename, struct config_type* config){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Failed to open %s\\n", filename);
        exit(1);
    }
    hid_t thing_test_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    if( thing_test_var < 0){
        printf("Failed to find dataset thing\\n");
        exit(1);
    }
    H5Dclose(thing_test_var);
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
    hid_t part_count_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    hid_t space = H5Dget_space(part_count_var);
    int ndims = H5Sget_simple_extent_ndims(space);
    if( ndims != 1 ){
        printf("Don't yet support multidimensional datasets.\\n");
        exit(1);
    }
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space, dims, NULL);
    int num_parts = dims[0];
    struct part *parts = (struct part*) malloc(sizeof(struct part) * num_parts);
    config->space.nparts = num_parts;

    hsize_t shape[2], offsets[2];
    int rank = 2;
    shape[0] = num_parts;
    shape[1] = 1;
    offsets[0] = 0;
    offsets[1] = 0;
    hid_t memspace = H5Screate_simple(rank, shape, NULL);
    hid_t filespace;
    hid_t thing_read_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    filespace = H5Dget_space(thing_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    int* thing_temp_array = (int*) malloc(sizeof(int) * num_parts);
    H5Dread(thing_read_var, H5T_STD_I32LE, memspace, filespace, H5P_DEFAULT, thing_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].thing = thing_temp_array[i];
    }
    free(thing_temp_array);
    H5Dclose(thing_read_var);
    hid_t pos_x_read_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    filespace = H5Dget_space(pos_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* pos_x_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.position[0] = pos_x_temp_array[i];
    }
    free(pos_x_temp_array);
    H5Dclose(pos_x_read_var);
    hid_t pos_y_read_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    filespace = H5Dget_space(pos_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* pos_y_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.position[1] = pos_y_temp_array[i];
    }
    free(pos_y_temp_array);
    H5Dclose(pos_y_read_var);
    hid_t pos_z_read_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    filespace = H5Dget_space(pos_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* pos_z_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.position[2] = pos_z_temp_array[i];
    }
    free(pos_z_temp_array);
    H5Dclose(pos_z_read_var);
    hid_t vel_x_read_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    filespace = H5Dget_space(vel_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* vel_x_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.velocity[0] = vel_x_temp_array[i];
    }
    free(vel_x_temp_array);
    H5Dclose(vel_x_read_var);
    hid_t vel_y_read_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    filespace = H5Dget_space(vel_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* vel_y_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.velocity[1] = vel_y_temp_array[i];
    }
    free(vel_y_temp_array);
    H5Dclose(vel_y_read_var);
    hid_t vel_z_read_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    filespace = H5Dget_space(vel_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    double* vel_z_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        parts[i].core_part.velocity[2] = vel_z_temp_array[i];
    }
    free(vel_z_temp_array);
    H5Dclose(vel_z_read_var);

    H5Fclose(file_id);
    return parts;
}

'''
    x = a.gen_code_c(part)
    assert x == correct
    a = hdf5_IO.HDF5_IO()
    a.add_output("thing", "thing")
    a.add_output("pos_x", "core_part.position[0]")
    a.add_output("pos_y", "core_part.position.y")
    a.add_output("pos_z", "core_part.position.z")
    a.add_output("vel_x", "core_part.velocity[0]")
    a.add_output("vel_y", "core_part.velocity[1]")
    a.add_output("vel_z", "core_part.velocity[2]")
    x = a.gen_code_c(part)
    correct = '''
void hdf5_output(struct part* parts, struct config_type* config, const char* filename){
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Couldn't create HDF5 file\\n");
        exit(1);
    }

    hsize_t dims[1];
    dims[0] = config->space.nparts;
    hid_t dim = H5Screate_simple(1, dims, NULL);
    hid_t thing_output_field = H5Dcreate2(file_id, "thing", H5T_STD_I32LE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int* thing_output_array = (int *) malloc(sizeof(int) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        thing_output_array[i] = parts[i].thing;
    }

    H5Dwrite(thing_output_field, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thing_output_array);
    free(thing_output_array);
    H5Dclose(thing_output_field);
    hid_t pos_x_output_field = H5Dcreate2(file_id, "pos_x", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_x_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        pos_x_output_array[i] = parts[i].core_part.position[0];
    }

    H5Dwrite(pos_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_x_output_array);
    free(pos_x_output_array);
    H5Dclose(pos_x_output_field);
    hid_t pos_y_output_field = H5Dcreate2(file_id, "pos_y", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_y_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        pos_y_output_array[i] = parts[i].core_part.position[1];
    }

    H5Dwrite(pos_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_y_output_array);
    free(pos_y_output_array);
    H5Dclose(pos_y_output_field);
    hid_t pos_z_output_field = H5Dcreate2(file_id, "pos_z", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_z_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        pos_z_output_array[i] = parts[i].core_part.position[2];
    }

    H5Dwrite(pos_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_z_output_array);
    free(pos_z_output_array);
    H5Dclose(pos_z_output_field);
    hid_t vel_x_output_field = H5Dcreate2(file_id, "vel_x", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_x_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        vel_x_output_array[i] = parts[i].core_part.velocity[0];
    }

    H5Dwrite(vel_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_x_output_array);
    free(vel_x_output_array);
    H5Dclose(vel_x_output_field);
    hid_t vel_y_output_field = H5Dcreate2(file_id, "vel_y", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_y_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        vel_y_output_array[i] = parts[i].core_part.velocity[1];
    }

    H5Dwrite(vel_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_y_output_array);
    free(vel_y_output_array);
    H5Dclose(vel_y_output_field);
    hid_t vel_z_output_field = H5Dcreate2(file_id, "vel_z", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_z_output_array = (double *) malloc(sizeof(double) * config->space.nparts);
    for( int i = 0; i < config->space.nparts; i++){
        vel_z_output_array[i] = parts[i].core_part.velocity[2];
    }

    H5Dwrite(vel_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_z_output_array);
    free(vel_z_output_array);
    H5Dclose(vel_z_output_field);

    H5Fclose(file_id);
}

'''
    assert x == correct

def test_call_input_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_input_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_input(\"myfile.hdf5\", config);"
    assert x == correct

def test_call_output_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_output_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_output(parts, config, \"myfile.hdf5\");"
    assert x == correct

def test_gen_code_fdps():
    a = hdf5_IO.HDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_fdps(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    x = a.gen_code_fdps(part)
    correct = '''    void hdf5_input(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char *filename){
        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if( file_id < 0 ){
            printf("Failed to open %s\\n", filename);
            exit(1);
        }
        hid_t thing_test_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
        if( thing_test_var < 0){
            printf("Failed to find dataset thing\\n");
            exit(1);
        }
        H5Dclose(thing_test_var);
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
        hid_t part_count_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
        hid_t space = H5Dget_space(part_count_var);
        int ndims = H5Sget_simple_extent_ndims(space);
        if( ndims != 1 ){
            printf("Don't yet support multidimensional datasets.\\n");
            exit(1);
        }
        hsize_t[1] dims;

        H5Fclose(file_id);
    }

'''
    assert x == correct

def test_call_input_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_input_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_input(\"myfile.hdf5\", config);"
    assert x == correct

def test_call_output_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_output_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_output(parts, config, \"myfile.hdf5\");"
    assert x == correct

def test_gen_code_fdps():
    a = hdf5_IO.HDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_fdps(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    x = a.gen_code_fdps(part)
    correct = '''    void hdf5_input(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char *filename){
        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        if( file_id < 0 ){
            printf("Failed to open %s\\n", filename);
            exit(1);
        }
        hid_t thing_test_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
        if( thing_test_var < 0){
            printf("Failed to find dataset thing\\n");
            exit(1);
        }
        H5Dclose(thing_test_var);
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
        hid_t part_count_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
        hid_t space = H5Dget_space(part_count_var);
        int ndims = H5Sget_simple_extent_ndims(space);
        if( ndims != 1 ){
            printf("Don't yet support multidimensional datasets.\\n");
            exit(1);
        }
        hsize_t[1] dims;

        H5Fclose(file_id);
    }

'''
    assert x == correct

def test_call_input_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_input_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_input(\"myfile.hdf5\", config);"
    assert x == correct

def test_call_output_c():
    a = hdf5_IO.HDF5_IO()
    x = a.call_output_c(1000, "\"myfile.hdf5\"")
    correct = "hdf5_output(parts, config, \"myfile.hdf5\");"
    assert x == correct

def test_gen_code_cabana():
    a = hdf5_IO.HDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_cabana(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    x = a.gen_code_cabana(part)
    correct = '''template <class aosoa_class, class aosoa_host_class> void hdf5_input(aosoa_class particle_system, aosoa_host_class parts_host, config_type& config, const char *filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Failed to open %s\\n", filename);
        exit(1);
    }
    hid_t thing_test_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    if( thing_test_var < 0){
        printf("Failed to find dataset thing\\n");
        exit(1);
    }
    H5Dclose(thing_test_var);
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
    hid_t part_count_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    hid_t space = H5Dget_space(part_count_var);
    int ndims = H5Sget_simple_extent_ndims(space);
    if( ndims != 1 ){
        printf("Don't yet support multidimensional datasets.\\n");
        exit(1);
    }
    hsize_t[1] dims;
    H5Sget_simple_extent_dims(space, dims, NULL);
    int num_parts = dims[0];
    int new_size = static_cast<int>(num_parts);
    particle_system.resize(new_size);
    parts_host.resize(new_size);
    hsize_t shape[2], offsets[2];
    int rank = 2;
    shape[0] = num_parts;
    shape[1] = 1;
    offsets[0] = 0;
    offsets[1] = 0;
    hid_t memspace = H5Screate_simple(rank, shape, NULL);
    hid_t filespace;
    H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    hid_t thing_read_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    filespace = H5Dget_space(thing_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto thing_slice = Cabana::slice<thing>(parts_host);
    int* thing_temp_array = (int*) malloc(sizeof(int) * num_parts);
    H5Dread(thing_read_var, H5T_STD_I32LE, memspace, filespace, H5P_DEFAULT, thing_temp_array);
    for( int i = 0; i < num_parts; i++){
        thing_slice(i) = thing_temp_array[i];
    }
    free(thing_temp_array);
    H5Dclose(thing_read_var);
    hid_t pos_x_read_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    filespace = H5Dget_space(pos_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto pos_x_slice = Cabana::slice<core_part_position>(parts_host);
    double* pos_x_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        pos_x_slice(i, 0) = pos_x_temp_array[i];
    }
    free(pos_x_temp_array);
    H5Dclose(pos_x_read_var);
    hid_t pos_y_read_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    filespace = H5Dget_space(pos_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto pos_y_slice = Cabana::slice<core_part_position>(parts_host);
    double* pos_y_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        pos_y_slice(i, 1) = pos_y_temp_array[i];
    }
    free(pos_y_temp_array);
    H5Dclose(pos_y_read_var);
    hid_t pos_z_read_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    filespace = H5Dget_space(pos_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto pos_z_slice = Cabana::slice<core_part_position>(parts_host);
    double* pos_z_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(pos_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        pos_z_slice(i, 2) = pos_z_temp_array[i];
    }
    free(pos_z_temp_array);
    H5Dclose(pos_z_read_var);
    hid_t vel_x_read_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    filespace = H5Dget_space(vel_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto vel_x_slice = Cabana::slice<core_part_velocity>(parts_host);
    double* vel_x_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        vel_x_slice(i, 0) = vel_x_temp_array[i];
    }
    free(vel_x_temp_array);
    H5Dclose(vel_x_read_var);
    hid_t vel_y_read_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    filespace = H5Dget_space(vel_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto vel_y_slice = Cabana::slice<core_part_velocity>(parts_host);
    double* vel_y_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        vel_y_slice(i, 1) = vel_y_temp_array[i];
    }
    free(vel_y_temp_array);
    H5Dclose(vel_y_read_var);
    hid_t vel_z_read_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    filespace = H5Dget_space(vel_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    auto vel_z_slice = Cabana::slice<core_part_velocity>(parts_host);
    double* vel_z_temp_array = (double*) malloc(sizeof(double) * num_parts);
    H5Dread(vel_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        vel_z_slice(i, 2) = vel_z_temp_array[i];
    }
    free(vel_z_temp_array);
    H5Dclose(vel_z_read_var);

    H5Fclose(file_id);
    Cabana::deep_copy(particle_system, parts_host);
}

'''
    assert correct == x
    a = hdf5_IO.HDF5_IO()
    a.add_output("thing", "thing")
    a.add_output("pos_x", "core_part.position[0]")
    a.add_output("pos_y", "core_part.position.y")
    a.add_output("pos_z", "core_part.position.z")
    a.add_output("vel_x", "core_part.velocity[0]")
    a.add_output("vel_y", "core_part.velocity[1]")
    a.add_output("vel_z", "core_part.velocity[2]")
    x = a.gen_code_cabana(part)

    correct = '''
template <class aosoa_class> void hdf5_output(aosoa_class particle_aosoa, config_type& config, const char* filename){
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if( file_id < 0 ){
        printf(\"Couldn't create HDF5 file\\n\");
        exit(1);
    }

    hsize_t dims[1];
    dims[0] = particle_aosoa.size();
    hid_t dim = H5Screate_simple(1, dims, NULL);
    hid_t thing_output_field = H5Dcreate2(file_id, \"thing\", H5T_STD_I32LE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int* thing_output_array = (int *) malloc(sizeof(int) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        thing_output_array[i] = Cabana::get<thing>(part);
    }

    H5Dwrite(thing_output_field, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thing_output_array);
    free(thing_output_array);
    H5Dclose(thing_output_field);
    hid_t pos_x_output_field = H5Dcreate2(file_id, \"pos_x\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_x_output_array[i] = Cabana::get<core_part_position>(part, 0);
    }

    H5Dwrite(pos_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_x_output_array);
    free(pos_x_output_array);
    H5Dclose(pos_x_output_field);
    hid_t pos_y_output_field = H5Dcreate2(file_id, \"pos_y\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_y_output_array[i] = Cabana::get<core_part_position>(part, 1);
    }

    H5Dwrite(pos_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_y_output_array);
    free(pos_y_output_array);
    H5Dclose(pos_y_output_field);
    hid_t pos_z_output_field = H5Dcreate2(file_id, \"pos_z\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* pos_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        pos_z_output_array[i] = Cabana::get<core_part_position>(part, 2);
    }

    H5Dwrite(pos_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_z_output_array);
    free(pos_z_output_array);
    H5Dclose(pos_z_output_field);
    hid_t vel_x_output_field = H5Dcreate2(file_id, \"vel_x\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_x_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_x_output_array[i] = Cabana::get<core_part_velocity>(part, 0);
    }

    H5Dwrite(vel_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_x_output_array);
    free(vel_x_output_array);
    H5Dclose(vel_x_output_field);
    hid_t vel_y_output_field = H5Dcreate2(file_id, \"vel_y\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_y_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_y_output_array[i] = Cabana::get<core_part_velocity>(part, 1);
    }

    H5Dwrite(vel_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_y_output_array);
    free(vel_y_output_array);
    H5Dclose(vel_y_output_field);
    hid_t vel_z_output_field = H5Dcreate2(file_id, \"vel_z\", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double* vel_z_output_array = (double *) malloc(sizeof(double) * particle_aosoa.size());
    for( int i = 0; i < particle_aosoa.size(); i++){
        auto part = particle_aosoa.getTuple(i);
        vel_z_output_array[i] = Cabana::get<core_part_velocity>(part, 2);
    }

    H5Dwrite(vel_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_z_output_array);
    free(vel_z_output_array);
    H5Dclose(vel_z_output_field);

    H5Fclose(file_id);
}

'''
    assert correct == x

def test_gen_code_fdps():
    a = hdf5_IO.HDF5_IO()
    correct = ""
    part = Particle()
    x = a.gen_code_fdps(part)
    assert x == correct
    part.add_element("thing", "c_int")
    a.add_input("thing", "thing")
    a.add_input("pos_x", "core_part.position[0]")
    a.add_input("pos_y", "core_part.position.y")
    a.add_input("pos_z", "core_part.position.z")
    a.add_input("vel_x", "core_part.velocity[0]")
    a.add_input("vel_y", "core_part.velocity[1]")
    a.add_input("vel_z", "core_part.velocity[2]")
    x = a.gen_code_fdps(part)
    correct = '''void hdf5_input(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char *filename){
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Failed to open %s\\n", filename);
        exit(1);
    }
    hid_t thing_test_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    if( thing_test_var < 0){
        printf("Failed to find dataset thing\\n");
        exit(1);
    }
    H5Dclose(thing_test_var);
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
    hid_t part_count_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    hid_t space = H5Dget_space(part_count_var);
    int ndims = H5Sget_simple_extent_ndims(space);
    if( ndims != 1 ){
        printf("Don't yet support multidimensional datasets.\\n");
        exit(1);
    }
    hsize_t[1] dims;
    H5Sget_simple_extent_dims(space, dims, NULL);
    int num_parts = dims[0];
    particle_system.initialize();
    particle_system.setNumberOfParticleLocal(num_parts);
    config.space.nparts = num_parts;
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.setBoundaryCondition ( PS::BOUNDARY_CONDITION_PERIODIC_XYZ );
    dinfo.setPosRootDomain(PS::F64vec(config.space.box_dims.x_min,
                                         config.space.box_dims.y_min,
                                         config.space.box_dims.z_min),
                              PS::F64vec(config.space.box_dims.x_max,
                                         config.space.box_dims.y_max,
                                         config.space.box_dims.z_max));

    hsize_t shape[2], offsets[2];
    int rank = 2;
    shape[0] = num_parts;
    shape[1] = 1;
    offsets[0] = 0;
    offsets[1] = 0;
    hid_t memspace = H5Screate_simple(rank, shape, NULL);
    hid_t filespace;
    H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    hid_t thing_read_var = H5Dopen2(file_id, "thing", H5P_DEFAULT);
    filespace = H5Dget_space(thing_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::S32* thing_temp_array = (PS::S32*) malloc(sizeof(PS::S32) * num_parts);
    H5Dread(thing_read_var, H5T_STD_I32LE, memspace, filespace, H5P_DEFAULT, thing_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].thing = thing_temp_array[i];
    }
    free(thing_temp_array);
    H5Dclose(thing_read_var);
    hid_t pos_x_read_var = H5Dopen2(file_id, "pos_x", H5P_DEFAULT);
    filespace = H5Dget_space(pos_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* pos_x_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(pos_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.position.x = pos_x_temp_array[i];
    }
    free(pos_x_temp_array);
    H5Dclose(pos_x_read_var);
    hid_t pos_y_read_var = H5Dopen2(file_id, "pos_y", H5P_DEFAULT);
    filespace = H5Dget_space(pos_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* pos_y_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(pos_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.position.y = pos_y_temp_array[i];
    }
    free(pos_y_temp_array);
    H5Dclose(pos_y_read_var);
    hid_t pos_z_read_var = H5Dopen2(file_id, "pos_z", H5P_DEFAULT);
    filespace = H5Dget_space(pos_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* pos_z_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(pos_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, pos_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.position.z = pos_z_temp_array[i];
    }
    free(pos_z_temp_array);
    H5Dclose(pos_z_read_var);
    hid_t vel_x_read_var = H5Dopen2(file_id, "vel_x", H5P_DEFAULT);
    filespace = H5Dget_space(vel_x_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* vel_x_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(vel_x_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_x_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.velocity[0] = vel_x_temp_array[i];
    }
    free(vel_x_temp_array);
    H5Dclose(vel_x_read_var);
    hid_t vel_y_read_var = H5Dopen2(file_id, "vel_y", H5P_DEFAULT);
    filespace = H5Dget_space(vel_y_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* vel_y_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(vel_y_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_y_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.velocity[1] = vel_y_temp_array[i];
    }
    free(vel_y_temp_array);
    H5Dclose(vel_y_read_var);
    hid_t vel_z_read_var = H5Dopen2(file_id, "vel_z", H5P_DEFAULT);
    filespace = H5Dget_space(vel_z_read_var);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);
    PS::F64* vel_z_temp_array = (PS::F64*) malloc(sizeof(PS::F64) * num_parts);
    H5Dread(vel_z_read_var, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, vel_z_temp_array);
    for( int i = 0; i < num_parts; i++){
        particle_system[i].core_part.velocity[2] = vel_z_temp_array[i];
    }
    free(vel_z_temp_array);
    H5Dclose(vel_z_read_var);

    H5Fclose(file_id);
}

'''
    assert correct == x
    a = hdf5_IO.HDF5_IO()
    a.add_output("thing", "thing")
    a.add_output("pos_x", "core_part.position[0]")
    a.add_output("pos_y", "core_part.position.y")
    a.add_output("pos_z", "core_part.position.z")
    a.add_output("vel_x", "core_part.velocity[0]")
    a.add_output("vel_y", "core_part.velocity[1]")
    a.add_output("vel_z", "core_part.velocity[2]")
    x = a.gen_code_fdps(part)
    correct = '''
void hdf5_output(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char* filename){
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if( file_id < 0 ){
        printf("Couldn't create HDF5 file\\n");
        exit(1);
    }

    hsize_t dims[1];
    dims[0] = config.space.nparts;
    hid_t dim = H5Screate_simple(1, dims, NULL);
    hid_t thing_output_field = H5Dcreate2(file_id, "thing", H5T_STD_I32LE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::S32* thing_output_array = (PS::S32 *) malloc(sizeof(PS::S32) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        thing_output_array[i] = particle_system[i].thing;
    }

    H5Dwrite(thing_output_field, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thing_output_array);
    free(thing_output_array);
    H5Dclose(thing_output_field);
    hid_t pos_x_output_field = H5Dcreate2(file_id, "pos_x", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* pos_x_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        pos_x_output_array[i] = particle_system[i].core_part.position.x;
    }

    H5Dwrite(pos_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_x_output_array);
    free(pos_x_output_array);
    H5Dclose(pos_x_output_field);
    hid_t pos_y_output_field = H5Dcreate2(file_id, "pos_y", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* pos_y_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        pos_y_output_array[i] = particle_system[i].core_part.position.y;
    }

    H5Dwrite(pos_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_y_output_array);
    free(pos_y_output_array);
    H5Dclose(pos_y_output_field);
    hid_t pos_z_output_field = H5Dcreate2(file_id, "pos_z", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* pos_z_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        pos_z_output_array[i] = particle_system[i].core_part.position.z;
    }

    H5Dwrite(pos_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos_z_output_array);
    free(pos_z_output_array);
    H5Dclose(pos_z_output_field);
    hid_t vel_x_output_field = H5Dcreate2(file_id, "vel_x", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* vel_x_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        vel_x_output_array[i] = particle_system[i].core_part.velocity[0];
    }

    H5Dwrite(vel_x_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_x_output_array);
    free(vel_x_output_array);
    H5Dclose(vel_x_output_field);
    hid_t vel_y_output_field = H5Dcreate2(file_id, "vel_y", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* vel_y_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        vel_y_output_array[i] = particle_system[i].core_part.velocity[1];
    }

    H5Dwrite(vel_y_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_y_output_array);
    free(vel_y_output_array);
    H5Dclose(vel_y_output_field);
    hid_t vel_z_output_field = H5Dcreate2(file_id, "vel_z", H5T_NATIVE_DOUBLE, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    PS::F64* vel_z_output_array = (PS::F64 *) malloc(sizeof(PS::F64) * config.space.nparts);
    for( int i = 0; i < config.space.nparts; i++){
        vel_z_output_array[i] = particle_system[i].core_part.velocity[2];
    }

    H5Dwrite(vel_z_output_field, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel_z_output_array);
    free(vel_z_output_array);
    H5Dclose(vel_z_output_field);

    H5Fclose(file_id);
}

'''
    assert correct == x
