from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
        c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
from HartreeParticleDSL.backends.Cabana_PIR_backend.Cabana_PIR_IO_Mixin import Cabana_PIR_IO_Mixin
from HartreeParticleDSL.HartreeParticleDSLExceptions import UnsupportedCodeError
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL


# New style module to go with Cabana_PIR and MPI.
class PHDF5_IO(IO_Module, Cabana_PIR_IO_Mixin):
    '''Implementation of the Parallel HDF5 IO Module.

    This module currently only supports the ``Cabana_PIR`` backend.
    
    This module uses HDF5 over MPI to perform File I/O across multiple nodes.

    To use this module, it needs to know how and which elements of particles
    to output. It currently doesn't support outputting elements of the config.

    >>> import HartreeParticleDSL.IO_modules.PHDF5_IO.PHDF5_IO as PHDF5_IO
    >>> import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
    >>> io_module = PHDF5_IO.PHDF5_IO()
    >>> io_module.add_input("x_positions", "core_part.position.x")
    >>> io_module.add_output("x_positions", "core_part.position.x")
    >>> HartreeParticleDSL.set_io_modules(io_module, io_module)

    This would read the data in `x_positions` into the core particle's 
    x position elements, and write the same back at the appropriate time.

    At the moment, the PHDF5 module only supports the following types:
    `int`, `double`, `float, `int64_t`, `int32_t`, `int8_t`, and `bool`.

    If you need other types, please raise an issue on github and they can
    be easily added provided HDF5 supports them.

    '''
    type_map = {"int": "H5T_STD_I32LE",
                "double" : "H5T_NATIVE_DOUBLE",
                "float" : "H5T_NATIVE_FLOAT",
                "int64_t" : "H5T_STD_I64LE",
                "int32_t" : "H5T_STD_I32LE",
                "int8_t" : "H5T_STD_I8LE",
                "bool" : "H5T_NATIVE_HBOOL"}

    def __init__(self, indent=4):
        super().__init__()
        self._inputs = {}
        self._outputs = {}
        self._indent = indent
        self._current_indent = 0

    def get_linked_libraries(self):
        return ["${HDF5_C_LIBRARIES}"]

    def add_input(self, hdf_input: str, particle_input: str) -> None:
        '''
        Tells the PHDF5 IO module to read in a particular element
        from the HDF5 input file.

        :param str hdf_input: The name of the field in the HDF5 input file(s).
        :param str particle_input: The element of the particle data structure \
                                   to place the input data into.
        '''
        self._inputs[hdf_input] = particle_input

    def add_output(self, hdf_output: str, particle_output: str) -> None:
        '''
        Tells the PHDF5 IO m odule to output an element of the particle
        data into a field in a HDF5 output file.

        :param str hdf_output: The name of the field in the HDF5 output files.
        :param str particle_output: The element of the particle data structure \
                                    containing the data to be output.
        '''
        self._outputs[hdf_output] = particle_output

    def indent(self) -> str:
        '''
        :returns: the string representing the current indentation.
        :rtype: str
        '''
        return " " * self._current_indent

    def _increment_indent(self):
        '''
        Increments the current indentation level.
        '''
        self._current_indent = self._current_indent + self._indent

    def _decrement_indent(self):
        '''
        Decrements the current indentation level.
        '''
        self._current_indent = self._current_indent - self._indent

    def get_header_includes_cabana_pir(self):
        '''
        Returns the includes required for the header to use this IO module for
        Cabana PIR.

        :returns: The header includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"hdf5.h\"")
        if HartreeParticleDSL.get_mpi():
            includes.append("\"mpi.h\"")
        return includes

    def get_includes_cabana_pir(self):
        '''
        Returns the includes required to use this IO module for Cabana PIR.

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"hdf5.h\"")
        if HartreeParticleDSL.get_mpi():
            includes.append("\"mpi.h\"")
        return includes

    def gen_code_cabana_pir(self, part_type):
        '''
        Returns the Cabana_PIR code required for this IO module.

        :param part_type: The particle type used in the system.
        :type part_type: Particle
        '''
        code = ""
        if len(self._inputs) > 0:
            # Generate the function to get the box size from the hdf5 file first.
            # ============================ get_box_size ===========================
            code = code + self.indent() + "void get_box_size(boundary &box, const char* filename){\n"
            self._increment_indent()
            code = code + self.indent() + "hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self._increment_indent()
            code = code + self.indent() + "std::cout << \"Failed to open file \" << filename << \"\\n\";\n"
            code = code + self.indent() + "exit(1);\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n"
            code = code + self.indent() + "hid_t temp_space;\n"
            code = code + self.indent() + "hsize_t dims[1];\n\n"

            # Load the box size
            code = code + self.indent() + "hid_t boxsize = H5Dopen2(file_id, \"Box_Size\", H5P_DEFAULT);\n"
            code = code + self.indent() + "temp_space = H5Dget_space(boxsize);\n"
            code = code + self.indent() + "H5Sget_simple_extent_dims(temp_space, dims, NULL);\n"
            code = code + self.indent() + "int size_box = dims[0];\n"
            code = code + self.indent() + "double* box_temp = (double*) malloc(sizeof(double) * size_box);\n"
            code = code + self.indent() + "H5Dread(boxsize, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_temp);\n"
            code = code + self.indent() + "box.x_min = box_temp[0];\n"
            code = code + self.indent() + "box.x_max = box_temp[1];\n"
            code = code + self.indent() + "if(size_box > 2){\n"
            self._increment_indent()
            code = code + self.indent() + "box.y_min = box_temp[2];\n"
            code = code + self.indent() + "box.y_max = box_temp[3];\n"
            self._decrement_indent()
            code = code + self.indent() + "}else{\n"
            self._increment_indent()
            code = code + self.indent() + "box.y_min = 0.0;\n"
            code = code + self.indent() + "box.y_max = 0.0;\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n"
            code = code + self.indent() + "if(size_box > 4){\n"
            self._increment_indent()
            code = code + self.indent() + "box.z_min = box_temp[4];\n"
            code = code + self.indent() + "box.z_max = box_temp[5];\n"
            self._decrement_indent()
            code = code + self.indent() + "}else{\n"
            self._increment_indent()
            code = code + self.indent() + "box.z_min = 0.0;\n"
            code = code + self.indent() + "box.z_max = 0.0;\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n\n"

            # Clean up
            code = code + self.indent() + "free(box_temp);\n"
            code = code + self.indent() + "H5Dclose(boxsize);\n"
            code = code + self.indent() + "H5Fclose(file_id);\n"
            self._decrement_indent()
            code = code + "}\n\n"
            # ============================ get_box_size ===========================


            # Create the hdf5_input function.
            # ============================= HDF5_INPUT ============================
            code = code + self.indent() + "template <class aosoa_class, class aosoa_host_class> void hdf5_input(aosoa_class &non_host_aosoa,\n"
            code = code + self.indent() + "    aosoa_host_class &particle_aosoa, config_type &config, boundary &box, const char* filename){\n" #TODO
            self._increment_indent()
            # Open the file
            code = code + self.indent() + "hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self._increment_indent()
            code = code + self.indent() + "std::cout << \"Failed to open file \" << filename << \"\\n\";\n"
            code = code + self.indent() + "exit(1);\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n"
            code = code + self.indent() + "hid_t temp_space;\n"
            code = code + self.indent() + "hsize_t dims[1];\n\n"

            # Check all the fields we asked for are in the file.
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_test_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "if( {0}_test_var < 0)".format(key) + "{\n"
                self._increment_indent()
                code = code + self.indent() + "printf(\"Failed to find dataset {0}\\n\");\n".format(key)
                code = code + self.indent() + "exit(1);\n"
                self._decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "H5Dclose({0}_test_var);\n".format(key)

            # All of the fields are in the HDF5 file as requested, fetch the particle count.
            code = code + self.indent() + "hid_t part_count_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(list(self._inputs.keys())[0])
            code = code + self.indent() + "hid_t space = H5Dget_space(part_count_var);\n"
            code = code + self.indent() + "int ndims = H5Sget_simple_extent_ndims(space);\n"

            # We need to have the positions as one of the things to be loaded in the particles to be able to do this.
            positions = [None, None, None]
            for key, value in self._inputs.items():
                if "core_part.position.x" in value or "core_part.position[0]" in value:
                    positions[0] = key
                if "core_part.position.y" in value or "core_part.position[1]" in value:
                    positions[1] = key
                if "core_part.position.z" in value or "core_part.position[2]" in value:
                    positions[2] = key

            # Start by loading the positions
            conditions = []
            for index, key in enumerate(positions):
                if key is not None:
                    code = code + self.indent() + f"hid_t {key}_var = H5Dopen2" + "(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                    code = code + self.indent() + f"temp_space = H5Dget_space({key}_var);\n"
                    code = code + self.indent() + "H5Sget_simple_extent_dims(temp_space, dims, NULL);\n"
                    code = code + self.indent() + f"double* {key}_temp_array = (double*) malloc(sizeof(double) * dims[0]);\n"
                    code = code + self.indent() + f"H5Dread({key}_var, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, {key}_temp_array);\n"
                    code = code + self.indent() + f"H5Dclose({key}_var);\n"

                    # Condition is < X_max and >= X_min
                    dimension = "x"
                    if index == 1:
                        dimension = "y"
                    elif index == 2:
                        dimension = "z"
                    conditions.append(f"({key}_temp_array[i] < box.local_{dimension}_max && {key}_temp_array[i] >= box.local_{dimension}_min)")

            # Now we loaded the positions and set up the positions, lets copy the positions into the particles.
            
            condition = " && ".join(conditions)
            # First count how many we have and resize the AoSoA
            code = code + self.indent() + "int global_parts = dims[0];\n"
            code = code + self.indent() + "int num_parts = 0;\n"
            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "for(int i = 0; i < global_parts; i++){\n"
                self._increment_indent()
                code = code + self.indent() + "if(" + condition + "){\n"
                self._increment_indent()
                code = code + self.indent() + "num_parts++;\n"
                self._decrement_indent()
                code = code + self.indent() + "}\n"
                self._decrement_indent()
                code = code + self.indent() + "}\n"
            else:
                code = code + self.indent() + "num_parts = global_parts;\n"
            code = code + self.indent() + "config.config_host(0).space.nparts = global_parts;\n"
            code = code + self.indent() + "int new_size = static_cast<int>(num_parts);\n"
            code = code + self.indent() + "particle_aosoa.resize(new_size);\n"
            code = code + self.indent() + "non_host_aosoa.resize(new_size);\n"

            # Now copy in the positions
            code = code + self.indent() + "int counter = 0;\n"
            code = code + self.indent() + "auto pos_slice = Cabana::slice<core_part_position>(particle_aosoa);\n"
            code = code + self.indent() + "for(int i = 0; i < global_parts; i++){\n"
            self._increment_indent()
            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "if(" + condition + "){\n"
                self._increment_indent()
            for index, key in enumerate(positions):
                if key is not None:
                    code = code + self.indent() + f"pos_slice(counter, {index}) = {key}_temp_array[i];\n"
            code = code + self.indent() + "counter++;\n"
            if HartreeParticleDSL.get_mpi():
                self._decrement_indent()
                code = code + self.indent() + "}\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n"

            # Ok, now we did the positions we should do anything else that is needed to be read in
            for key in self._inputs.keys():
                # Skip the positions now.
                if key in positions:
                    continue
                code = code + self.indent() + "hid_t {0}_read_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "space = H5Dget_space({0}_read_var);\n".format(key)
#                code = code + self.indent() + "H5Sselect_hyperslab(space, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
                # Find the type of the variable from the particle structure
                part_elem = self._inputs[key]
                part_indexing = ""
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 0"
                        elem_type = "double"
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 1"
                        elem_type = "double"
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 2"
                        elem_type = "double"
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = PHDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    raise UnsupportedCodeError(f"{part_elem} element not supported in Parallel HDF5 IO.")
                elem_type = Cabana_PIR._type_map[elem_type]
                code = code + self.indent() + "auto {0}_slice = Cabana::slice<{1}>(particle_aosoa);\n".format(key, part_elem)
                code = code + self.indent() + "{0}* {1}_temp_array = ({0}*) malloc(sizeof({0}) * global_parts);\n".format(elem_type, key)
                code = code + self.indent() + "H5Dread({0}_read_var, {1}, H5S_ALL, space, H5P_DEFAULT, {0}_temp_array);\n".format(key, h5_type)
                code = code + self.indent() + "counter = 0;\n"
                code = code + self.indent() + "for(int i = 0; i < global_parts; i++){\n"
                self._increment_indent()
                if HartreeParticleDSL.get_mpi():
                    code = code + self.indent() + "if(" + condition + "){\n"
                    self._increment_indent()
                code = code + self.indent() + "{0}_slice(counter{2}) = {1}_temp_array[i];\n".format(key, key, part_indexing)
                code = code + self.indent() + "counter++;\n"
                if HartreeParticleDSL.get_mpi():
                    self._decrement_indent()
                    code = code + self.indent() + "}\n"
                self._decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "free({0}_temp_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_read_var);\n".format(key)

            # At the end we need to cleanup.
            for index, key in enumerate(positions):
                if key is not None:
                    code = code + self.indent() + f"free({key}_temp_array);\n"
            code = code + self.indent() + "H5Fclose(file_id);\n"
            code = code + self.indent() + "Cabana::deep_copy(non_host_aosoa, particle_aosoa);\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n\n"
            # ============================= HDF5_INPUT ============================


        if len(self._outputs) > 0:
            # Create the hdf5_output function.
            # ============================= HDF5_OUTPUT ============================
            code = code + self.indent() + "template <class aosoa_host_class> void hdf5_output(aosoa_host_class &particle_aosoa,\n"
            code = code + self.indent() + "    const char* filename, boundary &box, config_type &config, int myrank, int nranks){\n"
            self._increment_indent()

            # Create the HDF5 file
            code = code + self.indent() + "hid_t acc_template = H5Pcreate(H5P_FILE_ACCESS);\n"
            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "MPI_Info info; MPI_Info_create(&info);\n"
                code = code + self.indent() + "H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);\n"
                code = code + self.indent() + "H5Pset_coll_metadata_write(acc_template, 1); //metadata writes are collective\n"
            code = code + self.indent() + "hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);\n"
            code = code + self.indent() + "if(file_id < 0){\n"
            self._increment_indent()
            code = code + self.indent() + "std::cout << \"[\" << myrank << \"] failed to open \" << filename << \"\\n\";\n"

            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "MPI_Abort(MPI_COMM_WORLD, 1);\n"
            else:
                code = code + self.indent() + "abort();\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n\n"


            # Get the local size, and local offset in the file.
            code = code + self.indent() + "hsize_t h_dims[1];\n"
            code = code + self.indent() + "h_dims[0] = particle_aosoa.size();\n"

            code = code + self.indent() + "int my_offset = 0;\n"

            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "if(myrank == 0 && nranks > 1){\n"
                self._increment_indent()
                code = code + self.indent() + "int npart = particle_aosoa.size();"
                code = code + self.indent() + "MPI_Send(&npart, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);\n"
                self._decrement_indent()
                code = code + self.indent() + "}else if (myrank == nranks-1){\n"
                self._increment_indent()
                code = code + self.indent() + "MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);\n"
                self._decrement_indent()
                code = code + self.indent() + "}else{\n"
                self._increment_indent()
                code = code + self.indent() + "MPI_Recv(&my_offset, 1, MPI_INT, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);\n"
                code = code + self.indent() + "int npart = my_offset + particle_aosoa.size();\n"
                code = code + self.indent() + "MPI_Send(&npart, 1, MPI_INT, myrank+1, 0, MPI_COMM_WORLD);\n"
                self._decrement_indent()
                code = code + self.indent() + "}\n"

                # More setup stuff
                code = code + self.indent() + "hid_t xf_id = H5Pcreate(H5P_DATASET_XFER);\n"
                code = code + self.indent() + "H5Pset_dxpl_mpio(xf_id, H5FD_MPIO_COLLECTIVE);\n"
            else:
                code = code + self.indent() + "hid_t xf_id = H5P_DEFAULT;\n"

            # Get global dimensions
            code = code + self.indent() + "int global_part_count = 0;\n"
            if HartreeParticleDSL.get_mpi():
                code = code + self.indent() + "MPI_Allreduce( &h_dims[0], &global_part_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);\n"
            else:
                code = code + self.indent() + "global_part_count = h_dims[0];\n"
            code = code + self.indent() + "hsize_t gpc = (hsize_t) global_part_count;\n"
            code = code + self.indent() + "hid_t global_dim = H5Screate_simple(1, &gpc, NULL);\n"

            code = code + self.indent() + "hsize_t offset[2];\n"
            code = code + self.indent() + "hsize_t count[2];\n"
            code = code + self.indent() + "offset[0] = my_offset;\n"
            code = code + self.indent() + "offset[1] = 0;\n"
            code = code + self.indent() + "count[0] = h_dims[0];\n"
            code = code + self.indent() + "count[1] = 0;\n"
            code = code + self.indent() + "H5Sselect_hyperslab(global_dim, H5S_SELECT_SET, offset, NULL, count, NULL);\n\n"

            code = code + self.indent() + "hid_t memspace = H5Screate_simple(1, h_dims, NULL);\n"

            code = code + self.indent() + "//Offset in memspace is 0\n"
            code = code + self.indent() + "offset[0] = 0;\n"
            code = code + self.indent() + "H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);\n"

            # Now we output the variables according to the offsets
            # Create dataset and output the requested values
            for key in self._outputs.keys():
                # Find the type of the variable from the particle structure
                part_elem = self._outputs[key]
                elem_type = None
                h5_type = None
                part_indexing = ""
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part_position"
                        part_indexing = ", 0"
                        elem_type = "double"
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part_position"
                        part_indexing = ", 1"
                        elem_type = "double"
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part_position"
                        part_indexing = ", 2"
                        elem_type = "double"
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 0"
                        elem_type = "double"
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 1"
                        elem_type = "double"
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 2"
                        elem_type = "double"
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = PHDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    raise UnsupportedCodeError(f"{part_elem} element not supported in HDF5 IO.")
                elem_type = Cabana_PIR._type_map[elem_type]
                code = code + self.indent() + "hid_t {0}_output_field = H5Dcreate2(file_id, \"{0}\", {1}, global_dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n".format(
                        key, h5_type)
                code = code + self.indent() + "{0}* {1}_output_array = ({0} *) malloc(sizeof({0}) * particle_aosoa.size());\n".format(elem_type, key)
                code = code + self.indent() + "for( int i = 0; i < particle_aosoa.size(); i++){\n"
                self._increment_indent()
                code = code + self.indent() + "auto part = particle_aosoa.getTuple(i);\n"
                code = code + self.indent() + "{0}_output_array[i] = Cabana::get<{1}>(part{2});\n".format(key, part_elem, part_indexing)
                self._decrement_indent()
                code = code + self.indent() + "}\n\n"
                code = code + self.indent() + "H5Dwrite({0}_output_field, {1}, memspace, global_dim, xf_id, {0}_output_array);\n".format(key, h5_type)
                code = code + self.indent() + "free({0}_output_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_output_field);\n".format(key)

            code = code + "\n\n"
            code = code + self.indent() + "H5Fclose(file_id);\n"
            self._decrement_indent()
            code = code + self.indent() + "}\n"
            # ============================= HDF5_OUTPUT ============================

        return code


    def call_get_box_size_pir(self, part_count: int, filename: str, current_indent=4) -> str:
        '''
        Returns the code required for the Cabana_PIR backend to use this IO module for reading
        in the box dimensions from the HDF5 file.

        :param int part_count: Unused particle count value.
        :param str filename: Filename to read the data from.
        :param int current_indent: The code current indentation level. Default is 4.

        :returns: The function call used to get the box size from the file.
        :rtype: str
        '''
        rval = " " * current_indent 
        rval = rval + f"get_box_size(config.config_host(0).space.box_dims, {filename});\n"
        return rval

    def call_input_cabana_pir(self, part_count, filename, current_indent=4):
        '''
        Returns the Cabana_PIR call required to use this IO module for input.
        :returns: The function call used for this module.
        :rtype: str
        '''
        indentation = " " * current_indent
        rval = "/* Create structures of size 1 to initialise, the HDF5 function will resize them */\n"
        rval = rval + indentation + f"Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", 1);\n"
        rval = rval + indentation + f"Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", 1);\n"
        rval = rval + indentation + f"hdf5_input<decltype(particle_aosoa), decltype(particle_aosoa_host)>(particle_aosoa, particle_aosoa_host, config, config.config_host(0).space.box_dims, {filename});\n"
        return rval

    def call_output_cabana_pir(self, part_count, filename, variable=None, current_indent=4, indentation=4):
        '''
        Returns the Cabana_PIR call required to use this IO module for output.

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        code = "{\n"
        current_indent = current_indent + indentation
        if variable is not None:
            code = code + current_indent * " " + "char filename[300];\n"
            if "\"" in filename:
                rep_filename = filename.replace("\"", "")
                code = code + current_indent * " " + "sprintf(filename, \"" + f"{rep_filename}%.4d.hdf5" + "\""+ f", {variable});\n"
            else:
                code = code + current_indent * " " + "sprintf(filename, \"" + f"{filename}%.4d.hdf5" + "\", " + f"{variable});\n"
        else:
            if "\"" in filename:
                code = code + current_indent * " " + "char filename[300]" + f" = {filename};\n"
            else:
                code = code + current_indent * " " + "char filename[300]" + f" = \"{filename}\";\n"
        code = code + current_indent * " " + "Cabana::deep_copy(particle_aosoa_host, particle_aosoa);\n"
        code = code + current_indent * " " + "int myrank, nranks;\n"
        if HartreeParticleDSL.get_mpi():
            code = code + current_indent * " " + "MPI_Comm_rank( MPI_COMM_WORLD, &myrank );\n"
            code = code + current_indent * " " + "MPI_Comm_size( MPI_COMM_WORLD, &nranks );\n"
        else:
            code = code + current_indent * " " + "myrank = 0;\n"
            code = code + current_indent * " " + "nranks = 1;\n"
        code = code + current_indent * " " + "hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, filename, config.config_host(0).space.box_dims, config, myrank, nranks);\n"
        current_indent = current_indent - indentation
        code = code + current_indent * " " + "}\n"

        return code
