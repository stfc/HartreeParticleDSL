from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.backends.C_AOS.C_AOS import C_AOS
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin
from HartreeParticleDSL.backends.FDPS_backend.FDPS import FDPS
from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin
from HartreeParticleDSL.backends.Cabana_backend.Cabana import Cabana
from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin

class HDF5_IO(IO_Module, C_AOS_IO_Mixin, FDPS_IO_Mixin, Cabana_IO_Mixin):
    '''Implementation of the HDF5 IO Module'''
    type_map = {c_int : "H5T_STD_I32LE",
                c_double : "H5T_NATIVE_DOUBLE",
                c_float : "H5T_NATIVE_FLOAT",
                c_int64_t : "H5T_STD_I64LE",
                c_int32_t : "H5T_STD_I32LE",
                c_int8_t : "H5T_STD_I8LE",
                c_bool : "H5T_NATIVE_HBOOL"}

    def __init__(self, indent=4):
        super().__init__()
        self._inputs = {}
        self._outputs = {}
        self._indent = indent
        self._current_indent = 0

    def add_input(self, hdf_input, particle_input):
        self._inputs[hdf_input] = particle_input

    def add_output(self, hdf_output, particle_output):
        self._outputs[hdf_output] = particle_output

    def indent(self):
        return " " * self._current_indent

    def increment_indent(self):
        self._current_indent = self._current_indent + self._indent

    def decrement_indent(self):
        self._current_indent = self._current_indent - self._indent

    # C_AOS_IO_Mixin functions

    def get_includes_c(self):
        '''
        :returns: The includes required for this IO module.
        :rtype: List of str

        '''
        includes = []
        includes.append("<stdlib.h>")
        includes.append("\"hdf5.h\"")
        return includes

    def gen_code_c(self, part_type):
        '''
        Returns the C code required for this IO module.

        :param part_type: The particle type used in the system.
        :type part_type: Particle
        '''
        code = ""

        if len(self._inputs) > 0:
            # Create the hdf5_input function
            code = code + self.indent() + "struct part* hdf5_input(const char* filename, struct config_type* config){\n"
            self.increment_indent()
            # Check that the asked for fields are in the file, and pull out the number of particles
            code = code + self.indent() + "hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Failed to open %s\\n\", filename);\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_test_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "if( {0}_test_var < 0)".format(key) + "{\n"
                self.increment_indent()
                code = code + self.indent() + "printf(\"Failed to find dataset {0}\\n\");\n".format(key)
                code = code + self.indent() + "exit(1);\n"
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "H5Dclose({0}_test_var);\n".format(key)
            # All of the fields are in the HDF5 file as requested, fetch the particle count.
            code = code + self.indent() + "hid_t part_count_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(list(self._inputs.keys())[0])
            code = code + self.indent() + "hid_t space = H5Dget_space(part_count_var);\n"
            code = code + self.indent() + "int ndims = H5Sget_simple_extent_ndims(space);\n"
            # Can't handle multi dimensional datasets for now
            code = code + self.indent() + "if( ndims != 1 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Don't yet support multidimensional datasets.\\n\");\n".format(key)
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            # Get the num parts
            code = code + self.indent() + "hsize_t dims[1];\n"
            code = code + self.indent() + "H5Sget_simple_extent_dims(space, dims, NULL);\n"
            code = code + self.indent() + "int num_parts = dims[0];\n"

            # Init the data structures
            code = code + self.indent() + "struct part *parts = (struct part*) malloc(sizeof(struct part) * num_parts);\n"
            code = code + self.indent() + "config->space.nparts = num_parts;\n\n"

            # Read the values from the file into the particles
            code = code + self.indent() + "hsize_t shape[2], offsets[2];\n"
            code = code + self.indent() + "int rank = 2;\n"
            code = code + self.indent() + "shape[0] = num_parts;\n"
            code = code + self.indent() + "shape[1] = 1;\n"
            code = code + self.indent() + "offsets[0] = 0;\n"
            code = code + self.indent() + "offsets[1] = 0;\n"
            code = code + self.indent() + "hid_t memspace = H5Screate_simple(rank, shape, NULL);\n"
            code = code + self.indent() + "hid_t filespace;\n"

            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_read_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "filespace = H5Dget_space({0}_read_var);\n".format(key)
                code = code + self.indent() + "H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
                # Find the type of the variable from the particle structure
                part_elem = self._inputs[key]
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part.position[0]"
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part.position[1]"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part.position[2]"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part.velocity[0]"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part.velocity[1]"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part.velocity[2]"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    assert False
                elem_type = C_AOS._type_map[elem_type]
                code = code + self.indent() + "{0}* {1}_temp_array = ({0}*) malloc(sizeof({0}) * num_parts);\n".format(elem_type, key)
                code = code + self.indent() + "H5Dread({0}_read_var, {1}, memspace, filespace, H5P_DEFAULT, {0}_temp_array);\n".format(key, h5_type)
                code = code + self.indent() + "for( int i = 0; i < num_parts; i++){\n"
                self.increment_indent()
                code = code + self.indent() + "parts[i].{0} = {1}_temp_array[i];\n".format(part_elem, key)
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "free({0}_temp_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_read_var);\n".format(key)

            # All data read in
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"

            code = code + self.indent() + "return parts;\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"

        if len(self._outputs) > 0:
            # Create the hdf5 output function
            code = code + "\n"
            code = code + self.indent() + "void hdf5_output(struct part* parts, struct config_type* config, const char* filename){\n"
            self.increment_indent()
            # Create the HDF5 file
            code = code + self.indent() + "hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Couldn't create HDF5 file\\n\");\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"
            # Setup dimensions of outputs
            code = code + self.indent() + "hsize_t dims[1];\n"
            code = code + self.indent() + "dims[0] = config->space.nparts;\n"
            code = code + self.indent() + "hid_t dim = H5Screate_simple(1, dims, NULL);\n"

            # Create dataset and output the requested values
            for key in self._outputs.keys():
                # Find the type of the variable from the particle structure
                part_elem = self._outputs[key]
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part.position[0]"
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part.position[1]"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part.position[2]"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part.velocity[0]"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part.velocity[1]"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part.velocity[2]"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    assert False
                elem_type = C_AOS._type_map[elem_type]
                code = code + self.indent() + "hid_t {0}_output_field = H5Dcreate2(file_id, \"{0}\", {1}, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n".format(
                        key, h5_type)
                # TODO Check if we successfully made the dataset
                code = code + self.indent() + "{0}* {1}_output_array = ({0} *) malloc(sizeof({0}) * config->space.nparts);\n".format(elem_type, key)
                code = code + self.indent() + "for( int i = 0; i < config->space.nparts; i++){\n"
                self.increment_indent()
                code = code + self.indent() + "{0}_output_array[i] = parts[i].{1};\n".format(key, part_elem)
                self.decrement_indent()
                code = code + self.indent() + "}\n\n"
                code = code + self.indent() + "H5Dwrite({0}_output_field, {1}, H5S_ALL, H5S_ALL, H5P_DEFAULT, {0}_output_array);\n".format(key, h5_type)
                code = code + self.indent() + "free({0}_output_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_output_field);\n".format(key)


            # End of output function
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"

        return code

    def call_input_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for input.

        :returns: The function call used for this module.
        :rtype: str
        '''

        return f"hdf5_input({filename}, config);"

    def call_output_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for output.

        :returns: The function call used for this module.
        :rtype: str
        '''

        return f"hdf5_output(parts, config, {filename});"

    # FDPS_IO_Mixin functions
    def gen_code_fdps(self, part_type):
        '''
        Returns the FDPS C++ code required for this IO module.

        :param part_type: The particle type used in the system.
        :type part_type: Particle
        '''
        code = ""

        if len(self._inputs) > 0:
            # Create the hdf5_input function
            code = code + self.indent() + "void hdf5_input(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char *filename){\n"
            self.increment_indent()
            code = code + self.indent() + "hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Failed to open %s\\n\", filename);\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_test_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "if( {0}_test_var < 0)".format(key) + "{\n"
                self.increment_indent()
                code = code + self.indent() + "printf(\"Failed to find dataset {0}\\n\");\n".format(key)
                code = code + self.indent() + "exit(1);\n"
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "H5Dclose({0}_test_var);\n".format(key)
            # All of the fields are in the HDF5 file as requested, fetch the particle count.
            code = code + self.indent() + "hid_t part_count_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(list(self._inputs.keys())[0])
            code = code + self.indent() + "hid_t space = H5Dget_space(part_count_var);\n"
            code = code + self.indent() + "int ndims = H5Sget_simple_extent_ndims(space);\n"
            # Can't handle multi dimensional datasets for now
            code = code + self.indent() + "if( ndims != 1 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Don't yet support multidimensional datasets.\\n\");\n".format(key)
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            # Get the num parts
            code = code + self.indent() + "hsize_t[1] dims;\n"
            code = code + self.indent() + "H5Sget_simple_extent_dims(space, dims, NULL);\n"
            code = code + self.indent() + "int num_parts = dims[0];\n"
       
            # Initialise the particle system
            code = code + self.indent() + "particle_system.initialize();\n"
            # Single node for now.
            code = code + self.indent() + "particle_system.setNumberOfParticleLocal(num_parts);\n"
            code = code + self.indent() + "config.space.nparts = num_parts;\n"
            # Setup the domain
            code = code + self.indent() + "PS::DomainInfo dinfo;\n"
            # Assume box sizes is setup for now.
            code = code + self.indent() + "dinfo.initialize();\n"
            # Assume periodic only for now.
            code = code + self.indent() + "dinfo.setBoundaryCondition ( PS::BOUNDARY_CONDITION_PERIODIC_XYZ );\n"
            code = code + self.indent() + \
            '''dinfo.setPosRootDomain(PS::F64vec(config.space.box_dims.x_min,
                                         config.space.box_dims.y_min,
                                         config.space.box_dims.z_min),
                              PS::F64vec(config.space.box_dims.x_max,
                                         config.space.box_dims.y_max,
                                         config.space.box_dims.z_max));\n\n'''

            # Read the values from the file into the particles
            code = code + self.indent() + "hsize_t shape[2], offsets[2];\n"
            code = code + self.indent() + "int rank = 2;\n"
            code = code + self.indent() + "shape[0] = num_parts;\n"
            code = code + self.indent() + "shape[1] = 1;\n"
            code = code + self.indent() + "offsets[0] = 0;\n"
            code = code + self.indent() + "offsets[1] = 0;\n"
            code = code + self.indent() + "hid_t memspace = H5Screate_simple(rank, shape, NULL);\n"
            code = code + self.indent() + "hid_t filespace;\n"
            code = code + self.indent() + "H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_read_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "filespace = H5Dget_space({0}_read_var);\n".format(key)
                code = code + self.indent() + "H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
                # Find the type of the variable from the particle structure
                part_elem = self._inputs[key]
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part.position.x"
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part.position.y"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part.position.z"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part.velocity[0]"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part.velocity[1]"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part.velocity[2]"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    assert False
                elem_type = FDPS._type_map[elem_type]
                code = code + self.indent() + "{0}* {1}_temp_array = ({0}*) malloc(sizeof({0}) * num_parts);\n".format(elem_type, key)
                code = code + self.indent() + "H5Dread({0}_read_var, {1}, memspace, filespace, H5P_DEFAULT, {0}_temp_array);\n".format(key, h5_type)
                code = code + self.indent() + "for( int i = 0; i < num_parts; i++){\n"
                self.increment_indent()
                code = code + self.indent() + "particle_system[i].{0} = {1}_temp_array[i];\n".format(part_elem, key)
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "free({0}_temp_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_read_var);\n".format(key)

            # All data read in
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"

        if len(self._outputs) > 0:
            # Create the hdf5 output function
            code = code + "\n"
            code = code + self.indent() + "void hdf5_output(PS::ParticleSystem<FullParticle>& particle_system, config_type& config, const char* filename){\n"
            self.increment_indent()
            # Create the HDF5 file
            code = code + self.indent() + "hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Couldn't create HDF5 file\\n\");\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"
            # Setup dimensions of outputs
            code = code + self.indent() + "hsize_t dims[1];\n"
            code = code + self.indent() + "dims[0] = config.space.nparts;\n"
            code = code + self.indent() + "hid_t dim = H5Screate_simple(1, dims, NULL);\n"

            # Create dataset and output the requested values
            for key in self._outputs.keys():
                # Find the type of the variable from the particle structure
                part_elem = self._outputs[key]
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part.position.x"
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part.position.y"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part.position.z"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part.velocity[0]"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part.velocity[1]"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part.velocity[2]"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    assert False
                elem_type = FDPS._type_map[elem_type]
                code = code + self.indent() + "hid_t {0}_output_field = H5Dcreate2(file_id, \"{0}\", {1}, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n".format(
                        key, h5_type)
                # TODO Check if we successfully made the dataset
                code = code + self.indent() + "{0}* {1}_output_array = ({0} *) malloc(sizeof({0}) * config.space.nparts);\n".format(elem_type, key)
                code = code + self.indent() + "for( int i = 0; i < config.space.nparts; i++){\n"
                self.increment_indent()
                code = code + self.indent() + "{0}_output_array[i] = particle_system[i].{1};\n".format(key, part_elem)
                self.decrement_indent()
                code = code + self.indent() + "}\n\n"
                code = code + self.indent() + "H5Dwrite({0}_output_field, {1}, H5S_ALL, H5S_ALL, H5P_DEFAULT, {0}_output_array);\n".format(key, h5_type)
                code = code + self.indent() + "free({0}_output_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_output_field);\n".format(key)

            # End of output function
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"
        return code

    def call_input_fdps(self, part_count, filename, current_indent=4):
        '''
        Returns the FDPS call required to use this IO module for input.

        :returns: The function call used for this module.
        :rtype: str
        '''
        return f"hdf5_input(particle_system, config, {filename});"

    def call_output_fdps(self, part_count, filename):
        '''
        Returns the FDPS call required to use this IO module for output.

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        return f"hdf5_output(particle_system, config, {filename});"

    def get_includes_fdps(self):
        '''
        Returns the includes required to use this IO module for FDPS.

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"hdf5.h\"")
        return includes

    def gen_code_cabana(self, part_type):
        '''
        Returns the Cabana code required for this IO module.

        :param part_type: The particle type used in the system.
        :type part_type: Particle
        '''
        code = ""
        if len(self._inputs) > 0:
            # Create the hdf5_input function
            code = code + self.indent() + "template <class aosoa_class, class aosoa_host_class> void hdf5_input(aosoa_class particle_system, aosoa_host_class parts_host, config_type& config, const char *filename){\n"
            self.increment_indent()
            code = code + self.indent() + "hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Failed to open %s\\n\", filename);\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_test_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "if( {0}_test_var < 0)".format(key) + "{\n"
                self.increment_indent()
                code = code + self.indent() + "printf(\"Failed to find dataset {0}\\n\");\n".format(key)
                code = code + self.indent() + "exit(1);\n"
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "H5Dclose({0}_test_var);\n".format(key)
            # All of the fields are in the HDF5 file as requested, fetch the particle count.
            code = code + self.indent() + "hid_t part_count_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(list(self._inputs.keys())[0])
            code = code + self.indent() + "hid_t space = H5Dget_space(part_count_var);\n"
            code = code + self.indent() + "int ndims = H5Sget_simple_extent_ndims(space);\n"
            # Can't handle multi dimensional datasets for now
            code = code + self.indent() + "if( ndims != 1 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Don't yet support multidimensional datasets.\\n\");\n".format(key)
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n"
            # Get the num parts
            code = code + self.indent() + "hsize_t[1] dims;\n"
            code = code + self.indent() + "H5Sget_simple_extent_dims(space, dims, NULL);\n"
            code = code + self.indent() + "int num_parts = dims[0];\n"
            # Resize the AoSoA
            code = code + self.indent() + "int new_size = static_cast<int>(num_parts);\n"
            code = code + self.indent() + "particle_system.resize(new_size);\n"
            code = code + self.indent() + "parts_host.resize(new_size);\n"
           
            # Read the values from the file into the particles
            code = code + self.indent() + "hsize_t shape[2], offsets[2];\n"
            code = code + self.indent() + "int rank = 2;\n"
            code = code + self.indent() + "shape[0] = num_parts;\n"
            code = code + self.indent() + "shape[1] = 1;\n"
            code = code + self.indent() + "offsets[0] = 0;\n"
            code = code + self.indent() + "offsets[1] = 0;\n"
            code = code + self.indent() + "hid_t memspace = H5Screate_simple(rank, shape, NULL);\n"
            code = code + self.indent() + "hid_t filespace;\n"
            code = code + self.indent() + "H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
            for key in self._inputs.keys():
                code = code + self.indent() + "hid_t {0}_read_var = H5Dopen2(file_id, \"{0}\", H5P_DEFAULT);\n".format(key)
                code = code + self.indent() + "filespace = H5Dget_space({0}_read_var);\n".format(key)
                code = code + self.indent() + "H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);\n"
                # Find the type of the variable from the particle structure
                part_elem = self._inputs[key]
                part_indexing = ""
                elem_type = None
                h5_type = None
                if part_elem.startswith("core_part"):
                    temp_string = part_elem.replace("core_part.", "")
                    if temp_string == "position.x" or temp_string == "position[0]":
                        part_elem = "core_part_position"
                        part_indexing = ", 0"
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part_position"
                        part_indexing = ", 1"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part_position"
                        part_indexing = ", 2"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 0"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 1"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 2"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    assert False
                elem_type = Cabana._type_map[elem_type]
                code = code + self.indent() + "auto {0}_slice = Cabana::slice<{1}>(parts_host);\n".format(key, part_elem)
                code = code + self.indent() + "{0}* {1}_temp_array = ({0}*) malloc(sizeof({0}) * num_parts);\n".format(elem_type, key)
                code = code + self.indent() + "H5Dread({0}_read_var, {1}, memspace, filespace, H5P_DEFAULT, {0}_temp_array);\n".format(key, h5_type)
                code = code + self.indent() + "for( int i = 0; i < num_parts; i++){\n"
                self.increment_indent()
                code = code + self.indent() + "{0}_slice(i{2}) = {1}_temp_array[i];\n".format(key, key, part_indexing)
                self.decrement_indent()
                code = code + self.indent() + "}\n"
                code = code + self.indent() + "free({0}_temp_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_read_var);\n".format(key)

            # All data read in
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"
            code = code + self.indent() + "Cabana::deep_copy(particle_system, parts_host);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"

        if len(self._outputs) > 0:
            # Create the hdf5 output function
            code = code + "\n"
            code = code + self.indent() + "template <class aosoa_class> void hdf5_output(aosoa_class particle_aosoa, config_type& config, const char* filename){\n"
            self.increment_indent()
            # Create the HDF5 file
            code = code + self.indent() + "hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);\n"
            code = code + self.indent() + "if( file_id < 0 ){\n"
            self.increment_indent()
            code = code + self.indent() + "printf(\"Couldn't create HDF5 file\\n\");\n"
            code = code + self.indent() + "exit(1);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"
            # Setup dimensions of outputs
            code = code + self.indent() + "hsize_t dims[1];\n"
            code = code + self.indent() + "dims[0] = particle_aosoa.size();\n"
            code = code + self.indent() + "hid_t dim = H5Screate_simple(1, dims, NULL);\n"

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
                        elem_type = c_double
                    elif temp_string == "position.y" or temp_string == "position[1]":
                        part_elem = "core_part_position"
                        part_indexing = ", 1"
                        elem_type = c_double
                    elif temp_string == "position.z" or temp_string == "position[2]":
                        part_elem = "core_part_position"
                        part_indexing = ", 2"
                        elem_type = c_double
                    elif temp_string == "velocity.x" or temp_string == "velocity[0]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 0"
                        elem_type = c_double
                    elif temp_string == "velocity.y" or temp_string == "velocity[1]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 1"
                        elem_type = c_double
                    elif temp_string == "velocity.z" or temp_string == "velocity[2]":
                        part_elem = "core_part_velocity"
                        part_indexing = ", 2"
                        elem_type = c_double
                elif part_elem.startswith("neighbour_part"):
                    pass
                else:
                    elem_type = part_type.particle_type[part_elem]['type']
                h5_type = HDF5_IO.type_map.get(elem_type, None)
                if h5_type is None:
                    print(part_elem, elem_type, h5_type)
                    assert False
                elem_type = Cabana._type_map[elem_type]
                code = code + self.indent() + "hid_t {0}_output_field = H5Dcreate2(file_id, \"{0}\", {1}, dim, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n".format(
                        key, h5_type)
                # TODO Check if we successfully made the dataset
                # FIXME Access to elements is different for Cabana
                code = code + self.indent() + "{0}* {1}_output_array = ({0} *) malloc(sizeof({0}) * particle_aosoa.size());\n".format(elem_type, key)
                code = code + self.indent() + "for( int i = 0; i < particle_aosoa.size(); i++){\n"
                self.increment_indent()
                code = code + self.indent() + "auto part = particle_aosoa.getTuple(i);\n"
                code = code + self.indent() + "{0}_output_array[i] = Cabana::get<{1}>(part{2});\n".format(key, part_elem, part_indexing)
                self.decrement_indent()
                code = code + self.indent() + "}\n\n"
                code = code + self.indent() + "H5Dwrite({0}_output_field, {1}, H5S_ALL, H5S_ALL, H5P_DEFAULT, {0}_output_array);\n".format(key, h5_type)
                code = code + self.indent() + "free({0}_output_array);\n".format(key)
                code = code + self.indent() + "H5Dclose({0}_output_field);\n".format(key)

            # End of output function
            code = code + "\n" + self.indent() + "H5Fclose(file_id);\n"
            self.decrement_indent()
            code = code + self.indent() + "}\n\n"
        return code

    def call_input_cabana(self, part_count, filename, current_indent=4):
        '''
        Returns the Cabana call required to use this IO module for input.
        NYI
        :returns: The function call used for this module.
        :rtype: str
        '''
        indentation = " " * current_indent
        rval = "/* Create structures of size 1 to initialise, the HDF5 function will resize them */\n"
        rval = rval + indentation + f"Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", 1);\n"
        rval = rval + indentation + f"Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", 1);\n"
        rval = rval + indentation + f"hdf5_input<decltype(particle_aosoa), decltype(particle_aosoa_host)>(particle_aosoa, particle_aosoa_host, config, {filename});\n"

        return ""

    def call_output_cabana(self, part_count, filename, variable=None):
        '''
        Returns the Cabana call required to use this IO module for output.

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        code = "{\n"
        if variable is not None:
            code = code + "char filename[300];\n"
            code = code + "        sprintf(filename, \"" + f"{filename}%.4d.hdf5" + "\", " + f"{variable});\n        "
        else:
            code = code + "char filename[300]" + f" = \"{filename}\";\n"
        code = code + "Cabana::deep_copy(particle_aosoa_host, particle_aosoa);\n"
        code = code + f"        hdf5_output<decltype(particle_aosoa_host)>(particle_aosoa_host, config, filename);"
        code = code + "        }\n"
        return code

    def get_includes_cabana(self):
        '''
        Returns the includes required to use this IO module for Cabana.

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"hdf5.h\"")
        return includes
