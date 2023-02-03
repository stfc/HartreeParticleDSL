from __future__ import annotations

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
from HartreeParticleDSL.backends.Cabana_PIR_backend.Cabana_PIR_IO_Mixin import Cabana_PIR_IO_Mixin
from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import Cabana_PIR_Visitor
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, DataType
from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol


class Cabana_PIR(Backend):
    '''
    Cabana Backend class using PIR.
    Outputs from this class use the Cabana & Kokkos system to compute
    parallel particle simulations.
    '''
    _type_map = {c_int : "int",
            "int" : "int",
            c_double : "double", # Is this used?
            "double" : "double",
            c_float : "float",
            "float" : "float",
            c_int64_t : "int64_t",
            "int64_t" : "int64_t",
            c_int32_t : "int32_t",
            "int32_t" : "int32_t",
            c_int8_t : "int8_t",
            "int8_t" : "int8_t",
            c_bool : "bool",
            "bool" : "bool"
            }
    _type_sizes = {_type_map[c_int] : 4,
            _type_map[c_double] : 8,
            _type_map[c_float] : 4,
            _type_map[c_int64_t] : 8,
            _type_map[c_int32_t] : 4,
            _type_map[c_int8_t] : 1,
            _type_map[c_bool] : 4
            }


    # Variables used to manage where the cutoff radius is stored
    CONSTANT = 0
    PARTICLE = 1


    def __init__(self) -> None:
        self._includes = []
        self._includes.append("<Cabana_Core.hpp>")
        self._includes.append("<Kokkos_Core.hpp>")
        self._includes.append("<Kokkos_Random.hpp>")
        self._includes.append("<iostream>")
        self._includes.append("<cmath>")
        self._includes.append("<cstdio>")
        self._includes.append("<math.h>")
        self._includes.append("\"part.h\"")

        self._input_module = None
        self._output_module = None
        self._coupled_systems = []

        # We need to keep track of what kernels use what slices
        self._kernel_slices = {}
        self._kernels = {}
        self._kernels_pir = {}

        # The Cabana backend needs to store the particle type reference
        self._particle = None #TODO Remove & pull from PIR

        self._structures = {}

        self._global_values = {}
        self._boundary_condition = None
        self._boundary_condition_tree = None

        self._current_kernel = None

        self._require_random = False

    def set_boundary_condition(self, boundary_condition_kernel) -> None:
        '''
        Set the boundary condition kernel to use with this backend. This kernel
        will be automatically run whenever a kernel which updates particle
        positions is invoked. If you do not want this behaviour then do not
        use this function.

        :param boundary_condition_kernel: The kernel to use as a boundary condition.
        :type boundary_contidion_kernel: A perpart kernel - usually denoted with \
                @kernels.perpart_interaction .
        '''
        if not isinstance(boundary_condition_kernel, kernels.perpart_kernel_wrapper):
            raise TypeError("Cannot set boundary condition to a non perpart kernel.")

        self._boundary_condition = boundary_condition_kernel
        self._boundary_condition_tree = boundary_condition_kernel.get_kernel_tree()

    @property
    def boundary_condition(self) -> PerPartKernel:
        '''
        :returns: the Particle IR PerPartKernel representation of the current 
        boundary condition for this backend.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kernels.PerPartKernel` or None

        '''
        from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import \
                pir_perpart_visitor
        if self._boundary_condition is None:
            return None
        conv = pir_perpart_visitor()
        return conv.visit(self._boundary_condition_tree)

    @property
    def structures(self):
        '''
        :returns: the structures currently used in this backend by coupled systems.
        :rtype: Dict of str, StructureType pairs.
        '''
        return self._structures

    def register_kernel(self, kernel_name: str, kernel: Kern):
        '''
        Registers a kernel with the backend for code generation,
        and use in invoke statements.

        :param str kernel_name: The name of the kernel to register.
        :param kernel: The Kern object corresponding to the named kernel.
        :type kernel: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kern.Kern`

        :raises InvalidNameError: If a kernel with the kernel_name already exists.
        '''
        if kernel_name in self._kernels.keys():
            raise InvalidNameError(f"Kernel with name {kernel_name} already exists.") 
        self._kernels[kernel_name] = kernel


    def add_structure(self, structure_type: StructureType, structure_name: str):
        '''
        Adds a structure to the backend for code generation.

        :param structure_type: The StructureType representing the structure.
        :type structure_type: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.StructureType`
        :param str structure_name: The name of the structure.
        '''
        # Should this check the structure type is defined?
        self._structures[structure_name] = structure_type

    def add_kernel_slices(self, kernel_name : str, kernel_slice: str):
        '''
        Adds access to a particle slice for a specific kernel.

        :param str kernel_name: The name of the kernel to for the kernel_slice \
                                to be generated for.
        :param str kernel_slice: The slice of the particle arrays required for \
                                 the kernel given by kernel_name.
        '''
        if kernel_name not in self._kernel_slices.keys():
            self._kernel_slices[kernel_name] = []
        if kernel_slice not in self._kernel_slices[kernel_name]:
            self._kernel_slices[kernel_name].append(kernel_slice)

    def add_type(self, type_name: str, the_type: DataType) -> None:
        '''
        Adds a new type to Cabana & PIR.

        :param str type_name: The name of the type used in the code, e.g. int \
                              is the name of the C integer type.
        :param the_type: The Particle_IR datatype represented by the type_name.
        :type type_type: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.DataType`

        :raises UnsupportedTypeError: if the_type is not a PIR DataType object.
        :raises InvalidNameError: if the type name is already present in PIR
        '''
        if type_mapping_str.get(type_name) is not None:
            raise InvalidNameError(f"Attempting to add new type {type_name} "
                                   f"but a type with that name already exists.")
        if not isinstance(the_type, DataType):
            raise UnsupportedTypeError(
                    f"Attempting to add a new type but got {type(the_type)} "
                    "but expected a Particle IR DataType instead.")

        type_mapping_str[type_name] = the_type

    def create_global_variable(self, c_type: DataType, name: str, initial_value: Union[str, NoneType]=None):
        '''
        Function to create a global variable in the header file.

        :param c_type: The c_type of this variable.
        :type c_type: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.DataType.DataType`
        :param str name: The name of the variable.
        :param initial_value: The initial value of this variable
        :type initial_value: str or None.
        '''
        from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
        from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol
        from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol
        from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
        from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType, \
                PointerType, ArrayType, ScalarType
        if isinstance(c_type, StructureType):
            HartreeParticleDSL.global_symbol_table().new_symbol(name, c_type, StructureSymbol)
            self._global_values[name] = None
        elif isinstance(c_type, PointerType):
            HartreeParticleDSL.global_symbol_table().new_symbol(name, c_type, PointerSymbol)
            self._global_values[name] = initial_value
        elif isinstance(c_type, ArrayType):
            HartreeParticleDSL.global_symbol_table().new_symbol(name, c_type, ArraySymbol)
            self._global_values[name] = initial_value
        elif isinstance(c_type, ScalarType):
            HartreeParticleDSL.global_symbol_table().new_symbol(name, c_type, ScalarTypeSymbol)
            self._global_values[name] = initial_value
        else:
            raise TypeError("Attempting to create global variable but c_type argument is not "
                    f"a supported datatype. Got {type(c_type)}")


    def set_io_modules(self, input_module, output_module):
        '''
        Function to set the IO Modules. Called from HartreeParticleDSL
        helper functions instead of directly by the user.

        :param input_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.Cabana_backend.XCXX` \
                            or None
        :param ouput_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.Cabana_backend.XCXX` \
                             or None
        :raises InvalidIOModuleError: if either input_module or output_module \
                                      are not implementing the XCXXXX \
                                      class
        '''
        # TODO
        if input_module is not None and not isinstance(input_module, Cabana_PIR_IO_Mixin):
            raise InvalidIOModuleError(f"{self.__class__.__name__} backend does not support "
                                       f"{input_module.__class__.__name__} IO "
                                       f"module at this time.")
        if output_module is not None and not isinstance(output_module, Cabana_PIR_IO_Mixin):
            raise InvalidIOModuleError(f"{self.__class__.__name__} backend does not support "
                                       f"{output_module.__class__.__name__} IO "
                                       f"module at this time.")
        self._input_module = input_module
        self._output_module = output_module

    def add_include(self, include_string):
        '''
        Function to add an include to the Cabana include list.

        :param include_string: The string to add to the includes, e.g.
                               "\"myfile.h\""
        :type include_string: str
        '''
        if include_string not in self._includes:
            self._includes.append(include_string)

    def generate_includes(self):
        '''
        Generates the includes required for computation of the system.
        Used during code generation to ensure correct headers are included
        to ensure code compiles correctly.

        :returns: String containing the includes required.
        :rtype: str
        '''
        if self._input_module is not None:
            input_includes = self._input_module.get_includes_cabana_pir() #FIXME
            for include in input_includes:
                if include not in self._includes:
                    self._includes.append(include)
        if self._output_module is not None:
            output_includes = self._output_module.get_includes_cabana_pir() #FIXME
            for include in output_includes:
                if include not in self._includes:
                    self._includes.append(include)

        include_string = ""
        for include in self._includes:
            include_string = include_string + f"#include {include}\n"
        return include_string

    def mpi_headers(self, config, particle):
        '''
        Generates the extra header file code required for this backend
        if MPI is enabled.
        '''
        #Templated class on aosoa type
        # Needs a space for each member of the particle, e.g.:
#        template<class aosoa> class Migrator{
#
#    private:
        # Setup indentation
        current_indent = 0
        indent = 4

        def get_indent():
            nonlocal current_indent
            return current_indent * " "

        def inc_indent():
            nonlocal current_indent
            nonlocal indent
            current_indent = current_indent + indent

        def dec_indent():
            nonlocal current_indent
            nonlocal indent
            current_indent = current_indent - indent

        rval = "\n\ntemplate<class aosoa> class Migrator{\n"
        inc_indent()
        rval = rval + "\n"
        rval = rval + get_indent() + "private:\n"
        inc_indent()
        element_info = {} # element name: [datatype, array_dimension (i.e. 0-N), array_indices]
        for element in particle.particle_type:
            if element == "core_part" or element == "neighbour_part":
                continue
            c_type = particle.particle_type[element]['type']
            is_array = particle.particle_type[element]['is_array']
            if is_array:
                x = c_type.index("[")
                dimensionality = 0
                array_indices = []
                array_str = c_type[x:]
                while array_str.find("[") >= 0:
                    dimensionality += 1
                    val = int(array_str[array_str.index("[")+1:array_str.index("]")])
                    array_indices.append(val)
                    array_str = array_str[array_str.index("]")+1:]
                element_info[element] = [c_type[:x], dimensionality, array_indices]
            else:
                element_info[element] = [c_type, 0, []]
#        Kokkos::View<int**, MemorySpace> id_space;
#        Kokkos::View<double**, MemorySpace> weight_space;
#        Kokkos::View<double**, MemorySpace> mass_space;
#        Kokkos::View<double**, MemorySpace> charge_space;
#        Kokkos::View<double**, MemorySpace> pos_space;
#        Kokkos::View<double**, MemorySpace> px_space;
#        Kokkos::View<double**, MemorySpace> py_space;
#        Kokkos::View<double**, MemorySpace> pz_space;
#        int _buffer_size;
        for element in element_info.keys():
            rval = rval + get_indent() + "Kokkos::View<"
            rval = rval + element_info[element][0]
            rval = rval + "**" + element_info[element][1]*"*"
            rval = rval + ", MemorySpace> " + element+ "_space;\n"
        # Add position, velocity and neighbour spaces separately
        rval = rval + get_indent() + "Kokkos::View<double***, MemorySpace> pos_space;\n"
        rval = rval + get_indent() + "Kokkos::View<double***, MemorySpace> vel_space;\n"
        rval = rval + get_indent() + "Kokkos::View<double**, MemorySpace> cutoff_space;\n"
        rval = rval + get_indent() + "int _buffer_size;\n"
        rval = rval + "\n"
        dec_indent()



#
#    public:
#        Migrator(int buffer_size, int nr_neighbours){
#            _buffer_size = buffer_size;
#            id_space = Kokkos::View<int**, MemorySpace>("temp_id", nr_neighbours, buffer_size);
#            weight_space = Kokkos::View<double**, MemorySpace>("temp_weight", nr_neighbours, buffer_size);
#            mass_space = Kokkos::View<double**, MemorySpace>("temp_mass", nr_neighbours, buffer_size);
#            charge_space = Kokkos::View<double**, MemorySpace>("temp_charge", nr_neighbours, buffer_size);
#            pos_space = Kokkos::View<double**, MemorySpace>("temp_pos", nr_neighbours, buffer_size);
#            px_space = Kokkos::View<double**, MemorySpace>("temp_px", nr_neighbours, buffer_size);
#            py_space = Kokkos::View<double**, MemorySpace>("temp_py", nr_neighbours, buffer_size);
#            pz_space = Kokkos::View<double**, MemorySpace>("temp_pz", nr_neighbours, buffer_size);
#        }
        rval = rval + get_indent() + "public:\n"
        inc_indent()
        rval = rval + get_indent() + "Migrator(int buffer_size, int nr_neighbours){\n"
        inc_indent()
        rval = rval + get_indent() + "_buffer_size = buffer_size;\n"
        for element in element_info.keys():
            rval = rval + get_indent() + element + "_space = "
            rval = rval + "Kokkos::View<" + element_info[element][0]
            rval = rval + "**" + element_info[element][1]*"*"
            rval = rval + ", MemorySpace>(\"" + element + "_id\", nr_neighbours, buffer_size"
            #Extra indices?
            for index in element_info[element][2]:
                rval = rval + f", {index}"
            rval = rval + ");\n"
        rval = rval + get_indent() + "pos_space = Kokkos::View<double***, MemorySpace>(\"temp_pos\", nr_neighbours, buffer_size, 3);\n"
        rval = rval + get_indent() + "vel_space = Kokkos::View<double***, MemorySpace>(\"temp_velocity\", nr_neighbours, buffer_size, 3);\n"
        rval = rval + get_indent() + "cutoff_space = Kokkos::View<double**, MemorySpace>(\"temp_cutoff\", nr_neighbours, buffer_size);\n"
        dec_indent()
        rval = rval + get_indent() + "}\n\n"
        dec_indent()
#       Final function is exchange_data call.
        rval = rval + get_indent() + "void exchange_data( aosoa &particle_aosoa, std::vector<int> neighbors, int myrank, int npart){\n\n"
        inc_indent()
        rval = rval + get_indent() + "auto rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, \"rank\");\n"
        rval = rval + get_indent() + "auto last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, \"last_pos\");\n"
        rval = rval + get_indent() + "auto pos_s = Cabana::slice<core_part_position>(particle_aosoa, \"position\");\n"
        rval = rval + get_indent() + "auto vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, \"velocity\");\n"
        rval = rval + get_indent() + "auto cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, \"cutoff\");\n"
        for element in element_info.keys():
            rval = rval + get_indent() + "auto " + element + "_s = Cabana::slice<" + element + ">(particle_aosoa, \"" + element + "\");\n"

        rval = rval + get_indent() + "int *send_count = (int*) malloc(sizeof(int) * neighbors.size());\n"
        rval = rval + get_indent() + "int count_neighbours = 0;\n"
        rval = rval + get_indent() + "int end = particle_aosoa.size() - 1;\n"
        rval = rval + get_indent() + "for(int i = 0; i < neighbors.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "    send_count[i] = 0;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n\n"
        rval = rval + get_indent() + "for(int i = particle_aosoa.size()-1; i>=0; i--){\n"
        inc_indent()
        rval = rval + get_indent() + "if(rank_slice(i) != myrank && rank_slice(i) >= 0){\n"
        inc_indent()
        rval = rval + get_indent() + "int therank = rank_slice(i);\n"
        rval = rval + get_indent() + "for(int k = 0; k < neighbors.size(); k++){\n"
        inc_indent()
        rval = rval + get_indent() + "if(therank == neighbors[k]){\n"
        inc_indent()
        rval = rval + get_indent() + "therank = k;\n"
        rval = rval + get_indent() + "break;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        dec_indent()
        rval = rval + get_indent() + "int pos = send_count[therank];\n"

        # Fill in the slice data
        rval = rval + get_indent() + "pos_space(therank, pos, 0) = pos_s(i, 0);\n"
        rval = rval + get_indent() + "pos_space(therank, pos, 1) = pos_s(i, 1);\n"
        rval = rval + get_indent() + "pos_space(therank, pos, 2) = pos_s(i, 2);\n"
        rval = rval + get_indent() + "vel_space(therank, pos, 0) = vel_s(i, 0);\n"
        rval = rval + get_indent() + "vel_space(therank, pos, 1) = vel_s(i, 1);\n"
        rval = rval + get_indent() + "vel_space(therank, pos, 2) = vel_s(i, 2);\n"
        rval = rval + get_indent() + "cutoff_space(therank, pos) = cutoff_s(i);\n"
        
        for element in element_info.keys():
            if element_info[element][1] == 0:
                rval = rval + get_indent() + element + "_space(therank, pos) = " + element + "_s(i);\n"
            else:
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + get_indent() + "for(int " + indexer + "=0; " + indexer + f" < {element_info[element][2][index]}; "
                    rval = rval + indexer + "++){\n"
                    inc_indent()
                rval = rval + get_indent() + element + "_space(therank, pos"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ") = " + element + "_s(i"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ");\n"
                for index in range(element_info[element][1]):
                    dec_indent()
                    rval = rval + get_indent() + "}\n"

        rval = rval + get_indent() + "send_count[therank]++;\n\n"
        # Slice data should be filled now.
        # Move from end
        rval = rval + get_indent() + "while(rank_slice(end) != myrank && end > 0){\n"
        inc_indent()
        rval = rval + get_indent() + "end--;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "if(end > i){\n"
        inc_indent()
        rval = rval + get_indent() + "rank_slice(i) = rank_slice(end);\n"
        rval = rval + get_indent() + "pos_s(i, 0) = pos_s(end, 0);\n"
        rval = rval + get_indent() + "pos_s(i, 1) = pos_s(end, 1);\n"
        rval = rval + get_indent() + "pos_s(i, 2) = pos_s(end, 2);\n"
        rval = rval + get_indent() + "vel_s(i, 0) = vel_s(end, 0);\n"
        rval = rval + get_indent() + "vel_s(i, 1) = vel_s(end, 1);\n"
        rval = rval + get_indent() + "vel_s(i, 2) = vel_s(end, 2);\n"
        rval = rval + get_indent() + "cutoff_s(i) = cutoff_s(end);\n"
        for element in element_info.keys():
            if element_info[element][1] == 0:
                rval = rval + get_indent() + element + "_s(i) = " + element + "_s(end);\n"
            else:
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + get_indent() + "for(int " + indexer + "=0; " + indexer + f" < {element_info[element][2][index]}; "
                    rval = rval + indexer + "++){\n"
                    inc_indent()
                rval = rval + get_indent() + element + "_s(i"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ") = " + element + "_s(end"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ");\n"
                for index in range(element_info[element][1]):
                    dec_indent()
                    rval = rval + get_indent() + "}\n"

        rval = rval + get_indent() + "rank_slice(end) = -1;\n"
        dec_indent()
        rval = rval + get_indent() + "}else{\n"
        inc_indent()
        rval = rval + get_indent() + "rank_slice(i) = -1;\n"
        rval = rval + get_indent() + "end++;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "continue;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n\n"

#                    // If we've moved too far from the boundary we can stop.
#                if( i < end && (part_pos_s(i, 0) < (region_max - ( 2.0 * max_movement + sorting_size))) &&
#                    (part_pos_s(i, 0) > (region_min + (2.0 * max_movement + 2.0*sorting_size)))){
#                        break;
#                    }
        dec_indent()
        rval = rval + get_indent() + "}\n"



        # At this stage data is collected, so send the data.
        rval = rval + get_indent() + "// Data collected, need to send information to neighbours to know what to expect\n"      
        rval = rval + get_indent() + "int *recv_count = (int*) malloc(sizeof(int) * neighbors.size());\n"
        rval = rval + get_indent() + "MPI_Request *requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2);\n"
        rval = rval + get_indent() + "int req_num = 0;\n"
        rval = rval + get_indent() + "for(int i = 0; i < neighbors.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "recv_count[i] = 0;\n"
        rval = rval + get_indent() + "if(neighbors[i] == myrank){\n"
        inc_indent()
        rval = rval + get_indent() + "continue;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "MPI_Irecv(&recv_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);\n"
        rval = rval + get_indent() + "MPI_Isend(&send_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);\n"
        rval = rval + get_indent() + "MPI_Barrier(MPI_COMM_WORLD);\n"
        rval = rval + get_indent() + "int total_size = 0;\n"
        rval = rval + get_indent() + "for(int i = 0; i < neighbors.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "     total_size += recv_count[i];\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"


        # Construct additional buffers
        # Looks like this is easy generally, can just use size but might be complex depending on 
        # what the dimensions represent, will check.
        # .extent(dimension) gets the number of elements in each dimension, so if we have views
        # of dimension neighbours.size(), fixed_size, [extra dimension sizes] then we can do 
        # sends of .extent(1) for each neighbour.
        rval = rval + get_indent() + "Kokkos::View<double***, MemorySpace> r_pos_space(\"temp_pos\", neighbors.size(), total_size, 3);\n"
        rval = rval + get_indent() + "Kokkos::View<double***, MemorySpace> r_vel_space(\"temp_vel\", neighbors.size(), total_size, 3);\n"
        rval = rval + get_indent() + "Kokkos::View<double**, MemorySpace> r_cutoff_space(\"temp_cutoff\", neighbors.size(), total_size);\n"

        for element in element_info.keys():
            if element_info[element][1] == 0:
                rval = rval + get_indent() + "Kokkos::View<" + element_info[element][0]
                rval = rval + "**"
                rval = rval + ", MemorySpace> r_" + element + "_space(\"temp_" + element + "\", neighbors.size(), total_size);\n"
            else:
                rval = rval + get_indent() + "Kokkos::View<" + element_info[element][0]
                rval = rval + "**" + element_info[element][1]*"*"
                rval = rval + ", MemorySpace> r_" + element + "_space(\"temp_" + element + "\", neighbors.size(), total_size"
                for i in range(element_info[element][1]):
                    rval = rval + f", {element_info[element][2][i]}"
                rval = rval + ");\n"
        rval = rval + "\n"

        # end of slice creation

        rval = rval + get_indent() + "free(requests);\n"
        val = len(element_info.keys()) + 3
        rval = rval + get_indent() + "requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2 * " + str(val) + ");\n"
        rval = rval + get_indent() + "req_num = 0;\n"
        rval = rval + get_indent() + "int tag = 0;\n"

        # Loop over neighbours and do sends and receives - 2262
        rval = rval + get_indent() + "for(int i = 0; i < neighbors.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "if(neighbors[i] != myrank){\n"
        inc_indent()
        rval = rval + get_indent() + "tag = 0;\n"
        rval = rval + get_indent() + "MPI_Irecv(&r_pos_space.data()[r_pos_space.extent(1)*i], recv_count[i]*3,"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        rval = rval + get_indent() + "MPI_Irecv(&r_vel_space.data()[r_vel_space.extent(1)*i], recv_count[i]*3,"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        rval = rval + get_indent() + "MPI_Irecv(&r_cutoff_space.data()[r_cutoff_space.extent(1)*i], recv_count[i],"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        for element in element_info.keys():
            mpi_dtype = None
            if element_info[element][0] == "int":
                mpi_dtype = "MPI_INT"
            elif element_info[element][0] == "double":
                mpi_dtype = "MPI_DOUBLE"
            elif element_info[element][0] == "float":
                mpi_dtype = "MPI_FLOAT"
            elif element_info[element][0] == "int64_t":
                mpi_dtype = "MPI_INT64_T"
            if mpi_dtype is None:
                raise NotImplementedError("Don't know currently how to support element "
                                          + element + " with datatype "
                                          + element_info[element][0] + " for MPI")
            if element_info[element][1] == 0:
                rval = rval + get_indent() + "MPI_Irecv(&r_"+element+"_space.data()[r_"+element+"_space.extent(1)*i], recv_count[i], "
                rval = rval + mpi_dtype + ", neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
            else:
                rval = rval + get_indent() + "MPI_Irecv(&r_"+element+"_space.data()[r_"+element+"_space.extent(1)"
                for i in range(element_info[element][1]):
                    rval = rval + f" * r_" + element + "_space.extent(" + str(i+2) + ")"
                rval = rval + " * i], recv_count[i]"
                for i in range(element_info[element][1]):
                    rval = rval + f" * {element_info[element][2][i]}"
                rval = rval + ", " + mpi_dtype + ", neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"

        # Now do sends
        rval = rval + get_indent() + "tag = 0;\n"
        rval = rval + get_indent() + "MPI_Isend(&pos_space.data()[pos_space.extent(1)*i], send_count[i]*3,"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        rval = rval + get_indent() + "MPI_Isend(&vel_space.data()[vel_space.extent(1)*i], send_count[i]*3,"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        rval = rval + get_indent() + "MPI_Isend(&cutoff_space.data()[cutoff_space.extent(1)*i], send_count[i],"
        rval = rval + " MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
        for element in element_info.keys():
            mpi_dtype = None
            if element_info[element][0] == "int":
                mpi_dtype = "MPI_INT"
            elif element_info[element][0] == "double":
                mpi_dtype = "MPI_DOUBLE"
            elif element_info[element][0] == "float":
                mpi_dtype = "MPI_FLOAT"
            elif element_info[element][0] == "int64_t":
                mpi_dtype = "MPI_INT64_T"
            if element_info[element][1] == 0:
                rval = rval + get_indent() + "MPI_Isend(&"+element+"_space.data()["+element+"_space.extent(1)*i], send_count[i], "
                rval = rval + mpi_dtype + ", neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"
            else:
                rval = rval + get_indent() + "MPI_Isend(&"+element+"_space.data()["+element+"_space.extent(1)"
                for i in range(element_info[element][1]):
                    rval = rval + f" * r_" + element + "_space.extent(" + str(i+2) + ")"
                rval = rval + " * i], send_count[i]"
                for i in range(element_info[element][1]):
                    rval = rval + f" * {element_info[element][2][i]}"
                rval = rval + ", " + mpi_dtype + ", neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);\n"

        dec_indent()
        rval = rval + get_indent() + "}\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"

        rval = rval + get_indent() + "MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);\n"
        rval = rval + get_indent() + "free(requests);\n"

        # Data is here, need to put it in AOSOA

        rval = rval + get_indent() + "int recvd = 0;\n"
        rval = rval + get_indent() + "int sent = 0;\n"
        rval = rval + get_indent() + "for(int i = 0; i < neighbors.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "recvd += recv_count[i];\n"
        rval = rval + get_indent() + "sent += send_count[i];\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "int size_change = recvd - sent;\n"
        rval = rval + get_indent() + "int current_size =  particle_aosoa.size();\n"
        rval = rval + get_indent() + "if(size_change != 0){\n"
        inc_indent()
        rval = rval + get_indent() + "particle_aosoa.resize(current_size+size_change);\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "auto new_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, \"new_rank\");\n"
        rval = rval + get_indent() + "for(int i = particle_aosoa.size() - 1; i > end; i--){\n"
        inc_indent()
        rval = rval + get_indent() + "new_rank_slice(i) = -1;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "if(size_change > 0){\n"
        inc_indent()
        rval = rval + get_indent() + "if(sent = 0){\n"
        inc_indent()
        rval = rval + get_indent() + "end = current_size;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "while(end < current_size && end < particle_aosoa.size() && new_rank_slice(end) != -1 ) end++;\n"
        rval = rval + get_indent() + "for(int i = 0; i < particle_aosoa.size(); i++){\n"
        inc_indent()
        rval = rval + get_indent() + "new_rank_slice(i) = -1;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "auto new_last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, \"new_last_pos\");\n"
        rval = rval + get_indent() + "auto new_pos_s = Cabana::slice<core_part_position>(particle_aosoa, \"new_position\");\n"
        rval = rval + get_indent() + "auto new_vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, \"new_velocity\");\n"
        rval = rval + get_indent() + "auto new_cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, \"new_cutoff\");\n"
        for element in element_info.keys():
            rval = rval + get_indent() + "auto new_" + element + "_s = Cabana::slice<" + element + ">(particle_aosoa, \"new_" + element + "\");\n"

        rval = rval + get_indent() + "int x = 0;\n"
        rval = rval + get_indent() + "for(int j = 0; j < neighbors.size(); j++){\n"
        inc_indent()
        rval = rval + get_indent() + "for(int i = 0; i < recv_count[j]; i++){\n"
        inc_indent()
        rval = rval + get_indent() + "new_pos_s(end+x, 0) = r_pos_space(j,i,0);\n"
        rval = rval + get_indent() + "new_pos_s(end+x, 1) = r_pos_space(j,i,1);\n"
        rval = rval + get_indent() + "new_pos_s(end+x, 2) = r_pos_space(j,i,2);\n"
        rval = rval + get_indent() + "new_vel_s(end+x, 0) = r_vel_space(j,i,0);\n"
        rval = rval + get_indent() + "new_vel_s(end+x, 1) = r_vel_space(j,i,1);\n"
        rval = rval + get_indent() + "new_vel_s(end+x, 2) = r_vel_space(j,i,2);\n"
        rval = rval + get_indent() + "new_cutoff_s(end+x) = r_cutoff_space(j,i);\n"
        for element in element_info.keys():
            if element_info[element][1] == 0:
                rval = rval + get_indent() + "new_"+element+"_s(end+x) = r_"+element+"_space(j,i);\n"
            else:
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + get_indent() + f"for(int {indexer}=0; {indexer} < {element_info[element][2][index]};"
                    rval = rval + indexer + "++){\n"
                    inc_indent()
                rval = rval + get_indent() + "new_" + element + "_s(end+x"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ") = " + "r_"+element + "_space(j,i"
                indexer = "i"
                for index in range(element_info[element][1]):
                    indexer = indexer + "i"
                    rval = rval + ", " + indexer
                rval = rval + ");\n"
                for index in range(element_info[element][1]):
                    dec_indent()
                    rval = rval + get_indent() + "}\n"
                pass

        rval = rval + get_indent() + "x++;\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        dec_indent()
        rval = rval + get_indent() + "}\n"
        rval = rval + get_indent() + "free(recv_count);\n"
        rval = rval + get_indent() + "free(send_count);\n"
        dec_indent()
        rval = rval + get_indent() + "}};\n\n"


        rval = rval + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + "int get_oneD_rank(int x_r, int y_r, int z_r, int x_ranks, int y_ranks, int z_ranks){\n"
        rval = rval + "    int oneD_rank = z_r*(x_ranks*y_ranks) + y_r*(x_ranks) + x_r;\n"
        rval = rval + "    return oneD_rank;\n"
        rval = rval + "}\n"

        rval = rval + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + "void get_threeD_rank(int rank, int *x, int *y, int *z, int x_ranks, int y_ranks, int z_ranks){\n"
        rval = rval + "    int z_r = rank / (x_ranks*y_ranks);\n"
        rval = rval + "    int y_r = (rank - z_r*x_ranks*y_ranks) / x_ranks;\n"
        rval = rval + "    int x_r = rank - z_r*x_ranks*y_ranks - y_r*x_ranks;\n"
        rval = rval + "    *x = x_r;\n"
        rval = rval + "    *y = y_r;\n"
        rval = rval + "    *z = z_r;\n"
        rval = rval + "};\n"


        rval = rval + "\n"
        rval = rval + "template<class PartPosSlice, class RankSlice>\n"
        rval = rval + "struct _rank_update_functor{\n"
        rval = rval + "    boundary _box;\n"
        rval = rval + "    PartPosSlice _part_pos;\n"
        rval = rval + "    RankSlice _rank;\n"
        rval = rval + "    int _myrank;\n"
        rval = rval + "    int _xranks;\n"
        rval = rval + "    int _yranks;\n"
        rval = rval + "    int _zranks;\n"
        rval = rval + "    int _local_x_rank;\n"
        rval = rval + "    int _local_y_rank;\n"
        rval = rval + "    int _local_z_rank;\n"
        rval = rval + "    int _nranks;\n\n"
        rval = rval + "    KOKKOS_INLINE_FUNCTION\n"
        rval = rval + "    _rank_update_functor(boundary box, PartPosSlice pos, RankSlice rank, int xranks, int yranks, int zranks, int myrank):\n"
        rval = rval + "        _box(box), _part_pos(pos), _rank(rank), _xranks(xranks), _yranks(yranks), _zranks(zranks), _myrank(myrank){\n"
        rval = rval + "        get_threeD_rank(myrank, &_local_x_rank, &_local_y_rank, &_local_z_rank, xranks, yranks, zranks);\n"
        rval = rval + "    }\n\n"
        rval = rval + "    KOKKOS_INLINE_FUNCTION\n"
        rval = rval + "    void operator()(const int ix, const int ij) const{\n"
        rval = rval + "        int xr, yr, zr;\n"
        rval = rval + "        get_threeD_rank(_myrank, &xr, &yr, &zr, _xranks, _yranks, _zranks);\n"
        rval = rval + "        xr = _local_x_rank;\n"
        rval = rval + "        yr = _local_y_rank;\n"
        rval = rval + "        zr = _local_z_rank;\n"
        # Need to compute new ranks BEFORE doing the boundary condition
        rval = rval + "        if(_part_pos.access(ix, ij, 0) >= _box.local_x_max){\n"
        rval = rval + "            xr = xr + 1;\n"
        rval = rval + "            if( xr >= _xranks ) xr = 0;\n"
        rval = rval + "        }\n"
        rval = rval + "        if(_part_pos.access(ix, ij, 0) < _box.local_x_min){\n"
        rval = rval + "            xr = xr - 1;\n"
        rval = rval + "            if( xr < 0 ) xr = _xranks-1;\n"
        rval = rval + "        }\n"
        rval = rval + "        if(_part_pos.access(ix, ij, 1) >= _box.local_y_max){\n"
        rval = rval + "            yr = yr + 1;\n"
        rval = rval + "            if( yr >= _yranks ) yr = 0;\n"
        rval = rval + "        }\n"
        rval = rval + "        if(_part_pos.access(ix, ij, 1) < _box.local_y_min){\n"
        rval = rval + "            yr = yr - 1;\n"
        rval = rval + "            if( yr < 0 ) yr = _yranks-1;\n"
        rval = rval + "        }\n"
        rval = rval + "        if(_part_pos.access(ix, ij, 2) >= _box.local_z_max){\n"
        rval = rval + "            zr = zr + 1;\n"
        rval = rval + "            if( zr >= _zranks ) zr = 0;\n"
        rval = rval + "        }\n"
        rval = rval + "        if(_part_pos.access(ix, ij, 2) < _box.local_z_min){\n"
        rval = rval + "            zr = zr - 1;\n"
        rval = rval + "            if( zr < 0 ) zr = _zranks-1;\n"
        rval = rval + "        }\n"
        # Turn 3D rank index to 1D index.
        rval = rval + "        _rank.access(ix, ij) = get_oneD_rank(xr, yr, zr, _xranks, _yranks, _zranks);\n"
        rval = rval + "    }};\n\n"

        return rval
        # Ideally we just use the Cabana inbuilt function in the future. At that point this 
        # call doesn't need to do as much.

    def gen_headers(self, config, part_type):
        '''
        Generates the headers required for computation of the system.
        Used during code generation to output the header files used
        by the system.

        :param config: The config type used in the system.
        :type config: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Config`
        :param part_type: The part_type used in the system.
        :type part_type: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Particle`
        '''
        self._config = config
        self._part_type = part_type

        # We need to add an extra type to the part type for neighbours_part_deletion_flag
        if "neighbour_part_deletion_flag" not in part_type.particle_type:
            part_type.add_element("neighbour_part_deletion_flag", "int")
        config_output = self.gen_config(config)
        part_output = self.gen_particle(part_type)
        with open('part.h', 'w') as f:
            f.write("#ifndef PART_H\n")
            f.write("#define PART_H\n")
            f.write("#include <Kokkos_Core.hpp>\n")
            f.write("#include <Cabana_Core.hpp>\n")
            extra_includes = set()
            if self._input_module is not None:
                extra = self._input_module.get_header_includes_cabana_pir()
                extra_includes = extra_includes.union(extra)
            if self._output_module is not None:
                extra = self._output_module.get_header_includes_cabana_pir()
                extra_includes = extra_includes.union(extra)
            for coupled_system in self._coupled_systems:
                extra = coupled_system.get_includes_header()
                extra_includes = extra_includes.union(extra)
            for include in extra_includes:
                f.write(f"#include {include}\n")
            if HartreeParticleDSL.get_cuda():
                f.write("using MemorySpace = Kokkos::CudaSpace;\n")
            else:
                f.write("using MemorySpace = Kokkos::HostSpace;\n")
            f.write("using ExecutionSpace = Kokkos::DefaultExecutionSpace;\n")
            f.write("using DeviceType = Kokkos::Device<Kokkos::DefaultExecutionSpace, MemorySpace>;\n")
            f.write("using HostType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;\n")
            f.write("const int VectorLength = 16;\n\n") #This vector length should be chosen better in future

            for name in self._global_values:
                symbol = HartreeParticleDSL.global_symbol_table().lookup(name)
                dt_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)
                if HartreeParticleDSL.get_cuda():
                    if self._global_values[name] is None:
                        raise NotImplementedError("Non-constant global variables are"
                                " not supported in Cabana_PIR with CUDA enabled.")
                    f.write(f"#define {name} {self._global_values[name]}\n")
                else:
                    if self._global_values[name] is not None:
                        f.write(f"extern {dt_str} {name};\n")
                    else:
                        f.write(f"extern {dt_str} {name};\n")

            f.write(config_output)
            f.write(part_output)

            # If we are generating for MPI then create the MPI header code
            if HartreeParticleDSL.get_mpi():
                f.write(self.mpi_headers(config, part_type))

        with open('part.h', "a") as f:
            f.write("#endif")

    def gen_kernel(self, kernel):
        '''
        Generates the Cabana code for a given user-defined kernel.
        Used during code generation to output the kernel code in C++.

        :param kernel: Kernel object to generate code for.
        :type kernel: :py:class:`HartreeParticleDSL.kernel_types.kernels.kernel`
        '''
        # The input to this is the Python AST for this kernel.
        # First we need to convert the Python AST to Particle IR.
        # Store the resulting Particle IR node in this for later codegen.
        # In principle we could do this codegen now as well, but I sort of want
        # to wait and do all the codegen from print_main, and output to file
        # for this backend.
        from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import \
                pir_perpart_visitor, pir_pairwise_visitor, pir_source_boundary_visitor, \
                pir_sink_boundary_visitor
        pir = None
        if isinstance(kernel, kernels.perpart_kernel_wrapper):
            conv = pir_perpart_visitor()
            pir = conv.visit(kernel.get_kernel_tree())
        elif isinstance(kernel, kernels.pairwise_kernel_wrapper):
            conv = pir_pairwise_visitor()
            pir = conv.visit(kernel.get_kernel_tree())
        elif isinstance(kernel, kernels.source_boundary_kernel_wrapper):
            conv = pir_source_boundary_visitor(kernel)
            pir = conv.visit(kernel.get_kernel_tree())
        elif isinstance(kernel, kernels.sink_boundary_kernel_wrapper):
            conv = pir_sink_boundary_visitor()
            pir = conv.visit(kernel.get_kernel_tree())
        else:
            raise NotImplementedError("Cabana PIR backend doesn't yet support "
                                      f"kernel of type {type(kernel)}.")

        name = pir.name
        self._kernels_pir[name] = pir

    def get_current_kernel(self):
        '''
        Returns the current PIR kernel that has code generation in progress.
        This is required so coupled systems can create particle accesses correctly.

        :returns: The current PIR kernel
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kern.Kern`
        '''
        return self._current_kernel


    def print_main(self, function):
        '''
        Generates the Cabana code for the supplied main function.

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import \
                ast_to_pir_visitor, pir_main_visitor
        from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
        from HartreeParticleDSL.Particle_IR.nodes.kernels import MainKernel
        cabana_pir = Cabana_PIR_Visitor(self)
        # Plan for this:
        # Find where the output files should go
        # Generate the header file into that directory
        self.gen_headers(self._config, self._part_type)

        # Copy related files from the coupled systems
        for sys in self._coupled_systems:
            sys.copy_files()

        # TODO Generate makefile
        # self._gen_makefile()

        # Open an output C++ file
        # Find all the kernels used by this main function and 
        # generate those kernels.
        conv = pir_main_visitor() # ast_to_pir_visitor()
        pir = conv.visit(function)
#        if not isinstance(pir, MainKernel):
#            body = pir.body.detach()
#            sym_tabl = pir.symbol_table
#            pir = MainKernel("main")
#            pir.children[0] = body
#            pir._symbol_table = sym_tabl

        kernels = pir.body.walk(Invoke)
        names = []
        for kernel in kernels:
            for arg in kernel.children:
                names.append(arg.value)

        codes = []
        if self._boundary_condition is not None:
            self.gen_kernel(self._boundary_condition)
            names.append(self._boundary_condition_tree.body[0].name)

        for name in names:
            kern_pir = self._kernels_pir[name]
            self._current_kernel = kern_pir
            codes.append(cabana_pir(kern_pir))
        self._current_kernel = None

        main_code = cabana_pir(pir)

        includes = self.generate_includes()

        with open('code.cpp', 'w') as f:
            f.write(includes + "\n")
            for name in self._global_values:
                symbol = HartreeParticleDSL.global_symbol_table().lookup(name)

# This code is covered in gen_headers so unreachable here.
#                if HartreeParticleDSL.get_cuda():
#                    if self._global_values[name] is None:
#                        raise NotImplementedError("Non-constant global variables are"
#                                " not supported in Cabana_PIR with CUDA enabled.")
                if not HartreeParticleDSL.get_cuda():
                    dt_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)
                    if self._global_values[name] is not None:
                        f.write(f"{dt_str} {name} = {self._global_values[name]};\n")
                    else:
                        f.write(f"{dt_str} {name};\n")
            input_module_header = ""
            if self._input_module is not None:
                input_module_header = self._input_module.gen_code_cabana_pir(self._part_type) #FIXME
            if input_module_header != "":
                f.write(input_module_header)
                f.write("\n")

            output_module_header = ""
            if self._output_module is not None and self._input_module is not self._output_module:
                output_module_header = self._output_module.gen_code_cabana_pir(self._part_type) #FIXME
            if output_module_header != "":
                # Do something later
                f.write(output_module_header)
                f.write("\n")
            for code in codes:
                f.write(code + "\n")
            f.write(main_code)


    def set_cutoff(self, cutoff, var_type=CONSTANT):
        '''
        Set the cutoff for pairwise interactions. NYI

        :raises: NotImplementedError
        '''
        raise NotImplementedError("Cabana PIR backend doesn't yet support pairwise interactions")

    def initialisation_code(self, particle_count, filename):
        '''
        Get the initialisation code for this code. This is done by the input IO module.

        :param int particle_count: The number of particles to create for this testcase. \
                Note that depending on the IO module this value may be ignored, e.g. \
                for IO modules that read the input from a file this is likely ignored.
        :param str filename: The filename to read the input from. Note that depending \
                on the IO module this value may be ignored, e.g. for the random particle \
                IO module this is ignored.
        '''
        return self._input_module.call_input_cabana_pir(particle_count, filename)

    def gen_particle(self, particle) -> str:
        '''
        Generate the Cabana-compatible code for the particle structure provided.

        Will generate extra sections if required for MPI.

        :param particle: The particle to generate the Cabana code for.
        :type particle: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Particle`

        :returns: The code string for the particle.
        :rtype: str
        '''
        # Store the particle for later
        self._particle = particle
        output = ""


        # Sort particle by size of element
        sizes = []
        for element in particle.particle_type:
            if element == "core_part" or element == "neighbour_part":
                sizes.append(0)
                continue
            c_type = particle.particle_type[element]['type']
            is_array = particle.particle_type[element]['is_array']
            if is_array:
                x = c_type.index("[")
                base_type = c_type[0:x]
                count = 1
                array_str = c_type[x:]
                while array_str.find("[") >= 0:
                    val = int(array_str[array_str.index("[")+1:array_str.index("]")])
                    count = count * val
                    array_str = array_str[array_str.index("]")+1:]
                size = Cabana_PIR._type_sizes[base_type]
                sizes.append(size * count)
            else:
                size = Cabana_PIR._type_sizes[c_type]
                sizes.append(size)
        fields = particle.particle_type.keys()
        sorted_fields = [x for _, x in sorted(zip(sizes, fields), reverse=True)]

        output = output + "enum FieldNames{"
        first = True
        for key in sorted_fields:
            if key != "core_part" and key != "neighbour_part" and "neighbour_part" not in key:
                if first:
                    output = output + f"{key} = 0"
                    first = False
                else:
                    output = output + f",\n                 {key}"
        # Extra field names for core part and neighbour part things
        if first:
            output = output + "core_part_velocity = 0"
        else:
            output = output + ",\n                 core_part_velocity"
        output = output + ",\n                 core_part_position"
        output = output + ",\n                 neighbour_part_cutoff"
        # Delete particles info
        output = output + ",\n                 neighbour_part_deletion_flag"
        if HartreeParticleDSL.get_mpi():
            output = output + ",\n                 neighbour_part_rank"
            output = output + ",\n                 neighbour_part_old_position"
        output = output + "\n               };\n"

        output = output + "using DataTypes = Cabana::MemberTypes<"
        first = True
        for key in sorted_fields:
            if key != "core_part" and key != "neighbour_part" and "neighbour_part" not in key:
                if first:
                    c_type = particle.particle_type[key]['type']
                    output = output + c_type
                    first = False
                else:
                    c_type = particle.particle_type[key]['type']
                    output = output + ",\n    " + c_type
        # Extra types for core part and neighbour part things
        if first:
            output = output + "double[3]"
        else:
            output = output + ",\n    double[3]"
        output = output + ",\n    double[3]"
        output = output + ",\n    double"
        output = output + ",\n    int"
        if HartreeParticleDSL.get_mpi():
            output = output + ",\n    int"
            output = output + ",\n    double[3]"
        output = output + ">;\n"

        return output

    def gen_config(self, config) -> str:
        '''
        Generate the Cabana-compatible code for the config structure provided.

        Will generate extra sections if required for MPI.

        :param config: The particle to generate the Cabana code for.
        :type config: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Config`

        :returns: The code string for the config.
        :rtype: str
        '''
        output = ""
        output = output + "struct boundary{\n"
        output = output + "    double x_min, x_max;\n"
        output = output + "    double y_min, y_max;\n"
        output = output + "    double z_min, z_max;\n"
        if HartreeParticleDSL.get_mpi():
            output = output + "    double local_x_min, local_x_max;\n"
            output = output + "    double local_y_min, local_y_max;\n"
            output = output + "    double local_z_min, local_z_max;\n"
            output = output + "    int x_ranks, y_ranks, z_ranks;\n"
        output = output + "};\n\n"

        output = output + "struct space_type{\n"
        output = output + "    boundary box_dims;\n"
        output = output + "    int nparts;\n"
        output = output + "};\n\n"

        # Currently empty neiughbour config
        output = output + "struct neighbour_config_type{\n"
        output = output + "};\n\n"

        output = output + "struct config_view_type{\n"
#        output = output + "    space_type space;\n"
#        output = output + "    neighbour_config_type neighbour_config;\n"
        for key in config.config_type:
            is_array = config.config_type[key]['is_array']
            c_type = config.config_type[key]['type']
            varnam = key
            if is_array:
                # Move array index to the correct position for C++ definitions
                x = c_type.index("[")
                varnam = varnam + c_type[x:]
                c_type = c_type[0:x]
            output = output + f"        {c_type} {varnam};\n"
        output = output + "};\n\n"

        output = output + "using config_struct_type = Kokkos::View<struct config_view_type*, MemorySpace>;\n"
        output = output + "using config_struct_host = config_struct_type::HostMirror;\n"

        output = output + "struct config_type{\n"
        output = output + "    config_struct_type config;\n"
        output = output + "    config_struct_host config_host;\n"
        output = output + "};\n"

        # There are some complex things to work out here. How do we add other Kokkos views into the config_type struct
        # We probably need to do some extra type analysis in the config_type object.
        return output

    def cleanup(self, *args, current_indent, **kwargs) -> str:
        '''
        :param int current_indent: Current indentation level for this code. \
                                   Unused.

        :returns: The cleanup code for this backend. This doesn't do much as \
                  Kokkos' initialisation is done through scoping.
        :rtype: str
        '''
        rval = "\n}\n" # Choosing no indentation for this } to match the outer one.
        return rval

    def println(self, string, *args, **kwargs):
        '''
        Function to output the println string for the Cabana PIR module.
        Called via the visitors when reaching a println statement in user
        code.

        For the Cabana module this is a call to `std::cout`, using the
        string argument as the first value to be passed into cout, and
        the args contains any further values used in the output.

        :param string: The formatted string to use with cout.                                                                                                                      :type string: str
        :param int current_indent: The current indentation level
        :param args: A list of strings containing the other values
                      to output with cout. Any strings to add to cout
                      should be surrounded with \"\".
        :type args: str
        '''
        current_indent = kwargs.get("current_indent", 0)
        output = " "*current_indent + f"std::cout << {string}"
        for arg in args:
            x = arg[0].replace('"', '') + arg[1:]
            x = x[0:len(x)-1] + x[-1].replace('"', '')

            output = output + f" << {x}"
        output = output + " << \"\\n\";\n"
        return output

    def random_number(self, *args, **kwargs):
        '''
        TODO
        '''
        return "_generator.drand(0., 1.)"

    def pi(self, *args, **kwargs):
        '''
        TODO
        '''
        return "M_PI"

    def initialise(self,particle_count, filename, current_indent, **kwargs) -> str:
        '''
        :param int particle_count: Particle count for the simulation. Passed to \
                the IO module where it may be ignored. Check the IO module for info.
        :param str filename: The filename to read input data from. Passed to the \
                IO m odule where it may be ignored. Check the IO module for info.
        :param int current_indent: The current indentation level of code to start \
                generating at.

        :returns: the code to initialise Kokkos for this backend, and the \
                  MPI initialisation if required.
        :rtype: str
        '''
        space = " "
        rval = space*current_indent + "Kokkos::ScopeGuard scope_guard(argc, argv);\n"
        if HartreeParticleDSL.get_mpi():
            rval = rval + space*current_indent + "int _provided;\n"
            rval = rval + space*current_indent + "MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &_provided  );\n"

        rval = rval + "{\n"
        rval = rval + space*current_indent + "config_type config;\n"
        rval = rval + space*current_indent + "config.config = config_struct_type(\"config\", 1);\n"
        rval = rval + space*current_indent + "config.config_host = Kokkos::create_mirror_view(config.config);\n"
        if self.get_require_random():
            rval = rval + space*current_indent + "Kokkos::Random_XorShift64_Pool<> random_pool(/*seed=*/12345);\n"

        # Initialise structures
        for structure in self._structures:
            if self._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._structures[structure]][0]
            else:
                typename = structure
            rval = rval + space*current_indent + f"{typename} {structure};\n"

        # Load the base box dimensions
        rval = rval + self._input_module.call_get_box_size_pir(int(particle_count), filename, current_indent=current_indent)

        # If we have MPI then we load the rank info
        if HartreeParticleDSL.get_mpi():
            rval = rval + space*current_indent + "int myrank = 0; MPI_Comm_rank( MPI_COMM_WORLD, &myrank );\n"
            rval = rval + space*current_indent + "int nranks = 1; MPI_Comm_size( MPI_COMM_WORLD, &nranks );\n"
        else:
            rval = rval + space*current_indent + "int myrank = 0;\n"
            rval = rval + space*current_indent + "int nranks = 1;\n"

        # Create any extra structures that are needed
        for struct in self._structures.keys():
            val = rval + space*current_indent + "struct " + struct + " " + struct + ";\n"

        # If we have MPI then we need a domain decomposition.
        # First we check if a coupled system has a preference
        if HartreeParticleDSL.get_mpi():
            decomposition_done = False
            for coupler in self._coupled_systems:
                if coupler.has_preferred_decomposition():
                    # If the coupler has a preference on domain decomposition then
                    # we intialise the coupled system and then take its domain
                    # decomposition
                    if decomposition_done:
                        raise NotImplementedError("Can't handle multiple coupled systems with a preferred decomposition")
                    # Initialise the coupler
                    rval = rval + coupler.setup_testcase(filename, current_indent=current_indent)
                    rval = rval + coupler.get_preferred_decomposition("config.config_host(0).space.box_dims", current_indent=current_indent)
                    decomposition_done = True

            # If the decomposition wasn't done by the coupler we need to do it here
            if not decomposition_done:
                raise NotImplementedError("Can't yet do a decomposition when coupled systems didn't do it for us.")
        else:
            # If no MPI, just setup testcase
            for coupler in self._coupled_systems:
                if coupler.has_preferred_decomposition():
                    rval = rval + coupler.setup_testcase(filename, current_indent=current_indent)
            

        # Load the input file.
        rval = rval + space*current_indent + f"{self._input_module.call_input_cabana_pir(int(particle_count), filename, current_indent=current_indent)}\n"
        rval = rval + space*current_indent + "Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());\n"

        # For each other coupled system we initialise them now.
        for coupled in self._coupled_systems:
            if not coupled.has_preferred_decomposition():
                rval = rval + coupled.setup_testcase(filename, current_indent=current_indent)

        # Initialise neighbours list
        if HartreeParticleDSL.get_mpi():
            rval = rval + space*current_indent + "std::vector<int> neighbors = {myrank};\n"
            rval = rval + space*current_indent + "int local_x, local_y, local_z;\n"
            rval = rval + space*current_indent + "get_threeD_rank(myrank, &local_x, &local_y, &local_z, config.config_host(0).space.box_dims.x_ranks,\n"
            rval = rval + space*current_indent + "     config.config_host(0).space.box_dims.y_ranks,  config.config_host(0).space.box_dims.z_ranks);\n"
            rval = rval + space*current_indent + "for(int x = local_x - 1; x <= local_x + 1; x++){\n"
            rval = rval + space*current_indent + "    for(int y = local_y - 1; y <= local_y + 1; y++){\n"
            rval = rval + space*current_indent + "        for( int z = local_z - 1; z <= local_z + 1; z++){\n"
            rval = rval + space*current_indent + "            int xr = x; int yr = y; int zr = z;\n"
            rval = rval + space*current_indent + "            if(xr < 0) xr = config.config_host(0).space.box_dims.x_ranks-1;\n"
            rval = rval + space*current_indent + "            if(xr >=  config.config_host(0).space.box_dims.x_ranks) xr = 0;\n"
            rval = rval + space*current_indent + "            if(yr < 0) yr =  config.config_host(0).space.box_dims.y_ranks-1;\n"
            rval = rval + space*current_indent + "            if(yr >=  config.config_host(0).space.box_dims.y_ranks) yr = 0;\n"
            rval = rval + space*current_indent + "            if(zr < 0) zr =  config.config_host(0).space.box_dims.z_ranks-1;\n"
            rval = rval + space*current_indent + "            if(zr >=  config.config_host(0).space.box_dims.z_ranks) zr = 0;\n"
            rval = rval + space*current_indent + "            neighbors.push_back(get_oneD_rank(xr, yr, zr,  config.config_host(0).space.box_dims.x_ranks,\n"
            rval = rval + space*current_indent + "                 config.config_host(0).space.box_dims.y_ranks,\n"
            rval = rval + space*current_indent + "                 config.config_host(0).space.box_dims.z_ranks));\n"
            rval = rval + space*current_indent + "        }\n"
            rval = rval + space*current_indent + "    }\n"
            rval = rval + space*current_indent + "}\n"
            rval = rval + space*current_indent + "std::sort( neighbors.begin(), neighbors.end() );\n"
            rval = rval + space*current_indent + "auto unique_end = std::unique( neighbors.begin(), neighbors.end() );\n"
            rval = rval + space*current_indent + "neighbors.resize( std::distance( neighbors.begin(), unique_end ) );\n"

            # Initialise the Migrator
            rval = rval + space*current_indent + "Migrator<decltype(particle_aosoa_host)> _migrator(particle_aosoa.size() / 10, neighbors.size());\n"

        # Need to do something with each kernel now.

        # We need the particle type to be able to initialise correctly
        # The initialisation of the core part and neighbour cutoff stuff is a bit weird, needs to be fixed.
        rval = rval + space*current_indent + "auto core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);\n"
        rval = rval + space*current_indent + "auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);\n"
        rval = rval + space*current_indent + "auto neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);\n"
        rval = rval + space*current_indent + "auto neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);\n"
        if HartreeParticleDSL.get_mpi():
            rval = rval + "    auto neighbour_part_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa);\n"
            rval = rval + "    auto neighbour_part_old_position_slice = Cabana::slice<neighbour_part_old_position>(particle_aosoa);\n"
        for key in self._particle.particle_type:
            if key != "core_part" and key != "neighbour_part" and "neighbour_part" not in key:
                rval = rval + space*current_indent + f"auto {key}_slice = Cabana::slice<{key}>" + "(particle_aosoa);\n"


        # Generate the functors
        for key in self._kernel_slices.keys():
            rval = rval + space*current_indent + f"{key}_functor"
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"decltype({slices}_slice)")
            if len(slice_names) > 0:
                rval = rval + "<"
                rval = rval + ", ".join(slice_names) + ">"

            rval = rval + f" {key}("
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"{slices}_slice")
            rval = rval + ", ".join(slice_names) + ", config.config"
            slice_names = []
            for struct in self._structures.keys():
                slice_names.append(struct)
            if len(slice_names) > 0:
                rval = rval + ", "
            rval = rval + ", ".join(slice_names) + ");\n"

        # TODO Initial MPI Communication and binning of particles
        # See issue #62
        return rval

    def add_coupler(self, coupled_system):
        '''
        Adds a coupled system to the list of coupled systems in the backend

        :param coupled_system: The object to couple with.
        :type coupled_system: Object

        :raises UnsupportedTypeError: If the object to couple with is not \
                                      an instance of base_coupler
        '''
        if not isinstance(coupled_system, base_coupler):
            raise UnsupportedTypeError("Can only couple to base_coupler classes "
                                       "or subclasses. Found {0}".format(
                                           type(coupled_system).__name__))
        self._coupled_systems.append(coupled_system)
        extra_includes = coupled_system.get_includes()
        for include in extra_includes:
            self._includes.append(include)

    def write_output(self, filename: str, variable:str=None, **kwargs) -> str:
        '''
        Generates the code to write a file output using the selected output module.

        :param str filename: The filename to write the file to.

        '''
        current_indent = kwargs.get("current_indent", 0)
        code = " " * current_indent
        code = code + self._output_module.call_output_cabana_pir(0, filename, variable, current_indent=current_indent) + "\n"
        return code

    def call_language_function(self, func_call: str, *args, **kwargs) -> str:
        '''
        Calls the backend/coupled system inbuilt function with func_call name
        using the associated arguments.

        :param str func_call: The name of the inbuilt function to call, e.g. \
                println

        :raises AttributeError: If the function is not found in the backend or \
                any coupled system.

        :returns: The C++ code associated with the function call for the given \
                arguments.
        :rtype: str
        '''
        string = ""
        try:
            fn = getattr(self, func_call)
            string = fn( *args, **kwargs)
            return string
        except (AttributeError) as err:
            pass
        for system in self._coupled_systems:
            try:
                fn = getattr(system, func_call)
                string = fn(*args, **kwargs)
                return string
            except (AttributeError) as err:
                pass
        raise AttributeError

    def get_require_random(self) -> bool:
        return self._require_random

    def get_extra_symbols(self, function_list):
        '''
        Returns the list of symbols that need to be added to the symbol table
        to correctly resolve any function calls from this node.

        :param function_list: List of names of calls used in this function.
        :type function_list: List of str.

        :returns: The list of tuples (name, Symbol) required for this function.
        '''
        # This backend needs to add any structures created.
        results = []
        if "random_number" in function_list:
            self._require_random = True
            results.append( ("_generator", AutoSymbol("_generator", "_random_pool.get_state(); _random_pool.free_state(_generator);"))) # TODO Fix.
        for struct in self._structures.keys():
            results.append( (struct, StructureSymbol(struct, self._structures[struct])) )


        # Check the coupled systems as well.
        for system in self._coupled_systems:
            results.extend(system.get_extra_symbols(function_list))
        return results

    def gen_mpi_comm_before_bcs(self, current_indent=4, indent=4) -> str:
        rval = ""
        indent_str = current_indent * " " 
        rval = rval + indent_str + "{\n"
        indent_str = indent_str + " " * indent
        rval = rval + indent_str + "_rank_update_functor<decltype(core_part_position_slice), decltype(neighbour_part_rank_slice)> ruf(\n"
        rval = rval + indent_str + "    config.config_host(0).space.box_dims, core_part_position_slice, neighbour_part_rank_slice,\n"
        rval = rval + indent_str + "    config.config_host(0).space.box_dims.x_ranks,\n"
        rval = rval + indent_str + "    config.config_host(0).space.box_dims.y_ranks,\n"
        rval = rval + indent_str + "    config.config_host(0).space.box_dims.z_ranks, myrank);\n"
        rval = rval + indent_str + "Cabana::SimdPolicy<VectorLength, ExecutionSpace> temp_policy(0, particle_aosoa.size());\n"
        rval = rval + indent_str + "Cabana::simd_parallel_for( temp_policy, ruf, \"rank_update_functor\");\n"
        rval = rval + indent_str + "Kokkos::fence();\n"
        rval = rval + current_indent * " " + "}\n"
        return rval

    def gen_slices_and_functors(self, current_indent=4, indent=4) -> str:
        '''
        Returns the code required to regenerate slices and functors for
        this backend. Used if something resizes the particle arrays.

        :returns: The code to regenerate the slices and functors
        :rtype: str
        '''
        indent_str = current_indent * " "
        rval = ""
        rval = rval + indent_str + "core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);\n"
        rval = rval + indent_str + "core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);\n"
        rval = rval + indent_str + "neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);\n"
        rval = rval + indent_str + "neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);\n"
        if HartreeParticleDSL.get_mpi():
            rval = rval + indent_str + "neighbour_part_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa);\n"
            rval = rval + indent_str + "neighbour_part_old_position_slice = Cabana::slice<neighbour_part_old_position>(particle_aosoa);\n"
        if self._particle is not None:
            for key in self._particle.particle_type:
                if key != "core_part" and key != "neighbour_part":
                    rval = rval + indent_str + f"{key}_slice = Cabana::slice<{key}>" + "(particle_aosoa);\n"

        rval = rval + indent_str + f"simd_policy = Cabana::SimdPolicy<VectorLength, ExecutionSpace>(0, particle_aosoa.size());\n"

        # Generate the functors
        for key in self._kernel_slices.keys():
            rval = rval + indent_str + f"{key} = {key}_functor"
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"decltype({slices}_slice)")
            if len(slice_names) > 0:
                rval = rval + "<"
                rval = rval + ", ".join(slice_names) + ">"

            rval = rval + f"("
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"{slices}_slice")
            rval = rval + ", ".join(slice_names) + ", config.config"
            slice_names = []
            for struct in self._structures.keys():
                slice_names.append(struct)
            if len(slice_names) > 0:
                rval = rval + ","
            rval = rval + ", ".join(slice_names) + ");\n" 
        return rval

    def gen_mpi_comm_after_bcs(self, current_indent=4, indent=4) -> str:
        '''
        Returns the code required to do the MPI communication for this
        backend.

        Two steps required for this process.
        First is to update the rank that each particle should be on. This is
        done in parallel with a Cabana parallel loop.
        Second is to call the migrator - currently this is code generated by my
        code, but ideally we can use Cabana's inbuilt migrate code.

        Once the binning code is also implemented we need to do this code here.

        :returns: The code to do MPI communication for this backend.
        :rtype: str
        '''
        indent_str = current_indent * " "
        rval = ""
        rval = rval + indent_str + "Cabana::deep_copy(particle_aosoa_host, particle_aosoa);\n"
        rval = rval + indent_str + "_migrator.exchange_data(particle_aosoa_host, neighbors, myrank, particle_aosoa_host.size());\n"
        rval = rval + indent_str + "particle_aosoa.resize(particle_aosoa_host.size());\n"
        rval = rval + indent_str + "Cabana::deep_copy(particle_aosoa, particle_aosoa_host);\n"

        # Need to remake slices and SimdPolicy as well.
        rval = rval + indent_str + "core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);\n"
        rval = rval + indent_str + "core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);\n"
        rval = rval + indent_str + "neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);\n"
        rval = rval + indent_str + "neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);\n"
        if HartreeParticleDSL.get_mpi():
            rval = rval + indent_str + "neighbour_part_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa);\n"
            rval = rval + indent_str + "neighbour_part_old_position_slice = Cabana::slice<neighbour_part_old_position>(particle_aosoa);\n"
        if self._particle is not None:
            for key in self._particle.particle_type:
                if key != "core_part" and key != "neighbour_part":
                    rval = rval + indent_str + f"{key}_slice = Cabana::slice<{key}>" + "(particle_aosoa);\n\n"

        rval = rval + indent_str + f"simd_policy = Cabana::SimdPolicy<VectorLength, ExecutionSpace>(0, particle_aosoa.size());\n\n"


        # Generate the functors
        for key in self._kernel_slices.keys():
            rval = rval + indent_str + f"{key} = {key}_functor"
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"decltype({slices}_slice)")
            if len(slice_names) > 0:
                rval = rval + "<"
                rval = rval + ", ".join(slice_names) + ">"

            rval = rval + f"("
            slice_names = []
            for slices in self._kernel_slices[key]:
                slice_names.append(f"{slices}_slice")
            rval = rval + ", ".join(slice_names) + ", config.config"
            slice_names = []
            for struct in self._structures.keys():
                slice_names.append(struct)
            if len(slice_names) > 0:
                rval = rval + ","
            rval = rval + ", ".join(slice_names) + ");\n" 
        return rval

    def _gen_makefile(self):
        raise NotImplementedError()
#        with open("CMakeLists.txt", "w") as f:
#            #TODO set project
#            import os
#            project_name = os.path.split(os.getcwd())[-1]
#            #TODO
#            output = "cmake_minimum_required (VERSION 3.10)\n"
#            output = output + f"project ({project_name})\n"
#            output = output + "set(Kokkos_DIR \"$ENV{Kokkos_ROOT}\" CACHE STRING \"Kokkos root directory\")\n"
#            required_packages = []
#            required_packages.append("Kokkos")
#            required_packages.append("Cabana")
#            if get_mpi():
#                required_packages.append("MPI")
#            # TODO IF we have GPU we need CUDAToolkit
#            # Check with IO and coupled systems for their needs
#            for sys in self._coupled_systems:
#                sys_req_pack = sys.get_required_packages()
#                for req in sys_req_pack:
#                    if req not in required_packages:
#                        required_packages.append(req)
#                
#            for package in required_packages:
#                output = output + f"find_package({package}, REQUIRED)\n"
#
#            if get_mpi():
#                output = output + "include_directories(SYSTEM ${MPI_INCLUDE_PATH})\n"
#
#            compiled_files = ["code.cpp"]
#            for sys in self._coupled_systems:
#                sys_files = sys.compilation_files()
#                for fil in sys_files:
#                    if fil not in compiled_files:
#                        compiled_files.append(fil)
#            
#            output = output + "add_executable(" + project_name + ".exe "
#            output = output + " ".join(compiled_files) + ")\n"
#
##target_link_libraries(Cabana_PIC.exe ${HDF5_C_LIBRARIES})
#            link_libraries = ["Cabana::cabanacore", "Kokkos::kokkos"]
#            if get_mpi():
#                link_libraries.append("${MPI_C_LIBRARIES}")
#            # TODO If we have GPU we need CUDArt or something
#            # TODO Implement coupled system and IO requirements
#            for lib in link_libraries:
#                output = output + "target_link_libraries(" + project_name + ".exe" + lib + ")\n"
#
#            f.write(output)
