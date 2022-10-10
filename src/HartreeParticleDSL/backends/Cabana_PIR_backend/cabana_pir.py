from __future__ import annotations

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin
from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import Cabana_PIR_Visitor
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, DataType


class Cabana_PIR(Backend):
    '''
    Cabana Backend class using PIR.
    Outputs from this class use the Cabana & Kokkos system to compute
    parallel particle simulations.
    '''
    _type_map = {c_int : "int",
            c_double : "double",
            c_float : "float",
            c_int64_t : "int64_t",
            c_int32_t : "int32_t",
            c_int8_t : "int8_t",
            c_bool : "bool",
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

    @property
    def structures(self):
        return self._structures

    def register_kernel(self, kernel_name: str, kernel: Kern):
        if kernel_name in self._kernels.keys():
            raise InvalidNameError(f"Kernel with name {kernel_name} already exists.") 
        self._kernels[kernel_name] = kernel


    def add_structure(self, structure_type: StructureType, structure_name: str):
        # Should this check the structure type is defined?
        self._structures[structure_name] = structure_type

    def add_kernel_slices(self, kernel_name, kernel_slice):
        '''
        Adds access to a particle slice for a specific kernel.
        '''
        if kernel_name not in self._kernel_slices.keys():
            self._kernel_slices[kernel_name] = []
        if kernel_slice not in self._kernel_slices[kernel_name]:
            self._kernel_slices[kernel_name].append(kernel_slice)

    def add_type(self, type_name: str, the_type: DataType) -> None:
        '''
        Adds a new type to Cabana & PIR.

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
        #TODO
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
        if input_module is not None and not isinstance(input_module, Cabana_IO_Mixin):
            raise InvalidIOModuleError(f"{self.__class__.__name__} backend does not support "
                                       f"{input_module.__class__.__name__} IO "
                                       f"module at this time.")
        if output_module is not None and not isinstance(output_module, Cabana_IO_Mixin):
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
            input_includes = self._input_module.get_includes_cabana() #FIXME
            for include in input_includes:
                if include not in self._includes:
                    self._includes.append(include)
        if self._output_module is not None:
            output_includes = self._output_module.get_includes_cabana() #FIXME
            for include in output_includes:
                if include not in self._includes:
                    self._includes.append(include)

        include_string = ""
        for include in self._includes:
            include_string = include_string + f"#include {include}\n"
        return include_string

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
        config_output = self.gen_config(config)
        part_output = self.gen_particle(part_type)
        with open('part.h', 'w') as f:
            f.write("#ifndef PART_H\n")
            f.write("#define PART_H\n")
            f.write("#include <Kokkos_Core.hpp>\n")
            f.write("#include <Cabana_Core.hpp>\n")
            for coupled_system in self._coupled_systems:
                extra_includes = coupled_system.get_includes_header()
                for include in extra_includes:
                    f.write(f"#include {include}\n")
            f.write("/*using MemorySpace = Kokkos::CudaSpace;*/\n")
            f.write("using MemorySpace = Kokkos::HostSpace;\n")
            f.write("using ExecutionSpace = Kokkos::DefaultExecutionSpace;\n")
            f.write("using DeviceType = Kokkos::Device<Kokkos::DefaultExecutionSpace, MemorySpace>;\n")
            f.write("using HostType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;\n")
            f.write("const int VectorLength = 16;\n\n") #This vector length should be chosen better in future

            for name in self._global_values:
#                print(HartreeParticleDSL.global_symbol_table().lookup(name))
                symbol = HartreeParticleDSL.global_symbol_table().lookup(name)
                dt_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)
                if self._global_values[name] is not None:
                    f.write(f"{dt_str} {name} = {self._global_values[name]};\n")
                else:
                    f.write(f"{dt_str} {name};\n")

            f.write(config_output)
            f.write(part_output)
            f.write("#endif")
        input_module_header = ""
        if self._input_module is not None:
            input_module_header = self._input_module.gen_code_cabana(part_type) #FIXME
        if input_module_header != "":
            print(input_module_header)
            print("\n")

        output_module_header = ""
        if self._output_module is not None:
            output_module_header = self._output_module.gen_code_cabana(part_type) #FIXME
        if output_module_header != "":
            # Do something later
            print(output_module_header)
            print("\n")

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
                pir_perpart_visitor, pir_pairwise_visitor
        pir = None
        if isinstance(kernel, kernels.perpart_kernel_wrapper):
            conv = pir_perpart_visitor()
            pir = conv.visit(kernel.get_kernel_tree())
        elif isinstance(kernel, kernels.pairwise_kernel_wrapper):
            conv = pir_pairwise_visitor()
            pir = conv.visit(kernel.get_kernel_tree())
        else:
            raise NotImplementedError("Cabana PIR backend doesn't yet support "
                                      f"kernel of type {type(kernel)}.")

        name = pir.name
        self._kernels_pir[name] = pir

    def print_main(self, function):
        '''
        Generates the Cabana code for the supplied main function.

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import \
                ast_to_pir_visitor
        from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
        from HartreeParticleDSL.Particle_IR.nodes.kernels import MainKernel
        cabana_pir = Cabana_PIR_Visitor(self)
        # Plan for this:
        # Find where the output files should go
        # Generate the header file into that directory
        self.gen_headers(self._config, self._part_type)
        # Open an output C++ file
        # Find all the kernels used by this main function and 
        # generate those kernels.
        conv = ast_to_pir_visitor()
        pir = conv.visit(function)
        if not isinstance(pir, MainKernel):
            body = pir.body.detach()
            sym_tabl = pir.symbol_table
            pir = MainKernel("main")
            pir.children[0] = body
            pir._symbol_table = sym_tabl

        kernels = pir.body.walk(Invoke)
        names = []
        for kernel in kernels:
            for arg in kernel.children:
                names.append(arg.value)

        codes = []
        for name in names:
            kern_pir = self._kernels_pir[name]
            codes.append(cabana_pir(kern_pir))

        main_code = cabana_pir(pir)

        includes = self.generate_includes()

        with open('code.cpp', 'w') as f:
            f.write(includes + "\n")
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
        return self._input_module.call_input_cabana(particle_count, filename)

    def gen_particle(self, particle):
        # Store the particle for later
        self._particle = particle
        # I have performance concerns about putting a struct inside the AoSoA
        # but it may be ok, we will see.
        # Output the core part type
        output = ""
        output = output + "struct core_part_type{\n"
        output = output + "    double position[3];\n"
        output = output + "    double velocity[3];\n"
        output = output + "};\n\n"

        #Output the neighbour part type
        output = output + "struct neighbour_part_type{\n"
        output = output + "    double cutoff;\n"
        output = output + "};\n\n"

        #particle.add_element("core_part_position", "double[3]")
#        particle.add_element("core_part_velocity", "double[3]")
#        particle.add_element("neighbour_part_cutoff", "double")

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
            if key != "core_part" and key != "neighbour_part":
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
        output = output + "\n               };\n"

        output = output + "using DataTypes = Cabana::MemberTypes<"
        first = True
        for key in sorted_fields:
            if key != "core_part" and key != "neighbour_part":
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
        output = output + ">;\n"

        return output

    def gen_config(self, config):
        # FIXME
        output = ""
        output = output + "struct boundary{\n"
        output = output + "    double x_min, x_max;\n"
        output = output + "    double y_min, y_max;\n"
        output = output + "    double z_min, z_max;\n"
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

    def cleanup(self, *args, current_indent, **kwargs):
        rval = "}\n"
        return rval

    def println(self, string, *args, **kwargs):
        '''
        Function to output the println string for the Cabana PIR module.
        Called via the visitors when reaching a println statement in user
        code.

        For the Cabana module this is a call to `std::cout`, using the
        string argument as the first value to be passed into cout, and
        the *args contains any further values used in the output.

        :param string: The formatted string to use with cout.                                                                                                                      :type string: str
        :param int current_indent: The current indentation level
        :param *args: A list of strings containing the other values
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

    def initialise(self,particle_count, filename, current_indent, **kwargs):
        # FIXME
        space = " "
        rval = space*current_indent + "Kokkos::ScopeGuard scope_guard(argc, argv);\n"
        rval = rval + "{\n"
        rval = rval + space*current_indent + "config_type config;\n"
        rval = rval + space*current_indent + "config.config = config_struct_type(\"config\", 1);\n"
        rval = rval + space*current_indent + "config.config_host = Kokkos::create_mirror_view(config.config);\n"
        rval = rval + space*current_indent + f"{self._input_module.call_input_cabana(int(particle_count), filename, current_indent=current_indent)}\n"
        rval = rval + space*current_indent + "Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());\n"

#        self.variable_scope = variable_scope()
        for struct in self._structures.keys():
            # Add this to the variable scope with a special c_type.
 #           self.variable_scope.add_variable(struct, self._structures[struct], False)
            rval = rval + space*current_indent + self._structures[struct] + " " + struct + ";\n"
        # Need to do something with each kernel now.
#        rval = rval + space*current_indent + "auto core_part_slice = Cabana::slice<core_part_space>(particle_aosoa);\n"
#        rval = rval + space*current_indent + "auto neighbour_part_slice = Cabana::slice<neighbour_part_space>(particle_aosoa);\n"

        # We need the particle type to be able to initialise correctly
        # The initialisation of the core part and neighbour cutoff stuf fis a bit weird, needs to be fixed.
        rval = rval + space*current_indent + "auto core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);\n"
        rval = rval + space*current_indent + "auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);\n"
        rval = rval + space*current_indent + "auto neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);\n"
        for key in self._particle.particle_type:
            if key != "core_part" and key != "neighbour_part":
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
                rval = rval + ","
            rval = rval + ", ".join(slice_names) + ");\n"
        return rval

    def add_coupler(self, coupled_system):
        '''
        Adds a coupled system to the list of coupled systems in the backend

        :param coupled_system: The object to couple with.
        :type couple_system: Object

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
        code = code + self._output_module.call_output_cabana(0, filename, variable) + "\n"
        return code

    def call_language_function(self, func_call: str, *args, **kwargs) -> str:
        string = ""
        try:
            fn = getattr(self, func_call)
            string = fn( *args, **kwargs)
            return string
        except (AttributeError) as err:
            pass
        for system in self._coupled_systems:
            print(system)
            try:
                fn = getattr(system, func_call)
                string = fn(*args, **kwargs)
                return string
            except (AttributeError) as err:
                pass
        raise AttributeError
