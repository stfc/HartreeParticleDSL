from __future__ import annotations

from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str


class Cabana_PIR(Backend):
    '''
    Cabana Backend class using PIR.
    Outputs from this class use the Cabana & Kokkos system to compute
    parallel particle simulations.
    '''

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

        # The Cabana backend needs to store the particle type reference
        self._particle = None #TODO Remove & pull from PIR

        self._structures = {}

    def register_kernel(self, kernel_name: str, kernel: Kern):
        if kernel_name in self._kernels.keys():
            raise NotImplementedError()
        self._kernels[kernel_name] = kernel


    def add_structure(self, structure_type: StructureType, structure_name: str):
        assert False
        # I think this is handled by symbol tables now.
        # Needs to be done somehow

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
        if not isinstance(the_type, Datatype):
            raise UnsupportedTypeError(
                    f"Attempting to add a new type but got {type(the_type)} "
                    "but expected a Particle IR DataType instead.")

        type_mapping_str[type_name] = the_type

    def create_global_variable(self, c_type: str, name: str, initial_value: str):
        '''
        Function to create a global variable in the header file.

        :param str c_type: The c_type of this variable.
        :param str name: The name of the variable.
        :param initial_value: The initial value of this variable
        '''
        #TODO
        raise NotImplementedError()

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

            for name in self._globals:
                ctype = Cabana._type_map[self._globals[name][0]]
                value = self._globals[name][1]
                f.write("{0} {1} = {2};\n".format(ctype, name, value))

            f.write(config_output)
            f.write(part_output)
            f.write("#endif")
        input_module_header = ""
        if self._input_module is not None:
            input_module_header = self._input_module.gen_code_cabana(part_type) #FIXME
        if input_module_header is not "":
            print(input_module_header)
            print("\n")

        output_module_header = ""
        if self._output_module is not None:
            output_module_header = self._output_module.gen_code_cabana(part_type) #FIXME
        if output_module_header is not "":
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
        raise NotImplementedError()
        if isinstance(kernel, kernels.perpart_kernel_wrapper):
            print(self._per_part_visitor.visit(tree))
        elif isinstance(kernel, kernels.pairwise_kernel_wrapper):
            print(self._pairwise_visitor.visit(tree))

    def print_main(self, function):
        '''
        Generates the Cabana code for the supplied main function.

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        raise NotImplementedError()

    def gen_invoke_perpart(self, kernel_name, current_indent, indent, kernel_type):
        '''
        Generates a Cabana compatible per_part kernel
        '''
        raise NotImplementedError()

    def gen_invoke(self, kernel_name, current_indent, indent):
        '''
        '''
        raise NotImplementedError("gen_invoke not yet implemented")

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        '''
        Set the cutoff for pairwise interactions. NYI

        :raises: NotImplementedError
        '''
        raise NotImplementedError("Cabana PIR backend doesn't yet support pairwise interactions")

    def initialisation_code(self, particle_count, filename):
        return self._input_module.call_input_cabana(particle_count, filename)

    def gen_particle(self, particle):
        raise NotImplementedError()

    def gen_config(self, config):
        raise NotImplementedError()

    def cleanup(self, current_indent, *args, **kwargs):
        rval = "}\n"
        return rval

    def initialise(self,particle_count, filename, current_indent, **kwargs):
        # FIXME
        space = " "
        rval = space*current_indent + "Kokkos::ScopeGuard scope_guard(argc, argv);\n"
        rval = rval + "{\n"
        rval = rval + space*current_indent + "config_type config;\n"
        rval = rval + space*current_indent + "config.config = config_struct_type(\"config\", 1);\n"
        rval = rval + space*current_indent + "config.config_host = Kokkos::create_mirror_view(config.config);\n"
        rval = rval + space*current_indent + f"{self._input_module.call_input_cabana(particle_count, filename, current_indent=current_indent)}\n"
        rval = rval + space*current_indent + "Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());\n"

        self.variable_scope = variable_scope()
        for struct in self._structures.keys():
            # Add this to the variable scope with a special c_type.
            self.variable_scope.add_variable(struct, self._structures[struct], False)
            rval = rval + space*current_indent + self._structures[struct] + " " + struct + ";\n"
        # Need to do something with each kernel now.
#        rval = rval + space*current_indent + "auto core_part_slice = Cabana::slice<core_part_space>(particle_aosoa);\n"
#        rval = rval + space*current_indent + "auto neighbour_part_slice = Cabana::slice<neighbour_part_space>(particle_aosoa);\n"

        # We need the particle type to be able to initialise correctly
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
            string = fn(*args, **kwargs)
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
