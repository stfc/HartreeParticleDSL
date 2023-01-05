import re
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
import HartreeParticleDSL.backends.Cabana_backend.Cabana_visitors as Cabana_visitors
from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.language_utils.variable_scope import variable_scope, \
                                                             variable_access, \
                                                             variable

class Cabana(Backend):
    '''
    Cabana Backend class.
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
            "FULLPART" : "FULLPART", # This is a proxy type, used to convert to slices when required
            "FULLCONF" : "FULLCONF" # This is a proxy type, used to handle the config type for later
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

    def __init__(self):
        self._pairwise_visitor = Cabana_visitors.cabana_pairwise_visitor(self) #FIXME
        self._per_part_visitor = Cabana_visitors.cabana_perpart_visitor(self)
        self._main_visitor = Cabana_visitors.cabana_main_visitor(self)
        self._includes = []
        self._includes.append("<Cabana_Core.hpp>")
        self._includes.append("<Kokkos_Core.hpp>")
        self._includes.append("<Kokkos_Random.hpp>")
        self._includes.append("<iostream>")
        self._includes.append("<cmath>")
        self._includes.append("<cstdio>")
        self._includes.append("\"part.h\"") # FIXME
        self._cutoff_type = Cabana.PARTICLE #FIXME
        self._cutoff = "neighbour_part.cutoff" #FIXME
        self._input_module = None
        self._output_module = None
        self._variable_scope = variable_scope()
        self._enable_variable_checks = True
        self._coupled_systems = []
        self._in_kernel_code = False

        # We need to keep track of what kernels use what slices
        self._kernel_slices = {}

        # Allowing global variable declaration for now
        self._globals = {}
        # The Cabana backend needs to store the particle type reference
        self._particle = None

        self._structures = {}

    def add_structure(self, structure_type, structure_name):
        if structure_type not in  Cabana._type_map.keys():
            raise UnsupportedTypeError("Cabana backend doesn't support {0}"
            " as a type for structure.".format(structure_type))
        self._structures[structure_name] = structure_type

    def add_kernel_slice(self, kernel_name, kernel_slice):
        '''
        Adds access to a particle slice for a specific kernel.
        '''
        if kernel_name not in self._kernel_slices.keys():
            self._kernel_slices[kernel_name] = []
        if kernel_slice not in self._kernel_slices[kernel_name]:
            self._kernel_slices[kernel_name].append(kernel_slice)

    @property
    def in_kernel_code(self):
        return self._in_kernel_code

    @in_kernel_code.setter
    def in_kernel_code(self, in_kern):
        assert isinstance(in_kern, bool)
        self._in_kernel_code = in_kern

    @property
    def variable_scope(self):
        return self._variable_scope

    @variable_scope.setter
    def variable_scope(self, variable_scope):
        self._variable_scope = variable_scope

    def disable_variable_checks(self):
        '''
        Disables validity checking when creating strings
        from variable accesses. Usually used during internals
        and tests when non-Cabana variables can appear.
        '''
        self._enable_variable_checks = False

    def enable_variable_checks(self):
        '''
        Enables validity checking when creating strings
        from variable accesses. This is the default state, and
        is disabled occasionally internally or during testing.
        '''
        self._enable_variable_checks = True

    def add_type(self, type_name, type_string):
        '''
        Function to add a special type to the type map for this backend.

        :param str type_name: The name of the type (this can just be the type_string)
        :param str type_string: The Cabana type string to use for this type
        '''
        Cabana._type_map[type_name] = type_string

    def create_global_variable(self, c_type, name, initial_value):
        '''
        Function to create a global variable in the header file.

        :param str c_type: The c_type of this variable.
        :param str name: The name of the variable.
        :param initial_value: The initial value of this variable
        '''
        if name in self._globals.keys():
            raise InvalidNameError(f"{name} already declared as a global variable.")
        self._globals[name] = (c_type, initial_value)

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

    def println(self, string, *args, **kwargs):
        '''
        Function to output the println string for the Caba module.
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
        tree = kernel.get_kernel_tree()
        self._variable_scope = variable_scope()
        for name in self._globals:
            ctype = self._globals[name][0]
            self._variable_scope.add_variable(name, ctype, False)
        if isinstance(kernel, kernels.perpart_kernel_wrapper):
            print(self._per_part_visitor.visit(tree))
        elif isinstance(kernel, kernels.pairwise_kernel_wrapper):
            print(self._pairwise_visitor.visit(tree))

    def print_main(self, function):
        '''
        Generates the Cabana code for the supplied main function.
        Currently outputs to STDOUT

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        print(self._main_visitor.visit(function))

    def gen_invoke_perpart(self, kernel_name, current_indent, indent, kernel_type):
        '''
        Generates a Cabana compatible per_part kernel
        '''
        space = " "
        rval = "" + space*current_indent
        rval = rval + "Kokkos::deep_copy(config.config, config.config_host);\n"
        # If we have external structures make sure they're copied to the functor as well
        if len(self._structures) > 0:
            rval = rval + space*current_indent + f"{kernel_name}.update_structs("
            struct_list = []
            for struct in self._structures:
                struct_list.append(struct)
            rval = rval + ", ".join(struct_list) + ");\n"
        rval = rval + space*current_indent + f"Cabana::simd_parallel_for(simd_policy, {kernel_name}, " + "\"" + kernel_name + "\");\n"
        rval = rval + space*current_indent + "Kokkos::fence();\n"
        return rval

    def gen_invoke(self, kernel_name, current_indent, indent, kernel_type):
        '''
        '''
        space = " "
        rval = f"\n {space*current_indent}/* INVOKE generated for {kernel_name} */\n"
        if kernel_type == kernels.perpart_kernel_wrapper:
            rval = rval + self.gen_invoke_perpart(kernel_name, current_indent, indent, kernel_type)
        else:
            raise NotImplementedError("gen_invoke not yet implemented")
        rval = rval + f"{space*current_indent}/* End of INVOKE generated for {kernel_name} */\n\n"
        return rval

    def get_particle_access(self, index, field):
        '''
        Returns the code to access a particle of the given
        index for the field provided.

        :param str index: The index value to access.
        :param str field: The field name to access.

        TODO: We could check the field exists
        '''
        # Remove any extra " from the field from passing through the DSL
        # FIXME
        field = field.replace('"', '')
        
        # Remove all the array indices and do something with them
        extra_indices = ""
        arrays = re.findall(r"\[[0-9*]*\]", field)
        if arrays is not None:
            for ind in arrays:
                field = field.replace(ind, "")
                ind = ind.replace("[", "")
                ind = ind.replace("]", "")
                extra_indices = extra_indices + ", " + ind

        # Create a variable access is probably the only sane way to do this.
        self._pairwise_visitor.addSlice(field)
        self._per_part_visitor.addSlice(field)
        self._main_visitor.addSlice(field)
        assert index == "part1" # FIXME Not handling part2 accesses for pairwise yet.

        return "_" + field + ".access(i, a" + extra_indices + ")"

    def _get_particle_position_internal(self, dimension):
        '''
        Returns the index corresponding to a dimension accessed.
        For Cabana, x -> 0, y-> 1, z->2
        '''
        if dimension == "x":
            return "0"
        if dimension == "y":
            return "1"
        if dimension == "z":
            return "2"
        raise InvalidNameError("The dimension argument should be x, y, or z")

    def get_particle_position(self, dimension):
        '''
        Returns the code to access a particle's position
        for each dimension. Dimensions are x/y/z. For FDPS
        the positions are stored in a PS::F64vec, so we return
        the relevant vector element

        :param str dimension: The dimension ("x", "y" or "z") to access

        :raises InvalidNameError: If the dimension argument is not
                                  "x", "y" or "z".

        :returns: The string to access a particle's position variable.
        :rtype: str
        '''
        # FIXME
        if dimension == "x":
            return "core_part_position[0]"
        if dimension == "y":
            return "core_part_position[1]"
        if dimension == "z":
            return "core_part_position[2]"
        raise InvalidNameError("The dimension argument should be x, y, or z")

    def get_pointer(self, var_code, *args, **kwargs):
        '''
        Returns the code to access the pointer to the supplied var_code.
        The var_code should already be Cabana C++ code.

        :param str var_code: The Cabana C++ code to take pointer from.

        :returns: The string pointer to the supplied var_code.
        :rtype: str
        '''
        return "&(" + var_code + ")"

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        '''
        Set the cutoff for pairwise interactions. NYI

        :raises: NotImplementedError
        '''
        raise NotImplementedError("Cabana backend doesn't yet support pairwise interactions")

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

        particle.add_element("core_part_position", "double[3]")
        particle.add_element("core_part_velocity", "double[3]")
        particle.add_element("neighbour_part_cutoff", "double")

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
                size = Cabana._type_sizes[base_type]
                sizes.append(size * count)
            else:
                size = Cabana._type_sizes[c_type]
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

    def create_variable(self, c_type, name, initial_value=None, **kwargs):
        current_indent = kwargs.get("current_indent", 0)
        if Cabana._type_map.get(c_type) is None:
            raise UnsupportedTypeError("Cabana does not support type \"{0}\""
                                        " in created variables.".format(c_type))
        ##Check name is allowed in C++
        name = name.replace('"', '')
        a = re.match("[a-zA-Z_][a-zA-Z_0-9]*", name)
        if a is None or a.group(0) != name:
            raise InvalidNameError("Cabana does not support \"{0}\" as a name"
                                   " for variables.".format(name))
        end = ";\n"
        if initial_value is not None:
            end = f" = {initial_value};\n"
        rval = " " * current_indent + Cabana._type_map.get(c_type) + " " + name + end
        self.variable_scope.add_variable(name, c_type, False)
        return rval

    def call_language_function(self,func_call, *args, **kwargs):
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

        string = func_call + "( "
        arguments = []
        for arg in args:
            arguments.append(arg)
        for kwarg in kwargs:
            arguments.append(f"{kwarg}={kwargs[kwarg]}")
        arg_string = ", ".join(arguments)
        string = string + arg_string + " )"
        current_index = re.search("current_indent=[0-9]*", string)
        current_indent = 0
        if current_index is not None:
            current_indent = int(current_index[0][15:])
        string = " "*current_indent + string
        string = re.sub(", current_indent=[0-9]*, indent=[0-9]*", "", string)
        string = re.sub(" current_indent=[0-9]*, indent=[0-9]*", "", string)
        string = string + ";\n"
        return string

    def access_to_string(self, var_access, check_valid=False):
        '''
        Takes a variable_access and converts it to a Cabana string to output

        :param var_access: The variable access to output as a string
        :type var_access: variable_access

        :raises UnsupportedTypeError: If the variable is of a type not supported
                                      by Cabana and not a child of a structure.

        :returns: The Cabana string for this variable access
        :rtype: str
        '''
        code_str = ""
        name = var_access.variable.var_name
        array_access = (len(var_access.array_indices) != 0)
        check = check_valid and self._enable_variable_checks
        # Check for type existing
        if not var_access.is_child and check:
            if Cabana._type_map.get(var_access.variable.var_type) is None:
                raise UnsupportedTypeError("Accessing a variable of type "
                                          f"{var_access.variable.var_type} "
                                           "which is not supported by Cabana backend."
                                           f" Variable name is {var_access.variable.var_name}")
        # Do something special if its of type FULLPART
        if var_access.variable.var_type == "FULLPART":
            #FIXME
            is_part_i = False
            if name == self._pairwise_visitor._part1 or name == self._per_part_visitor._part1:
                is_part_i = True
            # For now only support accessing part i
            if not is_part_i:
                raise InternalError("Attempting to access a particle other than part i")
            child = var_access.child
            if child is None:
                raise InternalError("Attempting to access a particle with no child access")
            name = child.variable.var_name
            if name == "core_part":
                # If its core part we have to do something special
                var_access = child
                child_name = var_access.child.variable.var_name
                name = name + "_" + child_name
                var_access = var_access.child
                if is_part_i:
                    code_str = "_" + name + ".access(i, a"
                #else:
                #    pass #FIXME for pairwise we do something different
                if child_name == "position":
                    code_str = code_str + ", " + self._get_particle_position_internal( var_access.child.variable.var_name)
                if len(var_access.array_indices) > 0:
                    for index in var_access.array_indices:
                        if isinstance(index, str):
                            code_str = code_str + f", {index}"
                        if isinstance(index, variable_access):
                            code_str = code_str + f", " + self.access_to_string(index)
                # Nothing else can happen here
                self._pairwise_visitor.addSlice(name)
                self._per_part_visitor.addSlice(name)
                self._main_visitor.addSlice(name)
                return code_str + ")"
            if name == "neighbour_part":
                # If its neighbour_part we have to do something special
                var_access = child
                child_name = var_access.child.variable.var_name
                name = name + "_" + child_name
                var_access = var_access.child
                if is_part_i:
                    code_str = "_" + name + ".access(i, a"
                #else:
                #    pass #FIXME for pairwise we do something different
                if len(var_access.array_indices) > 0:
                    for index in var_access.array_indices:
                        if isinstance(index, str):
                            code_str = code_str + f", {index}"
                        if isinstance(index, variable_access):
                            code_str = code_str + f", " + self.access_to_string(index)
                self._pairwise_visitor.addSlice(name)
                self._per_part_visitor.addSlice(name)
                self._main_visitor.addSlice(name)
                # Nothing else can happen here
                return code_str + ")"
            self._pairwise_visitor.addSlice(name)
            self._per_part_visitor.addSlice(name)
            self._main_visitor.addSlice(name)
            code_str = ""
            if is_part_i:
                code_str = "_" + name + ".access(i, a"
            #else:
            #    pass #FIXME for pairwise we do something different
            if array_access:
                for index in var_access.array_indices:
                    if isinstance(index, str):
                        code_str = code_str + f", {index}"
                    if isinstance(index, variable_access):
                        code_str = code_str + f", " + self.access_to_string(index)
            if var_access.child is not None:
#                if var_access.variable.is_pointer and not array_access:
#                    child = var_access.child
#                    child_str = self.access_to_string(child)
#                    code_str = code_str + "->" + child_str
                if len(var_access.child.array_indices) > 0:
                    for index in var_access.child.array_indices:
                        if isinstance(index, str):
                            code_str = code_str + f", {index}"
                        if isinstance(index, variable_access):
                            code_str = code_str + f", " + self.access_to_string(index)
#                else:
#                    child = var_access.child
#                    child_str = self.access_to_string(child)
#                    code_str = code_str + "." + child_str
            code_str = code_str + ")"
            return code_str
        # Do something special if its of type FULLCONF
        # FIXME Do something different if in main vs kernel?
        if var_access.variable.var_type == "FULLCONF":
            name = ""
            code_str = ""
            if self._in_kernel_code:
                if self._pairwise_visitor._config == "":
                    name = self._per_part_visitor._config
                else:
                    name = self._pairwise_visitor._config
                code_str = "_" + name + "(0)"
            else:
                name = "config.config_host"
                code_str = name + "(0)"
            if array_access:
                raise InternalError("Currently can't handle array accesses for config type in Cabana backend")
#                for index in var_access.array_indices:
#                    if isinstance(index, str):
#                        code_str = code_str + f", {index}"
#                    if isinstance(index, variable_access):
#                        code_str = code_str + f", " + self.access_to_string(index)
#            code_str = code_str + ")"
            if var_access.variable.is_pointer and not array_access:
                child = var_access.child
                child_str = self.access_to_string(child)
                code_str = code_str + "->" + child_str
            else:
                child = var_access.child
                child_str = self.access_to_string(child)
                code_str = code_str + "." + child_str
            return code_str

        # Otherwise, generate a normal access
        code_str = code_str + name
        if array_access:
            for index in var_access.array_indices:
                if isinstance(index, str):
                    code_str = code_str + f"[{index}]"
                if isinstance(index, variable_access):
                    code_str = code_str + "[" + self.access_to_string(index) + "]"
        if var_access.child is not None:
            if var_access.variable.is_pointer and not array_access:
                child = var_access.child
                child_str = self.access_to_string(child)
                code_str = code_str + "->" + child_str
            else:
                child = var_access.child
                child_str = self.access_to_string(child)
                code_str = code_str + "." + child_str

        return code_str

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

    def write_output(self, filename, variable=None, **kwargs):
        '''
        Generates the code to write a file output using the selected output module.

        :param str filename: The filename to write the file to.

        '''
        current_indent = kwargs.get("current_indent", 0)
        code = " " * current_indent
        code = code + self._output_module.call_output_cabana(0, filename, variable) + "\n"
        return code
