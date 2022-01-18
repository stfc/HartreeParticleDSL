import re
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
#import HartreeParticleDSL.backends.FDPS_backend.FDPS_visitors as FDPS_visitors
#from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.language_utils.variable_scope import variable_scope, \
                                                             variable_access

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
            c_bool : "bool"}

    # Variables used to manage where the cutoff radius is stored
    CONSTANT = 0
    PARTICLE = 1

    def __init__(self):
        self._pairwise_visitor = None #FIXME
        self._per_part_visitor = None #FIXME
        self._main_visitor = None
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

    @property
    def variable_scope(self):
        return self._variable_scope

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
            f.write(config_output)
            f.write(part_output)
            f.write("#endif")
        input_module_header = ""
        if self._input_module is not None:
            input_module_header = self._input_module.gen_code_cabana(part_type) #FIXME
        if input_module_header is not "":
            # Do something later
            pass
        output_module_header = ""
        if self._output_module is not None:
            output_module_header = self._output_module.gen_code_cabana(part_type) #FIXME
        if output_module_header is not "":
            # Do something later
            pass

    def gen_kernel(self, kernel):
        '''
        Generates the Cabana code for a given user-defined kernel.
        Used during code generation to output the kernel code in C++.

        :param kernel: Kernel object to generate code for.
        :type kernel: :py:class:`HartreeParticleDSL.kernel_types.kernels.kernel`
        '''
        tree = kernel.get_kernel_tree()
        self._variable_scope = variable_scope()
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
        pass #FIXME

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
        assert False
        # FIXME
        field = field.replace('"', '')
        return None

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
        assert False
        # FIXME
        if dimension == "x":
            return "core_part.position.x"
        if dimension == "y":
            return "core_part.position.y"
        if dimension == "z":
            return "core_part.position.z"
        raise InvalidNameError("The dimension argument should be x, y, or z")

    def get_pointer(self, var_code, *args, **kwargs):
        '''
        Returns the code to access the pointer to the supplied var_code.
        The var_code should already be Cabana C++ code.

        :param str var_code: The Cabana C++ code to take pointer from.

        :returns: The string pointer to the supplied var_code.
        :rtype: str
        '''
        return "&" + var_code

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        '''
        Set the cutoff for pairwise interactions. NYI

        :raises: NotImplementedError
        '''
        raise NotImplementedError("Cabana backend doesn't yet support pairwise interactions")

    def initialisation_code(self, particle_count, filename):
        return self._input_module.call_input_cabana(particle_count, filename)

    def gen_particle(self, particle):
        # FIXME
        assert False
        return None

    def gen_config(self, config):
        # FIXME
        assert False
        return None

    def cleanup(self, current_indent, *args, **kwargs):
        rval = ""
        assert False
#        rval = " "*current_indent + "PS::Finalize();\n" #FIXME
        return rval

    def initialise(self,particle_count, filename, current_indent, **kwargs):
        # FIXME
        assert False
        ravl = ""
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
        code_str = code_str + name
        array_access = (len(var_access.array_indices) != 0)
        check = check_valid and self._enable_variable_checks
        # Check for type existing
        if not var_access.is_child and check:
            if FDPS._type_map.get(var_access.variable.var_type) is None:
                raise UnsupportedTypeError("Accessing a variable of type "
                                          f"{var_access.variable.var_type} "
                                           "which is not supported by Cabana backend."
                                           f" Variable name is {var_access.variable.var_name}")
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
