import re
from HartreeParticleDSL.backends.base_backend.backend import Backend
import HartreeParticleDSL.backends.FDPS_backend.FDPS_visitors as FDPS_visitors
from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.language_utils.variable_scope import variable_scope, \
                                                             variable_access

class FDPS(Backend):
    '''
    FDPS Backend class.
    Outputs from this class use the FDPS system to compute parallel
    particle simulations.
    '''
    _type_map = {c_int : "PS::S32",
                 c_double : "PS::F64",
                 c_float : "PS::F32",
                 c_int64_t : "PS::S64",
                 c_int32_t : "PS::S32",
                 c_int8_t : "char",
                 c_bool : "bool",
                 "FullParticle" : "FullParticle",
                 "config_type" : "config_type"}

    # Variables used to manage where the cutoff radius is stored
    CONSTANT = 0
    PARTICLE = 1

    def __init__(self):
        self._pairwise_visitor = FDPS_visitors.fdps_pairwise_visitor(self)
        self._per_part_visitor = FDPS_visitors.fdps_perpart_visitor(self)
        self._main_visitor = FDPS_visitors.fdps_main_visitor(self)
        self._includes = []
        self._includes.append("<cmath>")
        self._includes.append("<cstdio>")
        self._includes.append("<iostream>")
        self._includes.append("<vector>")
        self._includes.append("<particle_simulator.hpp>")
        self._includes.append("\"part.h\"")
        self._cutoff_type = FDPS.PARTICLE
        self._cutoff = "neighbour_part.cutoff"
        self._input_module = None
        self._output_module = None
        self._variable_scope = variable_scope()
        self._enable_variable_checks = True

    @property
    def variable_scope(self):
        return self._variable_scope

    def disable_variable_checks(self):
        self._enable_variable_checks = False

    def enable_variable_checks(self):
        self._enable_variable_checks = True

    def set_io_modules(self, input_module, output_module):
        '''
        Function to set the IO Modules. Called from HartreeParticleDSL
        helper functions instead of directly by the user.

        :param input_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.FDPS_backend.XCXX` \
                            or None
        :param ouput_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.FDPS_backend.XCXX` \
                             or None
        :raises InvalidIOModuleError: if either input_module or output_module \
                                      are not implementing the XCXXXX \
                                      class
        '''
        # TODO
        if input_module is not None and not isinstance(input_module, FDPS_IO_Mixin):
            raise InvalidIOModuleError(f"{self.__class__.__name__} backend does not support "
                                       f"{input_module.__class__.__name__} IO "
                                       f"module at this time.")
        if output_module is not None and not isinstance(output_module, FDPS_IO_Mixin):
            raise InvalidIOModuleError(f"{self.__class__.__name__} backend does not support "
                                       f"{output_module.__class__.__name__} IO "
                                       f"module at this time.")
        self._input_module = input_module
        self._output_module = output_module

    def println(self, string, *args, **kwargs):
        '''
        Function to output the println string for the FDPS module.
        Called via the visitors when reaching a println statement in user
        code.

        For the FDPS module this is a call to `std::cout`, using the
        string argument as the first value to be passed into cout, and
        the *args contains any further values used in the output.

        :param string: The formatted string to use with cout.
        :type string: str
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
        Function to add an include to the FDPS include list.

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
            input_includes = self._input_module.get_includes_fdps()
            for include in input_includes:
                if include not in self._includes:
                    self._includes.append(include)
        if self._output_module is not None:
            output_includes = self._output_module.get_includes_fdps()
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
            input_module_header = self._input_module.gen_code_fdps(part_type)
        if input_module_header is not "":
            # Do something later
            pass
        output_module_header = ""
        if self._output_module is not None:
            output_module_header = self._output_module.gen_code_fdps(part_type)
        if output_module_header is not "":
            # Do something later
            pass

    def gen_kernel(self, kernel):
        '''
        Generates the FDPS code for a given user-defined kernel.
        Used during code generation to output the kernel code in C.

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
        Generates the FDPS code for the supplied main function.
        Currently outputs to STDOUT

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        print(self._main_visitor.visit(function))

    def gen_invoke_perpart(self, kernel_name, current_indent, indent, kernel_type):
        '''
        Generates an FDPS compatible per_part kernel
        '''
        space = " "
        rval = space*current_indent
        rval = rval + "for(PS::S32 i = 0; i < particle_system.getNumberOfParticleLocal(); i++){\n"
        current_indent = current_indent + indent
        rval = rval + space*current_indent
        rval = rval + f"{kernel_name}(particle_system[i], config);\n"
        current_indent = current_indent - indent
        rval = rval + space*current_indent + "}\n"
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

    def get_particle_access(self, dimension):
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
        The var_code should already be FDPS C++ code.

        :param str var_code: The FDPS C++ code to take pointer from.

        :returns: The string pointer to the supplied var_code.
        :rtype: str
        '''
        return "&" + var_code

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        '''
        Set the cutoff for pairwise interactions. NYI

        :raises: NotImplementedError
        '''
        raise NotImplementedError("FDPS backend doesn't yet support pairwise interactions")

    def initialisation_code(self, particle_count, filename):
        return self._input_module.call_input_fdps(particle_count, filename)

    def gen_particle(self, particle):
        # Output the core part type
        output = ""
        output = output + "struct core_part_type{\n"
        output = output + "    PS::F64vec position;\n"
        output = output + "    PS::F64 velocity[3];\n"
        output = output + "};\n\n"

        # Output the neighbour part type
        # Currently this is empty
        output = output + "struct neighbour_part_type{\n"
        output = output + "    PS::F64 cutoff;\n"
        output = output + "};\n\n"

        output = output + "class FullParticle{\n"
        output = output + "    public:\n"
        for key in particle.particle_type:
            is_array = particle.particle_type[key]['is_array']
            c_type = particle.particle_type[key]['type']
            varnam = key
            if is_array:
                # Move array index to the correct position for C++ structs
                x = c_type.index("[")
                varnam = varnam + c_type[x:]
                c_type = c_type[0:x]
            output = output + f"        {c_type} {varnam};\n"
        output = output + "        PS::F64vec getPos(){\n"
        output = output + "            return this->core_part.position;\n"
        output = output + "        }\n"
        output = output + "        void setPos(const PS::F64vec pos_new){\n"
        output = output + "            this->core_part.position = pos_new;\n"
        output = output + "        }\n"
        output = output + "        void clear(){\n"
        output = output + "        }\n"
        output = output + "};\n\n"
        return output

    def gen_config(self, config):
        # Output the space
        output = ""
        output = output + "struct boundary{\n"
        output = output + "    PS::F64 x_min, x_max;\n"
        output = output + "    PS::F64 y_min, y_max;\n"
        output = output + "    PS::F64 z_min, z_max;\n"
        output = output + "};\n\n"
        output = output + "struct space_type{\n"
        output = output + "    boundary box_dims;\n"
        output = output + "    PS::S32 nparts;\n"
        output = output + "};\n\n"

        # Currently empty neiughbour config
        output = output + "struct neighbour_config_type{\n"
        output = output + "};\n\n"


        output = output + "class config_type{\n"
        output = output + "    public:\n"
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
        return output

    def cleanup(self, current_indent, *args, **kwargs):
        rval = ""
        rval = " "*current_indent + "PS::Finalize();\n"
        return rval

    def initialise(self,particle_count, filename, current_indent, **kwargs):
        rval = " "*current_indent + "char **argv = NULL;\n"
        rval = rval + " "*current_indent + "int args = 0;\n"
        rval = rval + " "*current_indent + "PS::Initialize(args,argv);\n"
        rval = rval + " "*current_indent + "config_type config;\n"
        rval = rval + " "*current_indent + f"{self._input_module.call_input_fdps(particle_count, filename, current_indent=current_indent)}\n"
        return rval

    def create_variable(self, c_type, name, initial_value=None, **kwargs):
        current_indent = kwargs.get("current_indent", 0)
        if FDPS._type_map.get(c_type) is None:
            raise UnsupportedTypeError("FDPS does not support type \"{0}\""
                                        " in created variables.".format(c_type))
        ##Check name is allowed in C++
        name = name.replace('"', '')
        a = re.match("[a-zA-Z_][a-zA-Z_0-9]*", name)
        if a is None or a.group(0) != name:
            raise InvalidNameError("FDPS does not support \"{0}\" as a name"
                                   " for variables.".format(name))
        end = ";\n"
        if initial_value is not None:
            end = f" = {initial_value};\n"
        rval = " " * current_indent + FDPS._type_map.get(c_type) + " " + name + end
        self.variable_scope.add_variable(name, c_type, False)
        return rval

    def call_language_function(self,func_call, *args, **kwargs):
        string = ""
        try:
            fn = getattr(self, func_call)
            string = fn(*args, **kwargs)
        except (SyntaxError, TypeError, AttributeError) as err:
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
        Takes a variable_access and converts it to a FDPS string to output

        :param var_access: The variable access to output as a string
        :type var_access: variable_access

        :raises UnsupportedTypeError: If the variable is of a type not supported
                                      by FDPS and not a child of a structure.

        :returns: The FDPS string for this variable access
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
                                           "which is not supported by FDPS backend."
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
