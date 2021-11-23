import re
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import base_coupler
import HartreeParticleDSL.backends.C_AOS.visitors as c_visitors
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError
from HartreeParticleDSL.c_types import c_int, c_double, c_float, c_int64_t, \
                                       c_int32_t, c_int8_t, c_bool
from HartreeParticleDSL.language_utils.variable_scope import variable_scope, \
                                                             variable_access


class C_AOS(Backend):
    '''
    C_AOS Backend class.
    Outputs from this class are Serial C code with Arrays of Struct particle
    representation.
    '''
    _type_map = {c_int : "int",
                 c_double : "double",
                 c_float : "float",
                 c_int64_t : "long long int",
                 c_int32_t : "int",
                 c_int8_t : "char",
                 c_bool : "_Bool",
                 "struct part": "struct part",
                 "struct config_type": "struct config_type"}

    #Variables used to manage where the cutoff radius is stored
    CONSTANT = "C_AOS.CONSTANT"
    PARTICLE = "C_AOS.PARTICLE"

    def __init__(self):
        self._pairwise_visitor = c_visitors.c_pairwise_visitor(self)
        self._per_part_visitor = c_visitors.c_perpart_visitor(self)
        self._main_visitor = c_visitors.c_main_visitor(self)
        self._includes = []
        self._includes.append("<math.h>")
        self._includes.append("<stdio.h>")
        self._includes.append("\"part.h\"")
        self._cutoff_type = C_AOS.CONSTANT
        self._cutoff = "config->neighbour_config.cutoff"
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
        and tests when non-C_AOS variables can appear.
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
        :param str type_string: The C type string to use for this type
        '''
        C_AOS._type_map[type_name] = type_string

    def set_io_modules(self, input_module, output_module):
        '''
        Function to set the IO Modules. Called from HartreeParticleDSL
        helper functions instead of directly by the user.

        :param input_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin.C_AOS_IO_Mixin`
        :param ouput_module: The IO module to use for input from file.
        :type input_module: :py:class:`HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin.C_AOS_IO_Mixin` \
                             or None
        :raises InvalidIOModuleError: if either input_module or output_module \
                                      are not implementing the C_AOS_IO_Mixin \
                                      class
        '''
        if input_module is not None and not isinstance(input_module, C_AOS_IO_Mixin):
            raise InvalidIOModuleError(f"C_AOS backend does not support "
                                       f"{input_module.__class__.__name__} IO "
                                       f"module at this time.")
        if output_module is not None and not isinstance(output_module, C_AOS_IO_Mixin):
            raise InvalidIOModuleError(f"C_AOS backend does not support "
                                       f"{output_module.__class__.__name__} IO "
                                       f"module at this time.")
        self._input_module = input_module
        self._output_module = output_module

    def println(self, string, *args, **kwargs):
        '''
        Function to output the println string for the C_AOS module.
        Called via the visitors when reaching a println statement in user
        code.

        For the C_AOS module this is a call to `printf`, using the 
        string argument as the formatted string with a line end, and
        the *args contains any C values used in the formatted string.

        :param string: The formatted string to use with printf.
        :type string: str
        :param int current_indent: The current indentation level
        :param *args: A list of strings containing the C values required
                      for output.
        :type args: str
        '''
        current_indent = kwargs.get("current_indent", 0)
        x = string[0].replace('"', '') + string[1:]
        x = x[0:len(x)-1] + x[-1].replace('"', '')
        output = " "*current_indent + f"printf(\"{x}\\n\""
        for arg in args:
            output = output + f", {arg}"
        output = output + ");\n"
        return output

    def add_include(self, include_string):
        '''
        Function to add an include to the C_AOS include list.

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
            input_includes = self._input_module.get_includes_c()
            for include in input_includes:
                if include not in self._includes:
                    self._includes.append(include)
        if self._output_module is not None:
            output_includes = self._output_module.get_includes_c()
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
            input_module_header = self._input_module.gen_code_c(part_type)
        if input_module_header is not "":
            # Do something later
            pass
        output_module_header = ""
        if self._output_module is not None:
            output_module_header = self._output_module.gen_code_c(part_type)
        if output_module_header is not "":
            # Do something later
            pass

    def gen_kernel(self, kernel):
        '''
        Generates the C_AOS code for a given user-defined kernel.
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
        Generates the C_AOS code for the supplied main function.
        Currently outputs to STDOUT

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        self._variable_scope = variable_scope()
        print(self._main_visitor.visit(function))

    def gen_invoke(self, kernel_name, current_indent, indent, kernel_type):
        '''
        Generates the code for a specific invoke call
        Currently ouputs to STDOUT

        :param str kernel_name: Name of the kernel name to be generated.
        :param int current_indent: The current indentation level of code.
        :param int indent: The amount to indent for further indentation levels.
        :param kernel_type: The type of the kernel to be invoked.
        :type kernel_type: type
        '''
        space = " "
        rval = f"\n {space*current_indent}/* INVOKE generated for {kernel_name} */\n"
        if kernel_type == kernels.perpart_kernel_wrapper:
            rval = rval + " "*current_indent
            rval = rval + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
            current_indent = current_indent + indent
            rval = rval + " "*current_indent
            rval = rval + f"{kernel_name}(&parts[part1], config);\n"
            current_indent = current_indent - indent
            rval = rval + " "*current_indent
            rval = rval + "}\n"
        elif kernel_type == kernels.pairwise_kernel_wrapper:
            rval = rval + " "*current_indent
            rval = rval + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
            current_indent = current_indent + indent
            rval = rval + " "*current_indent
            rval = rval + "for( int part2 = 0; part2 < config->space.nparts; part2++){\n"
            current_indent = current_indent + indent
            rval = rval + " "*current_indent
            rval = rval + "if(part1 == part2) continue;\n"
            rval = rval + " "*current_indent
            rval = rval + "double r2 = 0.0;\n"
            rval = rval + " "*current_indent
            rval = rval + "for(int k = 0; k < 3; k++){\n"
            current_indent = current_indent + indent
            rval = rval +" "*current_indent
            rval = rval + "r2 += (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]) * (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]);\n"
            current_indent = current_indent - indent
            rval = rval + " "*current_indent
            rval = rval + "}\n"
            rval = rval + " "*current_indent
            rval = rval + f"if(r2 < ({self._cutoff} * {self._cutoff}))" + "{\n"
            current_indent = current_indent + indent
            rval = rval + " "*current_indent
            rval = rval + f"{kernel_name}(&parts[part1], &parts[part2], r2, config);\n"
            current_indent = current_indent - indent
            rval = rval + " "*current_indent
            rval = rval + "}\n"
            current_indent = current_indent - indent
            rval = rval + " "*current_indent
            rval = rval + "}\n"
            current_indent = current_indent - indent
            rval = rval + " "*current_indent
            rval = rval + "}\n"
        rval = rval + f"{space*current_indent}/* End of INVOKE generated for {kernel_name} */\n\n"
        return rval

    def initialisation_code(self, particle_count, filename):
        '''
        Returns the initialisation code for this system.

        :param int particle_count: The particle count required for this job.
        :param str filename: The filename required for this job

        :returns: The initialisation code for this system.
        :rtype: str
        '''
        return self._input_module.call_input_c(particle_count, filename)

    def gen_particle(self, particle):
        '''
        Returns the code to setup the particle type in this module.

        :param particle: The particle type to be generated as code.
        :type particle: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Particle` 

        :returns: The particle code for this system.
        :rtype: str
        '''
        # Output the core part type
        output = ""
        output = output + "struct core_part_type{\n"
        output = output + "    double position[3];\n"
        output = output + "    double velocity[3];\n"
        output = output + "};\n\n"

        # Output the neighbour part type
        # Currently this is empty
        output = output + "struct neighbour_part_type{\n"
        output = output + "};\n\n"

        output = output + "struct part{\n"
        for key in particle.particle_type:
            is_array = particle.particle_type[key]['is_array']
            c_type = particle.particle_type[key]['type']
            varnam = key
            if is_array:
                # Move array index to the correct position for C structs
                x = c_type.index("[")
                varnam = varnam + c_type[x:]
                c_type = c_type[0:x]
            output = output + f"    {c_type} {varnam};\n"
        output = output + "};\n\n"
        return output

    def gen_config(self, config):
        '''
        Returns the code to setup the config type in this module.

        :param particle: The particle type to be generated as code.
        :type particle: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Config` 

        :returns: The config code for this system.
        :rtype: str
        '''
        # Output the space
        output = ""
        output = output + "struct space_type{\n"
        output = output + "    double box_dims[3];\n"
        output = output + "    int nparts;\n"
        output = output + "};\n\n"

        # Currently empty neiughbour config
        output = output + "struct neighbour_config_type{\n"
        output = output + "    double cutoff;\n"
        output = output + "};\n\n"

        output = output + "struct config_type{\n"
        for key in config.config_type:
            is_array = config.config_type[key]['is_array']
            c_type = config.config_type[key]['type']
            varnam = key
            if is_array:
                # Move array index to the correct position for C structs
                x = c_type.index("[")
                varnam = varnam + c_type[x:]
                c_type = c_type[0:x]
            output = output + f"    {c_type} {varnam};\n"
        output = output + "};\n\n"
        return output

    def cleanup(self, current_indent, *args, **kwargs):
        '''
        Returns the cleanup code for this system.

        :param int current_indent: The current indentation level

        :returns: Cleanup code for this system.
        :rtype: str
        '''
        rval = ""
        rval = " "*current_indent + "free(config);\n"
        rval = rval + " "*current_indent + "free(parts);\n"
        return rval

    def initialise(self,particle_count, filename, current_indent, **kwargs):
        '''
        Returns the initialisation code for this sytem.

        :param int particle_count: The particle count for this system.
        :param str filename: The filename for this system.
        :param int current_indent: The current indentation for this system.

        :returns: Initialisation code for this system.
        :rtype: str
        '''
        self.variable_scope.add_variable("config", "struct config_type", True)
        self.variable_scope.add_variable("parts", "struct part", True)
        rval = " "*current_indent + "struct config_type* config = malloc(sizeof(struct config_type));\n"
        rval = rval + " "*current_indent + f"struct part* parts = {self._input_module.call_input_c(particle_count, filename)}\n"
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
        field = field.replace('"', '')

        return f"{index}->{field}"


    def get_particle_position(self, dimension):
        '''
        Returns the code to access a particle's position
        for each dimension. Dimensions are x/y/z. For C_AOS
        the positions are stored in a double[3], so we return
        the relevant array element

        :param str dimension: The dimension ("x", "y" or "z") to access

        :raises InvalidNameError: If the dimension argument is not
                                  "x", "y" or "z".
        
        :returns: The string to access a particle's position variable.
        :rtype: str
        '''
        if dimension == "x":
            return "core_part.position[0]"
        if dimension == "y":
            return "core_part.position[1]"
        if dimension == "z":
            return "core_part.position[2]"
        raise InvalidNameError("The dimension argument should be x, y, or z")

    def get_pointer(self, var_code):
        '''
        Returns the code to access the pointer to the supplied var_code.
        The var_code should already be C_AOS code.

        :param str var_code: The C_AOS code to take pointer from.

        :returns: The string pointer to the supplied var_code.
        :rtype: str
        '''
        return "&" + var_code

    def create_variable(self, c_type, name, initial_value=None, **kwargs):
        current_indent = kwargs.get("current_indent", 0)
        if C_AOS._type_map.get(c_type) is None:
            raise UnsupportedTypeError("C_AOS does not support type \"{0}\""
                                        " in created variables.".format(c_type))
        ##Check name is allowed in C
        name = name.replace('"', '')
        a = re.match("[a-zA-Z_][a-zA-Z_0-9]*", name)
        if a is None or a.group(0) != name:
            raise InvalidNameError("C_AOS does not support \"{0}\" as a name"
                                   " for variables.".format(name))
        end = ";\n"
        if initial_value is not None:
            end = f" = {initial_value};\n"
        rval = " " * current_indent + C_AOS._type_map.get(c_type) + " " + name + end
        self.variable_scope.add_variable(name, c_type, False)
        return rval

    def call_language_function(self,func_call, *args, **kwargs):
        '''
        Calls a language function from the visitors.

        :param str func_call: The function call to be called for this language.
        
        :returns: The generated code from the language function call.
        :rtype: str
        '''
        string = ""
        try:
            # Any arguments that were python module accesses would be
            # converted to c pointer accesses here (so use ->).
            # We need to change this into a python module access for now
            # in a temporary variable
            fn = getattr(self, func_call)
            fixed_args = []
            for arg in args:
                fixed_args.append(arg.replace("->", "."))
            string = fn(*fixed_args, **kwargs)
            return string
        except (AttributeError) as err:
            pass
        for system in self._coupled_systems:
            try:
                fn = getattr(system, func_call)
                fixed_args = []
                for arg in args:
                    fixed_args.append(arg.replace("->", "."))
                string = fn(*fixed_args, **kwargs)
                return string
            except (AttributeError) as err:
                pass

        string = func_call + "( "
        arguments = []
        for arg in args:
            arguments.append(arg)
        for kwarg in kwargs:
            arguments.append(str(kwarg) + "=" + str(kwargs[kwarg]))
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

    def set_cutoff(self, cutoff, var_type=CONSTANT, current_indent=0, **kwargs):
        '''
        Set the cutoff radius for pairwise interactions

        :param cutoff: The variable name or value of the cutoff radius
        :type cutoff: str or int
        :param var_type: The type of variable to use for the cutoff radius.
                         If C_AOS.CONSTANT, the value in cutoff is some constant
                         value, used for all particles. This is the default.
                         If C_AOS.PARTICLE is used, the variable in `cutoff`
                         should be a member of the particle structure, and
                         will be used individually for each particle.
        :param int current_indent: The current indent level of the code.

        :raises UnsupportedTypeError: If an unsupported var_type is supplied.

        :returns: A string containing the code to set this value if C_AOS.CONSTANT
                  is used, otherwise "".
        :returns: str
        '''
        if var_type == C_AOS.CONSTANT:
            self._cutoff_type = C_AOS.CONSTANT
            self._cutoff = "config->neighbour_config.cutoff"
            string = " "*current_indent + f"config.neighbour_config.cutoff = {cutoff};\n"
            return string
        if var_type == C_AOS.PARTICLE:
            self._cutoff_type = C_AOS.PARTICLE
            self._cutoff = f"parts[part1].{cutoff}"
            return ""
        raise UnsupportedTypeError("Unsupported var_type used in set_cutoff")

    def access_to_string(self, var_access, check_valid=False):
        '''
        Takes a variable_access and converts it to a C_AOS string to output

        :param var_access: The variable access to output as a string
        :type var_access: variable_access or str

        :raises UnsupportedTypeError: If the variable is of a type not supported
                                      by C_AOS and not a child of a structure.

        :returns: The C_AOS string for this variable access
        :rtype: str
        '''
        code_str = ""
        name = var_access.variable.var_name
        code_str = code_str + name
        array_access = (len(var_access.array_indices) != 0)
        check = check_valid and self._enable_variable_checks
        # Check for type existing
        if not var_access.is_child and check:
            if C_AOS._type_map.get(var_access.variable.var_type) is None:
                raise UnsupportedTypeError("Accessing a variable of type "
                                          f"{var_access.variable.var_type} "
                                           "which is not supported by C_AOS backend."
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
    
    def per_particle_loop_start(self, index_name):
        '''
        Returns the code to start a per_particle loop.

        :param str index_name: The name to use to index the loop

        :returns: The code to start a per_particle loop.
        :rtype: str
        '''
        loop = f"for( int {index_name} = 0; {index_name} < config->space.nparts"
        loop = loop + f"; {index_name}++)" + "{\n"
        return loop

    def particle_access(self, index_name, element):
        '''
        Returns the code to access a particle for a given index and element

        :param str index_name: The name of the index to access the particle
        :param str element: The particle element to access

        :returns: The code to access a particle
        :rtype: str
        '''
        access = f"parts[{index_name}].{element}" 
        return access

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
