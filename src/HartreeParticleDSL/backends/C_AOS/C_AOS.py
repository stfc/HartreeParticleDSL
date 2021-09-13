from HartreeParticleDSL.backends.base_backend.backend import Backend
import HartreeParticleDSL.backends.C_AOS.visitors as c_visitors
import HartreeParticleDSL.kernel_types.kernels as kernels
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin

class C_AOS(Backend):
    '''
    C_AOS Backend class.
    Outputs from this class are Serial C code with Arrays of Struct particle
    representation.
    '''
    def __init__(self):
        self._pairwise_visitor = c_visitors.c_pairwise_visitor()
        self._per_part_visitor = c_visitors.c_perpart_visitor()
        self._main_visitor = c_visitors.c_main_visitor()
        self._includes = []
        self._includes.append("<math.h>")
        self._includes.append("<stdio.h>")
        self._includes.append("\"part.h\"")
        self._input_module = None
        self._output_module = None

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

    def println(self, string, *args):
        '''
        Function to output the println string for the C_AOS module.
        Called via the visitors when reaching a println statement in user
        code.

        For the C_AOS module this is a call to `printf`, using the 
        string argument as the formatted string with a line end, and
        the *args contains any C values used in the formatted string.

        :param string: The formatted string to use with printf.
        :type string: str
        :param *args: A list of strings containing the C values required
                      for output.
        :type args: str
        '''
        output = f"printf(\"{string}\\n\""
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
        input_module_header = self._input_module.gen_code_c(part_type)
        if input_module_header is not "":
            # Do something later
            pass
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
        if isinstance(kernel, kernels.perpart_kernel_wrapper):
            self._per_part_visitor.visit(tree)
        elif isinstance(kernel, kernels.pairwise_kernel_wrapper):
            self._pairwise_visitor.visit(tree)

    def print_main(self, function):
        '''
        Generates the C_AOS code for the supplied main function.
        Currently outputs to STDOUT

        :param function: AST.Module of the Main function to generate code for.
        :type function: :py:class:`ast.Module`
        '''
        self._main_visitor.visit(function)

    def gen_invoke(self, kernel_name, current_indent, indent, kernel_type):
        space = " "
        print(f"\n {space*current_indent}/* INVOKE generated for {kernel_name} */")
        if kernel_type == kernels.perpart_kernel_wrapper:
            print(" "*current_indent, end="")
            print("for( int part1 = 0; part1 < config->space.nparts; part1++){")
            current_indent = current_indent + 4
            print(" "*current_indent, end="")
            print(f"{kernel_name}(&parts[part1], config);")
            current_indent = current_indent - 4
            print(" "*current_indent, end="")
            print("}")
        elif kernel_type == kernels.pairwise_kernel_wrapper:
            print(" "*current_indent, end="")
            print("for( int part1 = 0; part1 < config->space.nparts; part1++){")
            current_indent = current_indent + 4
            print(" "*current_indent, end="")
            print("for( int part2 = 0; part2 < config->space.nparts; part2++){")
            current_indent = current_indent + 4
            print(" "*current_indent, end="")
            print("if(part1 == part2) continue;")
            print(" "*current_indent, end="")
            print("double r2 = 0.0;")
            print(" "*current_indent, end="")
            print("for(int k = 0; k < 3; k++){")
            current_indent = current_indent + 4
            print(" "*current_indent, end="")
            print(" r2 += (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]) * (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]);")
            current_indent = current_indent - 4
            print(" "*current_indent, end="")
            print("}")
            print(" "*current_indent, end="")
            print("if(r2 < config->cutoff*config->cutoff){")
            current_indent = current_indent + 4
            print(" "*current_indent, end="")
            print(f"{kernel_name}(&parts[part1], &parts[part2], r2, config);")
            current_indent = current_indent - 4
            print(" "*current_indent, end="")
            print("}")
            current_indent = current_indent - 4
            print(" "*current_indent, end="")
            print("}")
            current_indent = current_indent - 4
            print(" "*current_indent, end="")
            print("}")
        print(f"{space*current_indent}/* End of INVOKE generated for {kernel_name} */\n")

    def initialisation_code(self, particle_count, filename):
        return self._input_module.call_input_c(particle_count, filename)

    def gen_particle(self, particle):
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
        # Output the space
        output = ""
        output = output + "struct space_type{\n"
        output = output + "    double box_dims[3];\n"
        output = output + "    int nparts;\n"
        output = output + "};\n\n"

        # Currently empty neiughbour config
        output = output + "struct neighbour_config_type{\n"
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
