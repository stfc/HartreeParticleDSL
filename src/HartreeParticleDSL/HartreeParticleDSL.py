import os
from HartreeParticleDSL.backends.base_backend import backend
import HartreeParticleDSL.IO_modules.base_IO_module.IO_module as io_modules
import HartreeParticleDSL.IO_modules.random_IO.random_IO as random_io
from HartreeParticleDSL.HartreeParticleDSLExceptions import SingletonInstanceError, \
                                                            RepeatedNameError, \
                                                            NoBackendError

class _HartreeParticleDSL():
    '''
    Base class used for HartreeParticleDSL.

    This class is hidden from the user and other modules and accessed through
    a variety of utility functions available to the user and other modules.
    Only one instance of this class can ever exist, and this class should 
    always be accessed through the get_instance() function.

    :param Backend: The backend to use for this application. This can be\
                    changed later with the set_backend function. Default\
                    value is the C_AOS backend.
    :type Backend: :py:class:`HartreeParticleDSL.backends.base_backend.Backend`
    '''
    the_instance = None

    def __init__(self, Backend=backend.Backend()):
        if _HartreeParticleDSL.the_instance is not None:
            raise SingletonInstanceError("Only one instance of _HartreeParticleDSL"
                                         " is allowed")
        self._input_module = None
        self._output_module = None
        self._part_type = None
        self._config_type = None
        self._kernel_names = []
        self._kernels = []
        self._backend = Backend
        self._includes = []
        self._includes.append("<math.h>")
        self._includes.append("<stdio.h>")
        self._includes.append("\"part.h\"")
        self._STARTDIR = os.getcwd()
        self._outdir = "."
        _HartreeParticleDSL.the_instance = self

    def get_instance():
        '''
        Helper function used to get the instance of _HartreeParticleDSL
        used for this program.

        :returns: Instance of _HartreeParticleDSL
        :rtype: _HartreeParticleDSL
        '''
        if _HartreeParticleDSL.the_instance is None:
            _HartreeParticleDSL()
        return _HartreeParticleDSL.the_instance

    def set_backend(self, backend):
        '''
        Function that sets the backend used by HartreeParticleDSL

        :param backend: The backend to use for this application.
        :type backend: :py:class:`HartreeParticleDSL.backends.base_backend.Backend`
        '''
        self._backend = backend
        self._backend.set_io_modules(self._input_module, self._output_module)

    def get_backend(self):
        '''
        Returns the backend used by HartreeParticleDSL

        :returns: The current backend
        :type backend: :py:class:`HartreeParticleDSL.backends.base_backend.Backend`
        '''
        return self._backend

    def set_particle_type(self, part):
        '''
        Function to set the Particle type used by HartreeParticleDSL

        :param part: The Particle to use for this application.
        :type part: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Particle`
        '''
        self._part_type = part

    def set_config_type(self, config):
        '''
        Function to set the Config type used by HartreeParticleDSL

        :param part: The Config to use for this application.
        :type part: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Config`
        '''
        self._config_type = config

    def set_io_modules(self, input_module=None, output_module=None):
        '''
        Function to set the Input and Output modules used by HartreeParticleDSL.
        Both the input and output modules can be the same module, but may also
        be different modules.

        :param input_module: the IO module to use for input for this 
                             application.
        :type input_module: :py:class:`HartreeParticleDSL.IO_Module.base_IO_module.IO_module`
        :param output_module: the IO module to use for output for this 
                              application.
        :type input_module: :py:class:`HartreeParticleDSL.IO_Module.base_IO_module.IO_module`
        '''
        self._input_module = input_module
        self._output_module = output_module
        if self._backend is not None:
            self._backend.set_io_modules(self._input_module, self._output_module)

    def register_kernel(self, kernel_name, kernel):
        '''
        Function to register a kernel used in HartreeParticleDSL. This is
        usually accessed directly by the wrapper function and not needed for
        user-level code.

        :param str kernel_name: The name of the kernel.
        :param kernel: The Kernel object corresponding to the kernel_name.
        :type kernel: :py:class:`HartreeParticleDSL.kernel_types.kernels.kernel`
        '''
        assert kernel_name not in self._kernel_names
        self._kernel_names.append(kernel_name)
        self._kernels.append(kernel)

    # Don't support casting yet so adding this for now
    def gen_random_double(self):
        print("double random_double(){")
        print("    return (double)(rand()) / (double)(RAND_MAX);")
        print("}")

    def generate_code(self):
        '''
        Function that generates the first section of code for the chosen 
        kernels, particle and config.

        This calls the _backend and generates the include code, the required
        header files and the kernel code

        At the moment this is just output to standard out.
        '''
        # Add the io module information to the backend chosen just in case
        self._backend.set_io_modules(self._input_module, self._output_module)

        #Get the backend to generate the includes
        print(self._backend.generate_includes())
        #Dummy placeholder for a generic function I'm using for example
        self.gen_random_double()
        #Get backends to generate the headers required (e.g. particle and config types)
        self._backend.gen_headers(self._config_type, self._part_type)
        for kernel in self._kernels:
            self._backend.gen_kernel(kernel)

    def initialisation_code(self, particle_count, filename):
        '''
        Function that returns the initialisation code call for the chosen
        IO module and backend combination
        
        :returns: The initialisation call for the chosen configuration.
        :rtype: str
        '''
        return self._backend.initialisation_code(particle_count=particle_count, filename=filename)

    def generate_invoke(self, kernel_name, current_indent, indent):
        '''
        Generates the code used to invoke a specific kernel. This is accessed
        by visitors in the main function to create code for a specific invoke
        call.

        At current this is output to stdout.

        :param kernel_name: The name of the kernel to create invoke code for.
        :type kernel_name: str
        :param current_indent: The current indentation level in the code.
        :type current_indent: int
        :param indent: The amount to indent for each line of code.
        :type indent: int

        :returns: The code for this invocation
        :rval: str
        '''
        assert kernel_name in self._kernel_names
        return self._backend.gen_invoke(kernel_name, current_indent, indent, type(self._kernels[self._kernel_names.index(kernel_name)]) )

    def println(self, string, *args):
        '''
        Generates the code to perform a println call. This is accessed by
        visitors in the main function to create code for a specific println
        call.

        Work in progress.

        :param string: The base string to use in a println call.
        :type string: str
        :param args: Any number of strings used for formatted string generation
        :type args: [str]+
        '''
        return self._backend.println(string, *args)

    def print_main(self, function):
        '''
        Generates the main function code for the specified function and backend.

        This currently outputs for stdout.

        :param function: The main function to be generated by the backend.
        :type function: :py:class:`
        '''
        self._backend.print_main(function)

    def set_output_dir(self, directory):
        '''
        Sets the output directory for HartreeParticleDSL.
        Relative paths are always used from the location the DSL is run from.

        :param str directory: The directory to use for output by HartreeParticleDSL.
        '''
        # return to initial directory
        os.chdir(self._STARTDIR)
        # Save this directory
        self._outdir = directory
        # Make the directory if it doesn't exist
        if not os.path.exists(directory):
            os.mkdir(directory)
        os.chdir(directory)

def set_particle_type(part):
    '''
    Function to set the Particle type used by HartreeParticleDSL

    :param part: The Particle to use for this application.
    :type part: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Particle`
    '''
    _HartreeParticleDSL.get_instance().set_particle_type(part)

def set_config_type(config):
    '''
    Function to set the Config type used by HartreeParticleDSL

    :param part: The Config to use for this application.
    :type part: :py:class:`HartreeParticleDSL.HartreeParticleDSL.Config`
    '''
    _HartreeParticleDSL.get_instance().set_config_type(config)

def register_kernel(kernel_name, kernel):
    '''
    Function to register a kernel used in HartreeParticleDSL. This is
    usually accessed directly by the wrapper function and not needed for
    user-level code.

    :param kernel_name: The name of the kernel.
    :type kernel_name: str
    :param kernel: The Kernel object corresponding to the kernel_name.
    :type kernel: :py:class:`HartreeParticleDSL.kernel_types.kernels.kernel
    '''
    _HartreeParticleDSL.get_instance().register_kernel(kernel_name, kernel)

def set_backend(backend):
    '''
    Function that sets the backend used by HartreeParticleDSL

    :param backend: The backend to use for this application.
    :type backend: :py:class:`HartreeParticleDSL.backends.base_backend.Backend`
    '''
    _HartreeParticleDSL.get_instance().set_backend(backend)

def get_backend():
    '''
    Function to retrieve the currently used backend.

    :returns: The current backend
    :rtype: :py:class:`HartreeParticleDSL.backends.base_backend.Backend`
    '''
    return _HartreeParticleDSL.get_instance().get_backend()

def set_io_modules(input_mod, output_mod):
    '''
    Function to set the Input and Output modules used by HartreeParticleDSL.
    Both the input and output modules can be the same module, but may also
    be different modules.

    :param input_module: the IO module to use for input for this 
                         application.
    :type input_module: :py:class:`HartreeParticleDSL.IO_Module.base_IO_module.IO_module`
    :param output_module: the IO module to use for output for this 
                          application.
    :type input_module: :py:class:`HartreeParticleDSL.IO_Module.base_IO_module.IO_module`
    '''
    _HartreeParticleDSL.get_instance().set_io_modules(input_mod, output_mod)

def gen_code():
    '''
    Function that generates the first section of code for the chosen 
    kernels, particle and config.

    This calls the chosen backend and generates the include code, the required
    header files and the kernel code

    At the moment this is just output to stdout.
    '''
    _HartreeParticleDSL.get_instance().generate_code()

def initialise(particle_count=0, filename=""):
    '''
    Function that returns the initialisation code call for the chosen
    IO module and backend combination
    
    :returns: The initialisation call for the chosen configuration.
    :rtype: str
    '''
    return _HartreeParticleDSL.get_instance().initialisation_code(particle_count=particle_count, filename = filename)

def gen_invoke(kernel_name, current_indent, indent):
    '''
    Generates the code used to invoke a specific kernel. This is accessed
    by visitors in the main function to create code for a specific invoke
    call.

    At current this is output to stdout.

    :param kernel_name: The name of the kernel to create invoke code for.
    :type kernel_name: str
    :param current_indent: The current indentation level in the code.
    :type current_indent: int
    :param indent: The amount to indent for each line of code.
    :type indent: int

    :returns: The string to invoke this kernel
    :rval: str
    '''
    return _HartreeParticleDSL.get_instance().generate_invoke(kernel_name, current_indent, indent)

def println(string, *args):
    '''
    Generates the code to perform a println call. This is accessed by
    visitors in the main function to create code for a specific println
    call.

    Work in progress.

    :param string: The base string to use in a println call.
    :type string: str
    :param args: Any number of strings used for formatted string generation
    :type args: [str]+
    '''
    return _HartreeParticleDSL.get_instance().println(string, *args)

def print_main(function):
    '''
    Generates the main function code for the specified function and backend.

    This currently outputs for stdout.

    :param function: The main function to be generated by the backend.
    :type function: :py:class:`
    '''
    return _HartreeParticleDSL.get_instance().print_main(function)

def set_output_dir(directory):
    '''
    Sets the output dirctory for HartreeParticleDSL.

    :param str directory: The directory to use for output by HartreeParticleDSL.
    '''
    _HartreeParticleDSL.get_instance().set_output_dir(directory)

class Particle():
    '''Particle class used in HartreeParticleDSL'''
    def __init__(self):
        self.particle_type = {}
        self.particle_type['core_part'] = {'type' : 'struct core_part_type', 'is_array' : False }
        self.particle_type['neighbour_part'] = {'type' : 'struct neighbour_part_type', 'is_array' : False}

    def reset_particle(self):
        '''
        Reset the particle class to its initial state.
        '''
        self.particle_type = {}
        self.particle_type['core_part'] = {'type' : 'struct core_part_type', 'is_array' : False }
        self.particle_type['neighbour_part'] = {'type' : 'struct neighbour_part_type', 'is_array' : False}

    def add_element(self, variable_name, c_type):
        '''
        Add an element to the Particle type.
        
        :param variable_name: The name of the variable type. This must be unique.
        :type variable_name: str
        :param c_type: The c_type associated with the variable type. At current
                       this is a string but may change in future.
        :type c_type: str
        '''
        if variable_name in self.particle_type:
            raise RepeatedNameError(f"The variable name {variable_name} is already in the particle type")
        is_array = False
        if "[" in c_type:
            is_array = True
        self.particle_type[variable_name] = {'type' : c_type, 'is_array' : is_array}


class Config():
    '''
        Config class used in HartreeParticleDSL.
        The Config class can be used to store simulation specific data,
        including constants etc.
        Support for modifiable fields in the Config class is a work in progress.
    '''
    def __init__(self):
        self.config_type = {}
        self.config_type['space'] = {'type' : 'struct space_type', 'is_array' : False }
        self.config_type['neighbour_config'] = {'type' : 'struct neighbour_config_type' , 'is_array' : False}

    def reset_config(self):
        '''
        Reset the Config class to its initial state.
        '''
        self.config_type = {}
        self.config_type['space'] = {'type' : 'struct space_type', 'is_array' : False }
        self.config_type['neighbour_config'] = {'type' : 'struct neighbour_config_type' , 'is_array' : False}

    def add_element(self, variable_name, c_type):
        '''
        Add an element to the Config type.
        
        :param variable_name: The name of the variable type. This must be unique.
        :type variable_name: str
        :param c_type: The c_type associated with the variable type. At current
                       this is a string but may change in future.
        :type c_type: str

        :raises RepeatedNameError: If variable_name is already in the Config
        '''
        if variable_name in self.config_type:
            raise RepeatedNameError(f"The variable name {variable_name} is already in the config type")
        is_array = False
        if "[" in c_type:
            is_array = True
        self.config_type[variable_name] = {'type' : c_type, 'is_array' : is_array}
