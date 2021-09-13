from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin

class Random_Particles(IO_Module, C_AOS_IO_Mixin):
    ''' Implementation of the randomly generated particles IO module '''
    def __init__(self):
        super().__init__()

    def gen_code_c(self, part_type):
        '''
        Returns the C code required for this IO module.
        This IO module uses a prewritten C implementation so has no C output.
        '''
        return ""

    def call_input_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for input.
        This solely calls the random_io function defined in the
        required header file for this module.

        :returns: The function call used for this module.
        :rtype: str
        '''
        return f"random_io({part_count}, config);"

    def call_output_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for output.
        In this case output is not supported so an empty string is 
        returned.

        :returns: Empty string
        :rtype: str
        '''

        return ""

    def get_includes_c(self):
        '''
        :returns: The includes required for this IO module.
        rtype: List of str
        '''
        includes = []
        includes.append("<stdlib.h>")
        includes.append("\"random_io.h\"")
        return includes
