from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin
from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin

class Random_Particles(IO_Module, C_AOS_IO_Mixin, FDPS_IO_Mixin):
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

    def gen_code_fdps(self, part_type):
        '''
        Generates and returns the FDPS code required for this IO module.

        This module uses predefined C++ code for FDPS.

        :returns: The code FDPS code required for this IO module.
        :rtype: str
        '''
        return ""

    def call_input_fdps(self, part_count, filename):
        '''
        Returns the call required to use this IO module for input.

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        # For FDPS, the IO module needs to create the object and size it up
        # before calling the function
        rval = "PS::ParticleSystem<FullParticle> particle_system;\n"
        rval = rval + "    particle_system.initialize();\n"
        rval = rval + f"    particle_system.setNumberOfParticleLocal({part_count});" + "\n"
        rval = rval + "    random_io( particle_system, config);\n"
        return rval
