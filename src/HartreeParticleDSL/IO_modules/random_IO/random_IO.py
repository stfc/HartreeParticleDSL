from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin
from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin
from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin
from HartreeParticleDSL.backends.Cabana_PIR_backend.Cabana_PIR_IO_Mixin import Cabana_PIR_IO_Mixin

class Random_Particles(IO_Module, C_AOS_IO_Mixin, FDPS_IO_Mixin, Cabana_IO_Mixin, Cabana_PIR_IO_Mixin):
    ''' Implementation of the randomly generated particles IO module.

    This module supports all backends.

    '''
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
        :rtype: List of str

        '''
        includes = []
        includes.append("<stdlib.h>")
        includes.append("\"random_io.h\"")
        return includes

    def get_includes_fdps(self):
        '''
        :returns: The includes required for this IO module for FDPS
        :rtype: List of str

        '''
        includes = []
        includes.append("\"random_io.hpp\"")
        return includes

    def gen_code_fdps(self, part_type):
        '''
        Generates and returns the FDPS code required for this IO module.

        This module uses predefined C++ code for FDPS.

        :returns: The code FDPS code required for this IO module.
        :rtype: str
        '''
        return ""

    def call_input_fdps(self, part_count, filename, current_indent=4):
        '''
        Returns the call required to use this IO module for input.

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        # For FDPS, the IO module needs to create the object and size it up
        # before calling the function
        indentation = " " * current_indent
        rval = "PS::ParticleSystem<FullParticle> particle_system;\n"
        rval = rval + indentation+ "particle_system.initialize();\n"
        rval = rval + f"{indentation}particle_system.setNumberOfParticleLocal({part_count});" + "\n"
        rval = rval + indentation + "random_io( particle_system, config);\n"
        return rval

    def gen_code_cabana(self, part_type):
        '''
        Generates and returns the Cabana code required for this IO module.

        :returns: The code Cabana code required for this IO module.
        :rtype: str
        '''
        return ""

    def call_input_cabana(self, part_count, filename, current_indent=4):
        '''
        Returns the call required to use this IO module for input.

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        if type(part_count) is not int:
            part_count = part_count[part_count.index("=")+1:]
        indentation = " " * current_indent
        rval = f"Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", {part_count});\n"
        rval = rval + f"{indentation}Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", {part_count});\n"
        rval = rval + indentation + "random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);\n"
        rval = rval + indentation + "Cabana::deep_copy(particle_aosoa, particle_aosoa_host);\n"

        return rval

    def get_includes_cabana(self):
        '''
        Returns the includes required to use this IO module for Cabana.

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"random_io_cabana.hpp\"")
        return includes

    def gen_code_cabana_pir(self, part_type):
        '''
        Generates and returns the Cabana PIR code required for this IO module.

        :returns: The Cabana PIR code required for this IO module.
        :rtype: str
        '''
        return ""

    def call_input_cabana_pir(self, part_count, filename, current_indent=4):
        '''
        Returns the call required to use this IO module for input.

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        if type(part_count) is not int:
            part_count = part_count[part_count.index("=")+1:]
        indentation = " " * current_indent
        rval = f"Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", {part_count});\n"
        rval = rval + f"{indentation}Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", {part_count});\n"
        rval = rval + indentation + "random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);\n"
        rval = rval + indentation + "Cabana::deep_copy(particle_aosoa, particle_aosoa_host);\n"

        return rval

    def call_get_box_size_pir(self, part_count, filename, current_indent=4):
        '''
        TODO
        '''
        #TODO Documentation
        return ""

    def get_includes_cabana_pir(self):
        '''
        Returns the includes required to use this IO module for Cabana PIR.

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        includes = []
        includes.append("\"random_io_cabana.hpp\"")
        return includes

    def get_header_includes_cabana_pir(self):
        return []
