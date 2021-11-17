from abc import ABCMeta

class Backend(metaclass=ABCMeta):
    CONSTANT = 0
    def __init__(self):
        pass

    @property
    def variable_scope(self):
        return None

    def disable_variable_checks(self):
        pass

    def enable_variable_checks(self):
        pass

    def set_io_modules(self, input_module, output_module):
        pass

    def add_include(self, include_string):
        pass

    def generate_includes(self):
        pass

    def gen_headers(self, config, parts):
        pass

    def gen_config(self, config):
        pass

    def gen_particle(self, parts):
        pass

    def println(self, string, *args):
        pass

    def gen_kernel(self, kernel):
        pass

    def print_main(self, function):
        pass

    def gen_invoke(self, kernel):
        pass

    def get_particle_access(self, dimension):
        '''
        Returns the code to access a particle's position
        for each dimension. Dimensions are x/y/z
        '''
        pass

    def get_pointer(self, code):
        pass

    def initialisation_code(self, particle_count, filename):
        pass

    def create_variable(self, c_type, name, initial_value=None, **kwargs):
        pass

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        pass

    def access_to_string(self, variable_access, check_valid):
        pass
