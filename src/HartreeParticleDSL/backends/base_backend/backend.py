from abc import ABCMeta

class Backend(metaclass=ABCMeta):
    CONSTANT = 0
    def __init__(self):
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

    def initialisation_code(self, particle_count, filename):
        pass

    def create_variable(self, c_type, name, initial_value=None, **kwargs):
        pass

    def set_cutoff(self, cutoff, var_type=CONSTANT):
        pass
