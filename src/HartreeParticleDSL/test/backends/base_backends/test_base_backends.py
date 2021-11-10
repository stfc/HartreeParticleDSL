from HartreeParticleDSL.backends.base_backend.backend import *
from HartreeParticleDSL.backends.base_backend.visitors import *


def test_base_backend():
    '''Test the Backend abstract class for coverage'''
    class test_backend(Backend):
        pass
    a = test_backend()
    a.set_io_modules(None, None)
    a.add_include("")
    a.generate_includes()
    a.gen_headers(None, None)
    a.gen_config(None)
    a.gen_particle(None)
    a.println(None)
    a.gen_kernel(None)
    a.print_main(None)
    a.gen_invoke(None)
    a.initialisation_code(None, None)
    a.create_variable(None, None, None)
    a.set_cutoff(None)
    a.get_particle_access(None)
    a.get_pointer(None)

def test_base_visitors():
    a = baseVisitor(None)
    a = main_baseVisitor(None)
    a = pairwise_baseVisitor(None)
    a = per_part_baseVisitor(None)
