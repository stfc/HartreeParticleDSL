import pytest
import ast
import textwrap
import inspect
import os
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSLExceptions import SingletonInstanceError, \
                                                            RepeatedNameError, \
                                                            NoBackendError
from HartreeParticleDSL.backends.base_backend.backend import Backend
from HartreeParticleDSL.backends.C_AOS.C_AOS import C_AOS
from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.kernel_types.kernels import perpart_kernel_wrapper

def test_get_instance():
    '''Test get_instance function of _HartreeParticleDSL '''
    HartreeParticleDSL._HartreeParticleDSL.the_instance = None
    a = HartreeParticleDSL._HartreeParticleDSL.get_instance()
    assert isinstance(a, HartreeParticleDSL._HartreeParticleDSL)

def test_single_instance():
    '''Test only one instance of HartreeParticleDSL is allowed'''
    HartreeParticleDSL._HartreeParticleDSL.the_instance = None
    a = HartreeParticleDSL._HartreeParticleDSL()
    with pytest.raises(SingletonInstanceError) as excinfo:
        b = HartreeParticleDSL._HartreeParticleDSL()
    assert "Only one instance of _HartreeParticleDSL is allowed" in str(excinfo.value)

def test_get_backend():
    '''Test the get_backend fucntion of HartreeParticleDSL'''
    assert isinstance(HartreeParticleDSL.get_backend(), Backend)

def test_set_backend():
    '''Test the set_backend function of HartreeParticleDSL'''
    backend = C_AOS()
    HartreeParticleDSL._HartreeParticleDSL.get_instance().set_backend(backend)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._backend is backend

def test_set_particle_type():
    '''Test the set_particle_type function of HartreeParticleDSL'''
    part = HartreeParticleDSL.Particle()
    HartreeParticleDSL._HartreeParticleDSL.get_instance().set_particle_type(part)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._part_type is part

def test_set_config_type():
    '''Test the set_config_type function of HartreeParticleDSL'''
    conf = HartreeParticleDSL.Config()
    HartreeParticleDSL._HartreeParticleDSL.get_instance().set_config_type(conf)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._config_type is conf

def test_set_io_modules():
    '''Test the set_io_modules function of HartreeParticleDSL'''
    # Use None io modules for this test
    HartreeParticleDSL._HartreeParticleDSL.get_instance().set_io_modules(None, None)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._input_module is None
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._output_module is None

def test_gen_random_double(capsys):
    '''Temporary function to test gen_random_double'''
    HartreeParticleDSL._HartreeParticleDSL.get_instance().gen_random_double()
    captured = capsys.readouterr()
    assert "double random_double(){\n    return (double)(rand()) / (double)(RAND_MAX);\n}\n" == captured.out

def test_config_init():
    '''Test config initialisation'''
    conf = HartreeParticleDSL.Config()
    assert conf.config_type['space'] == {'type' : 'struct space_type', 'is_array' : False }
    assert conf.config_type['neighbour_config'] == {'type': 'struct neighbour_config_type' , 'is_array' : False}

def test_initialisation_code():
    '''Test initialisation_code function'''
    r = Random_Particles()
    HartreeParticleDSL._HartreeParticleDSL.get_instance().set_io_modules(r, None)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._input_module is not None
    a = HartreeParticleDSL._HartreeParticleDSL.get_instance().initialisation_code(100, "abc.def")
    assert a == "random_io(100, config);"

def test_set_output_dir():
    '''Test the set_output_dir function'''
    orig = os.getcwd()
    HartreeParticleDSL.set_output_dir("test_out_dir")
    assert "test_out_dir" == HartreeParticleDSL._HartreeParticleDSL.get_instance()._outdir 
    assert "test_out_dir" in os.getcwd()
    assert orig in HartreeParticleDSL._HartreeParticleDSL.get_instance()._STARTDIR
    HartreeParticleDSL.set_output_dir(".")
    assert "test_out_dir" not in os.getcwd()
    assert orig == os.getcwd()
    os.rmdir("test_out_dir")


def test_println():
    '''Test println function'''
    a = HartreeParticleDSL._HartreeParticleDSL.get_instance().println("%f", "val")
    assert a == "printf(\"%f\\n\", val);\n"

def test_register_kernel():
    '''Test register_kernel'''
    def kern(part1, config):
        s = s - 1
    kern = ast.parse(textwrap.dedent(inspect.getsource(kern)))
    parser = perpart_kernel_wrapper(kern)
    HartreeParticleDSL._HartreeParticleDSL.get_instance().register_kernel("kern", parser)
    assert "kern" in HartreeParticleDSL._HartreeParticleDSL.get_instance()._kernel_names
    assert parser in HartreeParticleDSL._HartreeParticleDSL.get_instance()._kernels

def test_generate_code(capsys):
    '''Test generate_code'''
    HartreeParticleDSL.get_backend().variable_scope.add_variable("s", "c_int", False)
    HartreeParticleDSL._HartreeParticleDSL.get_instance().generate_code()
    captured = capsys.readouterr()
    correct = "#include <math.h>\n"
    correct = correct + "#include <stdio.h>\n"
    correct = correct + "#include \"part.h\"\n"
    correct = correct + "#include <stdlib.h>\n"
    correct = correct + "#include \"random_io.h\"\n\n"
    correct = correct + "double random_double(){\n"
    correct = correct + "    return (double)(rand()) / (double)(RAND_MAX);\n"
    correct = correct + "}\n"
    correct = correct + "void kern( struct part *part1, struct config_type *config )\n"
    correct = correct + "{\n"
    correct = correct + "    s = ( s - 1 );\n"
    correct = correct + "}\n\n"
    assert correct == captured.out


def test_global_set_particle_type():
    '''Test the globally used set_particle_type function'''
    HartreeParticleDSL._HartreeParticleDSL.the_instance = None
    part = HartreeParticleDSL.Particle()
    HartreeParticleDSL.set_particle_type(part)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._part_type is part

def test_global_set_config_type():
    '''Test the globally used set_config_type function'''
    conf = HartreeParticleDSL.Config()
    HartreeParticleDSL.set_config_type(conf)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._config_type is conf

def test_global_set_backend():
    '''Test the globally used set_backend function'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._backend is aos

def test_global_set_io_modules():
    '''Test the globally used set_io_modules function'''
    HartreeParticleDSL.set_io_modules(None, None)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._input_module is None
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._output_module is None

def test_global_gen_code(capsys):
    '''Test the globally used gen_code function'''
    HartreeParticleDSL.gen_code()
    captured = capsys.readouterr()
    assert ("#include <math.h>\n#include <stdio.h>\n#include \"part.h\""
            "\n\ndouble random_double(){\n"
            "    return (double)(rand()) / (double)(RAND_MAX);\n}\n") == captured.out

def test_global_initalise():
    '''Test the globally used initalise function'''
    r = Random_Particles()
    HartreeParticleDSL.set_io_modules(r, None)
    assert HartreeParticleDSL._HartreeParticleDSL.get_instance()._input_module is r
    a = HartreeParticleDSL.initialise(100, "abc.def")
    assert a == "random_io(100, config);"

def test_global_println():
    '''Test the globally used println function'''
    HartreeParticleDSL.set_backend(C_AOS())
    a = HartreeParticleDSL.println("%f", "val")
    assert a == "printf(\"%f\\n\", val);\n"

def test_global_print_main(capsys):
    def main():
        create_variable(c_int, z)
        z = z + 1
    m_func = ast.parse(textwrap.dedent(inspect.getsource(main)))
    HartreeParticleDSL.print_main(m_func)
    captured = capsys.readouterr()
    assert "int main(  )\n{\n    int z;\n    z = ( z + 1 );\n}\n\n" == captured.out


def test_config_add_element():
    '''Test the add_element function of the config'''
    conf = HartreeParticleDSL.Config()
    conf.add_element("value", "double")
    assert conf.config_type["value"] == {'type' : 'double', 'is_array': False}
    conf.add_element("value2", "double[4]")
    assert conf.config_type["value2"] == {'type' : 'double[4]', 'is_array' : True}
    with pytest.raises(RepeatedNameError) as excinfo:
        conf.add_element("value", "int")
    assert "The variable name value is already in the config type" in str(excinfo.value)

def test_config_reset():
    conf = HartreeParticleDSL.Config()
    conf.add_element("value", "double")
    conf.reset_config()
    assert "value" not in conf.config_type.keys()
    assert conf.config_type['space'] == {'type' : 'struct space_type', 'is_array' : False }
    assert conf.config_type['neighbour_config'] == {'type': 'struct neighbour_config_type' , 'is_array' : False}

def test_particle_init():
    '''Test particle initialisation'''
    part = HartreeParticleDSL.Particle()
    part.add_element("value", "double")
    part.reset_particle()
    assert part.particle_type['core_part'] == {'type' : 'struct core_part_type', 'is_array' : False }
    assert part.particle_type['neighbour_part'] == {'type' : 'struct neighbour_part_type', 'is_array' : False}

def test_particle_add_element():
    '''Test the add_element function of the particle'''
    part = HartreeParticleDSL.Particle()
    part.add_element("value", "double")
    assert part.particle_type["value"] == {'type' : 'double', 'is_array': False}
    part.add_element("value2", "double[4]")
    assert part.particle_type["value2"] == {'type' : 'double[4]', 'is_array' : True}
    with pytest.raises(RepeatedNameError) as excinfo:
        part.add_element("value", "int")
    assert "The variable name value is already in the particle type" in str(excinfo.value)

