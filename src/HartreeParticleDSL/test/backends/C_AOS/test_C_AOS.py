from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import *
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config
from HartreeParticleDSL.backends.C_AOS.C_AOS import *
from HartreeParticleDSL.c_types import *
from HartreeParticleDSL.language_utils.variable_scope import *
import HartreeParticleDSL.kernel_types.kernels as kernels
import pytest
import os
import inspect
import ast

def test_set_io_modules():
    '''Test the set_io_modules function of C_AOS module'''
    mod = Random_Particles()
    backend = C_AOS()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "C_AOS backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "C_AOS backend does not support temp_module" in str(excinfo.value)

def test_println():
    '''Tests the println function of C_AOS module'''
    backend = C_AOS()
    statement = backend.println("%f string 123", "total", current_indent=0)
    assert statement == "printf(\"%f string 123\\n\", total);\n"

def test_addinclude():
    '''Tests the addinclude function of C_AOS module'''
    backend = C_AOS()
    backend.add_include("<string.h>")
    assert "<string.h>" in backend._includes
    backend.add_include("<string.h>")
    assert backend._includes.count("<string.h>") is 1

def test_generateincludes():
    '''Tests the generate_includes function of C_AOS module'''
    backend = C_AOS()
    mod = Random_Particles()
    backend.set_io_modules(mod, None)
    strin = backend.generate_includes()
    assert "<stdio.h>" in strin
    assert "random_io.h" in strin
    assert "<math.h>" in strin
    assert "\"part.h\"" in strin
    assert "<stdlib.h>" in strin
    backend = C_AOS()
    backend.set_io_modules(None, mod)
    strin = backend.generate_includes()
    assert "<stdio.h>" in strin
    assert "random_io.h" in strin
    assert "<math.h>" in strin
    assert "\"part.h\"" in strin
    assert "<stdlib.h>" in strin

def test_gen_headers():
    '''Test the gen_headers function of C_AOS module'''
    backend = C_AOS()
    part = Particle()
    config = Config()
    mod = Random_Particles()
    backend.set_io_modules(mod, mod)
    backend.gen_headers(config, part)
    f_str = ""
    with open('part.h', 'r') as f:
        f_str = f.readlines()
    assert f_str[0] == '#ifndef PART_H\n'
    assert f_str[1] == '#define PART_H\n'
    assert f_str[2] == 'struct space_type{\n'
    assert f_str[3] == '    double box_dims[3];\n'
    assert f_str[4] == '    int nparts;\n'
    assert f_str[5] == '};\n'
    assert f_str[6] == '\n'
    assert f_str[7] == 'struct neighbour_config_type{\n'
    assert f_str[8] == '    double cutoff;\n'
    assert f_str[9] == '};\n'
    assert f_str[10] == '\n'
    assert f_str[11] == 'struct config_type{\n'
    assert f_str[12] == '    struct space_type space;\n'
    assert f_str[13] == '    struct neighbour_config_type neighbour_config;\n'
    assert f_str[14] == '};\n'
    assert f_str[15] == '\n'
    assert f_str[16] == 'struct core_part_type{\n'
    assert f_str[17] == '    double position[3];\n'
    assert f_str[18] == '    double velocity[3];\n'
    assert f_str[19] == '};\n'
    assert f_str[20] == '\n'
    assert f_str[21] == 'struct neighbour_part_type{\n'
    assert f_str[22] == '};\n'
    assert f_str[23] == '\n'
    assert f_str[24] == 'struct part{\n'
    assert f_str[25] == '    struct core_part_type core_part;\n'
    assert f_str[26] == '    struct neighbour_part_type neighbour_part;\n'
    assert f_str[27] == '};\n'
    assert f_str[28] == '\n'
    assert f_str[29] == '#endif'
    os.remove("part.h")

def kern(part1, part2, r2, config):
    part1.a = part1.a + 2.0

def kern2(part, r2):
    part.a = part.a + 2.0

def test_gen_pairwise_kernel(capsys):
    '''Test the gen_kernel function of C_AOS module for a pairwise kernel'''
    backend = C_AOS()
    kernel = kernels.pairwise_interaction(kern)
    backend.gen_kernel(kernel)
    captured = capsys.readouterr()
    correct = 'void kern( struct part *part1, struct part *part2, double r2, struct config_type *config )\n'
    correct = correct + "{\n"
    correct = correct + "    part1->a = ( part1->a + 2.0 );\n"
    correct = correct + "}\n"
    assert correct in captured.out

def test_gen_perpart_kernel(capsys):
    '''Test the gen_kernel function of C_AOS module for a perpart kernel'''
    backend = C_AOS()
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    captured = capsys.readouterr()
    correct = 'void kern2( struct part *part, struct config_type *r2 )\n'
    correct = correct + "{\n"
    correct = correct + "    part->a = ( part->a + 2.0 );\n"
    correct = correct + "}\n"
    assert correct in captured.out

def main():
    a = a + 1

def test_print_main(capsys):
    '''Test the print_main function of C_AOS'''
    backend = C_AOS()
    m = ast.parse(inspect.getsource(main))
    backend.print_main(m)
    captured = capsys.readouterr()
    correct = 'int main(  )\n'
    correct = correct + "{\n"
    correct = correct + "    a = ( a + 1 );\n"
    correct = correct + "}\n"
    assert correct in captured.out

def test_gen_invoke_perpart():
    '''Test the gen_invoke function of C_AOS for perpart kernel'''
    backend = C_AOS()
    out = backend.gen_invoke("kern2", 0, 1, kernels.perpart_kernel_wrapper)
    correct = "\n /* INVOKE generated for kern2 */\n"
    correct = correct + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
    correct = correct + " kern2(&parts[part1], config);\n"
    correct = correct + "}\n"
    correct = correct + "/* End of INVOKE generated for kern2 */\n\n"
    assert correct in out

def test_gen_invoke_pairwise():
    '''Test the gen_invoke function of C_AOS for pairwise kernel'''
    backend = C_AOS()
    out = backend.gen_invoke("kern", 0, 1, kernels.pairwise_kernel_wrapper)
    correct = "\n /* INVOKE generated for kern */\n"
    correct = correct + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
    correct = correct + " for( int part2 = 0; part2 < config->space.nparts; part2++){\n"
    correct = correct + "  if(part1 == part2) continue;\n"
    correct = correct + "  double r2 = 0.0;\n"
    correct = correct + "  for(int k = 0; k < 3; k++){\n"
    correct = correct + "   r2 += (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]) * (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]);\n"
    correct = correct + "  }\n"
    correct = correct + "  if(r2 < (config->neighbour_config.cutoff * config->neighbour_config.cutoff)){\n"
    correct = correct + "   kern(&parts[part1], &parts[part2], r2, config);\n"
    correct = correct + "  }\n"
    correct = correct + " }\n"
    correct = correct + "}\n"
    correct = correct + "/* End of INVOKE generated for kern */"
    assert correct in out

def test_gen_invoke_fail(capsys):
    '''Test the gen_invoke function outputs nothing for unsupported kernel types'''
    pass

def test_initialisation_code(capsys):
    '''Test the initialisation_code function of C_AOS'''
    mod = Random_Particles()
    backend = C_AOS()
    backend.set_io_modules(mod, mod)
    val = backend.initialisation_code(1000, None)
    assert val == "random_io(1000, config);"

def test_gen_particle():
    '''Test the gen_particle function of C_AOS'''
    backend = C_AOS()
    part = Particle()
    part.add_element("thing", "double[4]")
    rval = backend.gen_particle(part)
    correct = "struct core_part_type{\n"
    correct = correct + "    double position[3];\n"
    correct = correct + "    double velocity[3];\n"
    correct = correct + "};\n\n"
    correct = correct + "struct neighbour_part_type{\n"
    correct = correct + "};\n\n"
    correct = correct + "struct part{\n"
    correct = correct + "    struct core_part_type core_part;\n"
    correct = correct + "    struct neighbour_part_type neighbour_part;\n"
    correct = correct + "    double thing[4];\n"
    correct = correct + "};\n\n"
    assert correct == rval

def test_gen_config():
    '''Test the gen_config function of C_AOS'''
    backend = C_AOS()
    conf = Config()
    conf.add_element("thing", "double[4]")
    rval = backend.gen_config(conf)
    correct = "struct space_type{\n"
    correct = correct + "    double box_dims[3];\n"
    correct = correct + "    int nparts;\n"
    correct = correct + "};\n\n"
    correct = correct + "struct neighbour_config_type{\n"
    correct = correct + "    double cutoff;\n"
    correct = correct + "};\n\n"
    correct = correct + "struct config_type{\n"
    correct = correct + "    struct space_type space;\n"
    correct = correct + "    struct neighbour_config_type neighbour_config;\n"
    correct = correct + "    double thing[4];\n"
    correct = correct + "};\n\n"
    assert correct == rval

def test_cleanup():
    '''Test the cleanup function of C_AOS'''
    backend = C_AOS()
    correct = "  free(config);\n  free(parts);\n"
    rval = backend.cleanup(current_indent=2)
    assert correct == rval

def test_initialise():
    '''Test the initialise function of C_AOS'''
    mod = Random_Particles()
    backend = C_AOS()
    backend.set_io_modules(mod, mod)

    correct = " struct config_type* config = malloc(sizeof(struct config_type));\n"
    correct = correct + " struct part* parts = {0}\n".format(backend._input_module.call_input_c(100, "abc.def"))
    rval = backend.initialise(particle_count=100, filename="abc.def", current_indent=1)
    assert correct in rval

def test_call_language_function():
    '''Test the test_call_language function of C_AOS'''
    backend = C_AOS()
    rval1 = backend.call_language_function("a_c_call", "*part", "20")
    assert rval1 == ("a_c_call( *part, 20 )" + ";\n")
    rval2 = backend.call_language_function("cleanup", current_indent=2, indent=1)
    assert rval2 == "  free(config);\n  free(parts);\n"
    rval3 = backend.call_language_function("a_c_call", "*part", "20", current_indent=4, indent=1)
    assert rval3 == "    a_c_call( *part, 20 );\n"
    rval4 = backend.call_language_function("set_cutoff", "0.5", "C_AOS->CONSTANT")
    assert rval4 == "config.neighbour_config.cutoff = 0.5;\n"


def test_get_particle_access():
    '''Test the get_particle_access function of C_AOS'''
    backend = C_AOS()
    a = "core_part.position[0]"
    b = "core_part.position[1]"
    c = "core_part.position[2]"
    assert backend.get_particle_access("x") == a
    assert backend.get_particle_access("y") == b
    assert backend.get_particle_access("z") == c
    with pytest.raises(InvalidNameError) as excinfo:
        backend.get_particle_access("abc")
    assert "The dimension argument should be x, y, or z" in str(excinfo.value) 


def test_create_variable():
    '''Test the create_variable function of C_AOS'''
    backend = C_AOS()

    # Test int
    out = backend.create_variable(c_int, "a")
    assert out == "int a;\n"

    # Test double with initial value
    out = backend.create_variable(c_double, "b", 0.0)
    assert out == "double b = 0.0;\n"

    # Test float with current_indent
    out = backend.create_variable(c_float, "c", current_indent=4)
    assert out == "    float c;\n"

    # Test int64
    out = backend.create_variable(c_int64_t, "d")
    assert out == "long long int d;\n"

    # Test int32
    out = backend.create_variable(c_int32_t, "e")
    assert out == "int e;\n"

    # Test int8
    out = backend.create_variable(c_int8_t, "foo")
    assert out == "char foo;\n"

    #Test bool
    out = backend.create_variable(c_bool, "g")
    assert out == "_Bool g;\n"

    #Test some illegal names
    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "space variable")
    assert ("C_AOS does not support \"space variable\" as a name"
           " for variables.") in str(excinfo.value)

    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "0abc")
    assert ("C_AOS does not support \"0abc\" as a name"
           " for variables.") in str(excinfo.value)

    #Test illegal type
    with pytest.raises(UnsupportedTypeError) as excinfo:
        out = backend.create_variable("not a type", "a")
    assert ("C_AOS does not support type \"not a type\" in"
            " created variables.") in str(excinfo.value)

def test_set_cutoff_constant():
    '''Test the set_cutoff member of the C_AOS backend with a constant
    cutoff'''
    backend = C_AOS()
    rval = backend.set_cutoff(3.5, var_type=C_AOS.CONSTANT)
    assert rval == "config.neighbour_config.cutoff = 3.5;\n"
    assert backend._cutoff_type == C_AOS.CONSTANT
    assert backend._cutoff == "config->neighbour_config.cutoff"
    out = backend.gen_invoke("kern", 0, 1, kernels.pairwise_kernel_wrapper)
    correct = "\n /* INVOKE generated for kern */\n"
    correct = correct + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
    correct = correct + " for( int part2 = 0; part2 < config->space.nparts; part2++){\n"
    correct = correct + "  if(part1 == part2) continue;\n"
    correct = correct + "  double r2 = 0.0;\n"
    correct = correct + "  for(int k = 0; k < 3; k++){\n"
    correct = correct + "   r2 += (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]) * (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]);\n"
    correct = correct + "  }\n"
    correct = correct + "  if(r2 < (config->neighbour_config.cutoff * config->neighbour_config.cutoff)){\n"
    correct = correct + "   kern(&parts[part1], &parts[part2], r2, config);\n"
    correct = correct + "  }\n"
    correct = correct + " }\n"
    correct = correct + "}\n"
    correct = correct + "/* End of INVOKE generated for kern */"
    assert correct in out

def test_set_cutoff_particle():
    '''Test the set_cutoff member of the C_AOS backend with a particle
    variable'''
    backend = C_AOS()
    rval = backend.set_cutoff("h", var_type=C_AOS.PARTICLE, current_indent=4)
    assert rval == ""
    assert backend._cutoff_type == C_AOS.PARTICLE
    assert backend._cutoff == "parts[part1].h"
    out = backend.gen_invoke("kern", 0, 1, kernels.pairwise_kernel_wrapper)
    correct = "\n /* INVOKE generated for kern */\n"
    correct = correct + "for( int part1 = 0; part1 < config->space.nparts; part1++){\n"
    correct = correct + " for( int part2 = 0; part2 < config->space.nparts; part2++){\n"
    correct = correct + "  if(part1 == part2) continue;\n"
    correct = correct + "  double r2 = 0.0;\n"
    correct = correct + "  for(int k = 0; k < 3; k++){\n"
    correct = correct + "   r2 += (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]) * (parts[part1].core_part.position[k] - parts[part2].core_part.position[k]);\n"
    correct = correct + "  }\n"
    correct = correct + "  if(r2 < (parts[part1].h * parts[part1].h)){\n"
    correct = correct + "   kern(&parts[part1], &parts[part2], r2, config);\n"
    correct = correct + "  }\n"
    correct = correct + " }\n"
    correct = correct + "}\n"
    correct = correct + "/* End of INVOKE generated for kern */"
    assert correct in out

def test_set_cutoff_error():
    '''Test the set_cutoff member of the C_AOS backend with an invalid request'''
    backend = C_AOS()
    with pytest.raises(UnsupportedTypeError) as excinfo:
        rval = backend.set_cutoff("fail", var_type="thing1")
    assert "Unsupported var_type used in set_cutoff" in str(excinfo.value)

def test_get_pointer():
    '''Test the get_pointer method of the C_AOS backend'''
    backend = C_AOS()
    strin = backend.get_pointer("variable")
    assert strin == "&variable"

def test_access_to_string():
    '''Test the access_to_string method of the C_AOS backend'''
    backend = C_AOS()
    HartreeParticleDSL.set_backend(backend)
    backend.variable_scope.add_variable("var1", "double", True)
    backend.variable_scope.add_variable("var2", "double", False)

    var1_access = variable_access(backend.variable_scope.get_variable("var1"))
    var2_access = variable_access(backend.variable_scope.get_variable("var2"))
    var1_access.child = var2_access
    var2_access.add_array_index("i")

    z = str(var1_access)
    assert z == "var1->var2[i]"

    var1_access.add_array_index("z")
    z = str(var1_access)
    assert z == "var1[z].var2[i]"
