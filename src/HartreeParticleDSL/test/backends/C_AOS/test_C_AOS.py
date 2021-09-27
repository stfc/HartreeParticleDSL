from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config
from HartreeParticleDSL.backends.C_AOS.C_AOS import *
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
    print(f_str)
    assert f_str[0] == '#ifndef PART_H\n'
    assert f_str[1] == '#define PART_H\n'
    assert f_str[2] == 'struct space_type{\n'
    assert f_str[3] == '    double box_dims[3];\n'
    assert f_str[4] == '    int nparts;\n', '};\n'
    assert f_str[5] == '};\n'
    assert f_str[6] == '\n'
    assert f_str[7] == 'struct neighbour_config_type{\n'
    assert f_str[8] == '};\n'
    assert f_str[9] == '\n'
    assert f_str[10] == 'struct config_type{\n'
    assert f_str[11] == '    struct space_type space;\n'
    assert f_str[12] == '    struct neighbour_config_type neighbour_config;\n'
    assert f_str[13] == '};\n'
    assert f_str[14] == '\n'
    assert f_str[15] == 'struct core_part_type{\n'
    assert f_str[16] == '    double position[3];\n'
    assert f_str[17] == '    double velocity[3];\n'
    assert f_str[18] == '};\n'
    assert f_str[19] == '\n'
    assert f_str[20] == 'struct neighbour_part_type{\n'
    assert f_str[21] == '};\n'
    assert f_str[22] == '\n'
    assert f_str[23] == 'struct part{\n'
    assert f_str[24] == '    struct core_part_type core_part;\n'
    assert f_str[25] == '    struct neighbour_part_type neighbour_part;\n'
    assert f_str[26] == '};\n'
    assert f_str[27] == '\n'
    assert f_str[28] == '#endif'
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
    correct = correct + "  if(r2 < config->cutoff*config->cutoff){\n"
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

def test_call_language_function():
    '''Test the test_call_language function of C_AOS'''
    backend = C_AOS()
    func1 = "a_c_call(*part, 20)"
    rval1 = backend.call_language_function(func1)
    assert rval1 == func1
    func2 = "cleanup(current_indent=2, indent=1)"
    rval2 = backend.call_language_function(func2)
    assert rval2 == "  free(config);\n  free(parts);\n"
