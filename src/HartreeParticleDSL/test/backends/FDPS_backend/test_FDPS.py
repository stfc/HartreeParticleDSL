from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import *
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config, _HartreeParticleDSL
from HartreeParticleDSL.backends.FDPS_backend.FDPS import *
from HartreeParticleDSL.c_types import *
import HartreeParticleDSL.kernel_types.kernels as kernels
import pytest
import os
import inspect
import ast

def test_set_io_modules():
    ''' Tests the set_io_modules function of the FDPS backend'''
    mod = Random_Particles()
    backend = FDPS()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "FDPS backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "FDPS backend does not support temp_module" in str(excinfo.value)

def test_println():
    ''' Tests the println module of the FDPS backend'''
    backend = FDPS()
    statement = backend.println("\"\"", "total", "\"\"string123\"\"", current_indent=0)
    assert statement == "std::cout << \"\" << total << \"string123\" << \"\\n\";\n"

def test_addinclude():
    ''' Tests the addinclude module of the FDPS backend'''
    backend = FDPS()
    backend.add_include("<string.h>")
    assert "<string.h>" in backend._includes
    backend.add_include("<string.h>")
    assert backend._includes.count("<string.h>") is 1

def test_generate_include():
    ''' Tests the generate_include module of the FDPS backend'''
    backend = FDPS()
    mod = Random_Particles()
    backend.set_io_modules(mod, None)
    strin = backend.generate_includes()
    assert "<cmath>" in strin
    assert "<cstdio>" in strin
    assert "<iostream>" in strin
    assert "<vector>" in strin
    assert "<particle_simulator.hpp>" in strin
    assert "\"part.h\"" in strin
    backend = FDPS()
    backend.set_io_modules(None, mod)
    string = backend.generate_includes()
    assert "<cmath>" in strin
    assert "<cstdio>" in strin
    assert "<iostream>" in strin
    assert "<vector>" in strin
    assert "<particle_simulator.hpp>" in strin
    assert "\"part.h\"" in strin

def test_gen_headers():
    ''' Tests the gen_headers module of the FDPS backend'''
    backend = FDPS()
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
    assert f_str[2] == 'struct boundary{\n'
    assert f_str[3] == '    PS::F64 x_min, x_max;\n'
    assert f_str[4] == '    PS::F64 y_min, y_max;\n'
    assert f_str[5] == '    PS::F64 z_min, z_max;\n'
    assert f_str[6] == '};\n'
    assert f_str[7] == '\n'
    assert f_str[8] == 'struct space_type{\n'
    assert f_str[9] == '    boundary box_dims;\n'
    assert f_str[10] == '    PS::S32 nparts;\n'
    assert f_str[11] == '};\n'
    assert f_str[12] == '\n'
    assert f_str[13] == 'struct neighbour_config_type{\n'
    assert f_str[14] == '};\n'
    assert f_str[15] == '\n'
    assert f_str[16] == 'class config_type{\n'
    assert f_str[17] == '    public:\n'
    assert f_str[18] == '        struct space_type space;\n'
    assert f_str[19] == '        struct neighbour_config_type neighbour_config;\n'
    assert f_str[20] == '};\n'
    assert f_str[21] == '\n'
    assert f_str[22] == 'struct core_part_type{\n'
    assert f_str[23] == '    PS::F64vec position;\n'
    assert f_str[24] == '    PS::F64 velocity[3];\n'
    assert f_str[25] == '};\n'
    assert f_str[26] == '\n'
    assert f_str[27] == 'struct neighbour_part_type{\n'
    assert f_str[28] == '    PS::F64 cutoff;\n'
    assert f_str[29] == '};\n'
    assert f_str[30] == '\n'
    assert f_str[31] == 'class FullParticle{\n'
    assert f_str[32] == '    public:\n'
    assert f_str[33] == '        struct core_part_type core_part;\n'
    assert f_str[34] == '        struct neighbour_part_type neighbour_part;\n'
    assert f_str[35] == '        PS::F64vec getPos(){\n'
    assert f_str[36] == '            return this->core_part.position;\n'
    assert f_str[37] == '        }\n'
    assert f_str[38] == '        void setPos(const PS::F64vec pos_new){\n'
    assert f_str[39] == '            this->core_part.position = pos_new;\n'
    assert f_str[40] == '        }\n'
    assert f_str[41] == '        void clear(){\n'
    assert f_str[42] == '        }\n'
    assert f_str[43] == '};\n'
    assert f_str[44] == '\n'
    assert f_str[45] == '#endif'
    os.remove("part.h")

def kern(part1, part2, r2, config):
    part1.a = part1.a + 2.0

def kern2(part, r2):
    part.a = part.a + 2.0

def test_gen_pairwise_kernel(capsys):
    '''Test the gen_kernel function of FDPS module for a pairwise kernel'''
    _HartreeParticleDSL.the_instance = None
    backend = FDPS()
    kernel = kernels.pairwise_interaction(kern)
    with pytest.raises(UnsupportedCodeError) as excinfo:
        backend.gen_kernel(kernel)
    assert "FDPS doesn't yet support pairwise operations" in str(excinfo.value)

def test_gen_perpart_kernel(capsys):
    '''Test the gen_kernel function of FDPS module for a perpart kernel'''
    backend = FDPS()
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    captured = capsys.readouterr()
    correct = 'void kern2( FullParticle& part, config_type& r2 )\n'
    correct = correct + "{\n"
    correct = correct + "    part.a = ( part.a + 2.0 );\n"
    correct = correct + "}\n"
    assert correct in captured.out

def main():
    a = a + 1

def test_print_main(capsys):
    '''Test the print_main function of FDPS'''
    backend = FDPS()
    m = ast.parse(inspect.getsource(main))
    backend.print_main(m)
    captured = capsys.readouterr()
    correct = 'int main(  )\n'
    correct = correct + "{\n"
    correct = correct + "    a = ( a + 1 );\n"
    correct = correct + "}\n"
    assert correct in captured.out

def test_gen_invoke_perpart():
    '''Test the gen_invoke function of FDPS for perpart kernel'''
    backend = FDPS()
    out = backend.gen_invoke("kern2", 0, 1, kernels.perpart_kernel_wrapper)
    correct = "\n /* INVOKE generated for kern2 */\n"
    correct = correct + "for(PS::S32 i = 0; i < particle_system.getNumberOfParticleLocal(); i++){\n"
    correct = correct + " kern2(particle_system[i], config);\n"
    correct = correct + "}\n"
    correct = correct + "/* End of INVOKE generated for kern2 */\n\n"
    assert correct in out

def test_gen_invoke_pairwise():
    '''Test the gen_invoke function of FDPS for pairwise kernel'''
    backend = FDPS()
    with pytest.raises(NotImplementedError) as excinfo:
        out = backend.gen_invoke("kern", 0, 1, kernels.pairwise_kernel_wrapper)
    assert "gen_invoke not yet implemented" in str(excinfo.value)

def test_initialisation_code(capsys):
    '''Test the initialisation_code function of FDPS'''
    mod = Random_Particles()
    backend = FDPS()
    backend.set_io_modules(mod, mod)
    val = backend.initialisation_code(1000, None)
    assert val == ("PS::ParticleSystem<FullParticle> particle_system;\n"
                   "    particle_system.initialize();\n"
                   "    particle_system.setNumberOfParticleLocal(1000);\n"
                   "    random_io( particle_system, config);\n")

def test_gen_particle():
    '''Test the gen_particle function of FDPS'''
    backend = FDPS()
    part = Particle()
    part.add_element("thing", "double[4]")
    rval = backend.gen_particle(part)
    correct = "struct core_part_type{\n"
    correct = correct + "    PS::F64vec position;\n"
    correct = correct + "    PS::F64 velocity[3];\n"
    correct = correct + "};\n\n"
    correct = correct + "struct neighbour_part_type{\n"
    correct = correct + "    PS::F64 cutoff;\n"
    correct = correct + "};\n\n"
    correct = correct + "class FullParticle{\n"
    correct = correct + "    public:\n"
    correct = correct + "        struct core_part_type core_part;\n"
    correct = correct + "        struct neighbour_part_type neighbour_part;\n"
    correct = correct + "        double thing[4];\n"
    correct = correct + "        PS::F64vec getPos(){\n"
    correct = correct + "            return this->core_part.position;\n"
    correct = correct + "        }\n"
    correct = correct + "        void setPos(const PS::F64vec pos_new){\n"
    correct = correct + "            this->core_part.position = pos_new;\n"
    correct = correct + "        }\n"
    correct = correct + "        void clear(){\n"
    correct = correct + "        }\n"
    correct = correct + "};\n\n"
    assert correct == rval

def test_gen_config():
    '''Test the gen_config function of FDPS'''
    backend = FDPS()
    conf = Config()
    conf.add_element("thing", "double[4]")
    rval = backend.gen_config(conf)
    correct = "struct boundary{\n"
    correct = correct + "    PS::F64 x_min, x_max;\n"
    correct = correct + "    PS::F64 y_min, y_max;\n"
    correct = correct + "    PS::F64 z_min, z_max;\n"
    correct = correct + "};\n\n"
    correct = correct + "struct space_type{\n"
    correct = correct + "    boundary box_dims;\n"
    correct = correct + "    PS::S32 nparts;\n"
    correct = correct + "};\n\n"
    correct = correct + "struct neighbour_config_type{\n"
    correct = correct + "};\n\n"
    correct = correct + "class config_type{\n"
    correct = correct + "    public:\n"
    correct = correct + "        struct space_type space;\n"
    correct = correct + "        struct neighbour_config_type neighbour_config;\n"
    correct = correct + "        double thing[4];\n"
    correct = correct + "};\n\n"
    assert correct == rval

def test_cleanup():
    '''Test the cleanup function of FDPS'''
    backend = FDPS()
    correct = "  PS::Finalize();\n"
    rval = backend.cleanup(current_indent=2)
    assert correct == rval

def test_call_language_function():
    '''Test the test_call_language function of FDPS'''
    backend = FDPS()
    func1 = "a_c_call( *part, 20 )"
    rval1 = backend.call_language_function("a_c_call", "*part", "20")
    assert rval1 == (func1 + ";\n")
    rval2 = backend.call_language_function("cleanup", current_indent=2, indent=1)
    assert rval2 == "  PS::Finalize();\n"
    func3 = "a_c_call(*part, 20, current_indent=4, indent=1)"
    rval3 = backend.call_language_function("a_c_call", "*part", "20", current_indent=4, indent=1)
    assert rval3 == "    a_c_call( *part, 20 );\n"

def test_create_variable():
    '''Test the create_variable function of FDPS'''
    backend = FDPS()

    # Test int
    out = backend.create_variable(c_int, "a")
    assert out == "PS::S32 a;\n"

    # Test double with initial value
    out = backend.create_variable(c_double, "b", 0.0)
    assert out == "PS::F64 b = 0.0;\n"

    # Test float with current_indent
    out = backend.create_variable(c_float, "c", current_indent=4)
    assert out == "    PS::F32 c;\n"

    # Test int64
    out = backend.create_variable(c_int64_t, "d")
    assert out == "PS::S64 d;\n"

    # Test int32
    out = backend.create_variable(c_int32_t, "e")
    assert out == "PS::S32 e;\n"

    # Test int8
    out = backend.create_variable(c_int8_t, "foo")
    assert out == "char foo;\n"

    #Test bool
    out = backend.create_variable(c_bool, "g")
    assert out == "bool g;\n"

    #Test some illegal names
    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "space variable")
    assert ("FDPS does not support \"space variable\" as a name"
           " for variables.") in str(excinfo.value)

    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "0abc")
    assert ("FDPS does not support \"0abc\" as a name"
           " for variables.") in str(excinfo.value)

    #Test illegal type
    with pytest.raises(UnsupportedTypeError) as excinfo:
        out = backend.create_variable("not a type", "a")
    assert ("FDPS does not support type \"not a type\" in"
            " created variables.") in str(excinfo.value)

def test_set_cutoff():
    '''Tests the set_cutoff function of FDPS'''
    backend = FDPS()
    with pytest.raises(NotImplementedError) as excinfo:
        backend.set_cutoff(None, None)
    assert "FDPS backend doesn't yet support pairwise interactions" in \
            str(excinfo.value)

def test_initialise():
    '''Test the initialise function of FDPS'''
    mod = Random_Particles()
    backend = FDPS()
    backend.set_io_modules(mod, mod)

    correct = " char **argv = NULL;\n"
    correct = correct + " int args = 0;\n"
    correct = correct + " PS::Initialize(args,argv);\n"
    correct = correct + " config_type config;\n"
    correct = correct + " {0}\n".format(backend._input_module.call_input_fdps(100, "abc.def", current_indent=1))
    rval = backend.initialise(particle_count=100, filename="abc.def", current_indent=1)
    assert correct in rval

def test_get_particle_access():
    '''Test the get_particle_access function of FDPS'''
    backend = FDPS()
    a = "core_part.position.x"
    b = "core_part.position.y"
    c = "core_part.position.z"
    assert backend.get_particle_access("x") == a
    assert backend.get_particle_access("y") == b
    assert backend.get_particle_access("z") == c
    with pytest.raises(InvalidNameError) as excinfo:
        backend.get_particle_access("abc")
    assert "The dimension argument should be x, y, or z" in str(excinfo.value)

def test_get_pointer():
    '''Test the get_pointer method of the FDPS backend'''
    backend = FDPS()
    strin = backend.get_pointer("variable")
    assert strin == "&variable"
