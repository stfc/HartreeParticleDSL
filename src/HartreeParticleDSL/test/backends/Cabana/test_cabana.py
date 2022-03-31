from HartreeParticleDSL.backends.Cabana_backend.Cabana_visitors import *
from HartreeParticleDSL.backends.Cabana_backend.Cabana import *
from HartreeParticleDSL.IO_modules.base_IO_module.IO_module import IO_Module
from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import *
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config, _HartreeParticleDSL
from HartreeParticleDSL.c_types import *
from HartreeParticleDSL.language_utils.variable_scope import *
from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import *
import HartreeParticleDSL.kernel_types.kernels as kernels
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
from HartreeParticleDSL.HartreeParticleDSL import set_backend

def test_set_io_modules():
    '''Test the set_io_modules function of Cabana module'''
    mod = Random_Particles()
    backend = Cabana()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "Cabana backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "Cabana backend does not support temp_module" in str(excinfo.value)

def test_println():
    '''Tests the println function of Cabana module'''
    backend = Cabana()
    statement = backend.println("\"string 123\"", "total", current_indent=0)
    assert statement == "std::cout << \"string 123\" << total << \"\\n\";\n"

def test_addinclude():
    '''Tests the addinclude function of Cabana module'''
    backend = Cabana()
    backend.add_include("<string.h>")
    assert "<string.h>" in backend._includes
    backend.add_include("<string.h>")
    assert backend._includes.count("<string.h>") is 1

def test_generateincludes():
    '''Tests the generate_includes function of Cabana module'''
    backend = Cabana()
    mod = Random_Particles()
    backend.set_io_modules(mod, None)
    strin = backend.generate_includes()
    assert "<Cabana_Core.hpp>" in strin
    assert "random_io_cabana.hpp" in strin
    assert "<Cabana_Core.hpp>" in strin
    assert "\"part.h\"" in strin
    assert "<Kokkos_Random.hpp>" in strin
    assert "<iostream>" in strin
    assert "<cmath>" in strin
    assert "<cstdio>" in strin
    backend = Cabana()
    backend.set_io_modules(None, mod)
    strin = backend.generate_includes()
    assert "<Cabana_Core.hpp>" in strin
    assert "random_io_cabana.hpp" in strin
    assert "<Cabana_Core.hpp>" in strin
    assert "\"part.h\"" in strin
    assert "<Kokkos_Random.hpp>" in strin
    assert "<iostream>" in strin
    assert "<cmath>" in strin
    assert "<cstdio>" in strin

def test_gen_headers(capsys):
    # TODO
    pass

def kern2(part, config):
    create_variable(c_double, a, 2.0)
    part.a = part.a + 2.0

def test_gen_perpart_kernel(capsys):
    '''Test the gen_kernel function of Cabana module for a perpart kernel'''
    _HartreeParticleDSL.the_instance = None
    backend = Cabana()
    HartreeParticleDSL.set_backend(backend)
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    captured = capsys.readouterr()
    correct = '''template < class A >
struct kern2_functor{
    config_struct_type _config;
    A _a;

    KOKKOS_INLINE_FUNCTION
     kern2_functor( A a, config_struct_type config) :
    _a(a), _config(config) {}

    void operator()(const int i, const int a) const{
        double a = 2.0;
        _a.access(i, a) = ( _a.access(i, a) + 2.0 );

    }
};

'''
    assert correct in captured.out

def main():
    create_variable(c_int, a, 0)
    a = a + 1

def test_print_main(capsys):
    '''Test the print_main function of Cabana'''
    backend = Cabana()
    HartreeParticleDSL.set_backend(backend)
    m = ast.parse(inspect.getsource(main))
    backend.print_main(m)
    captured = capsys.readouterr()
    correct = '''int main( int argc, char* argv[] )
{
    int a = 0;
    a = ( a + 1 );
}

'''
    assert correct in captured.out

def test_gen_invoke_perpart():
    '''Test the gen_invoke function of Cabana for perpart kernel'''
    backend = Cabana()
    out = backend.gen_invoke("kern2", 0, 1, kernels.perpart_kernel_wrapper)
    correct = ''' /* INVOKE generated for kern2 */
Kokkos::deep_copy(config.config, config.config_host);
Cabana::simd_parallel_for(simd_policy, kern2, "kern2");
Kokkos::fence();
/* End of INVOKE generated for kern2 */
'''
    assert correct in out

def test_initialisation_code(capsys):
    '''Test the initialisation_code function of Cabana'''
    mod = Random_Particles()
    backend = Cabana()
    backend.set_io_modules(mod, mod)
    val = backend.initialisation_code(1000, None)
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 1000);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 1000);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert val == correct

def test_gen_particle():
    '''Test the gen_particle function of Cabana'''
    backend = Cabana()
    part = Particle()
    part.add_element("thing", "double[4]")
    rval = backend.gen_particle(part)
    correct = '''struct core_part_type{
    double position[3];
    double velocity[3];
};

struct neighbour_part_type{
    double cutoff;
};

enum FieldNames{thing = 0,
                 core_part_velocity,
                 core_part_position,
                 neighbour_part_cutoff
               };
using DataTypes = Cabana::MemberTypes<double[4],
    double[3],
    double[3],
    double>;
'''
    assert correct == rval

def test_gen_config():
    '''Test the gen_config function of Cabana'''
    backend = Cabana()
    conf = Config()
    conf.add_element("thing", "double[4]")
    rval = backend.gen_config(conf)
    correct = '''struct boundary{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
};

struct space_type{
    boundary box_dims;
    int nparts;
};

struct neighbour_config_type{
};

struct config_view_type{
        struct space_type space;
        struct neighbour_config_type neighbour_config;
        double thing[4];
};

using config_struct_type = Kokkos::View<struct config_view_type*, MemorySpace>;
using config_struct_host = config_struct_type::HostMirror;
struct config_type{
    config_struct_type config;
    config_struct_host config_host;
};
'''
    assert correct == rval

def test_cleanup():
    '''Test the cleanup function of Cabana'''
    backend = Cabana()
    correct = "}\n"
    rval = backend.cleanup(current_indent=2)
    assert correct == rval

def test_call_language_function():
    '''Test the test_call_language function of Cabana'''
    backend = Cabana()
    rval1 = backend.call_language_function("a_c_call", "*part", "20")
    assert rval1 == ("a_c_call( *part, 20 )" + ";\n")
    rval2 = backend.call_language_function("cleanup", current_indent=2, indent=1)
    assert rval2 == "}\n"
    rval3 = backend.call_language_function("a_c_call", "*part", "20", current_indent=4, indent=1)
    assert rval3 == "    a_c_call( *part, 20 );\n"
    with pytest.raises(NotImplementedError) as excinfo:
        val4 = backend.call_language_function("set_cutoff", "0.5", "C_AOS.CONSTANT")
    assert "Cabana backend doesn't yet support pairwise interactions" in str(excinfo.value)

class coupler_test(base_coupler):
    def __init__(self):
        pass

    def a_function(self):
        return "test_string"

    def b_function(self, arg):
        return arg

    def get_includes(self):
        return ["\"test.h\""]

def test_add_coupler():
    '''Test the add_coupler function of Cabana'''
    backend = Cabana()
    coupler = coupler_test()
    backend.add_coupler(coupler)
    assert coupler in backend._coupled_systems
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_coupler(32)
    assert ("Can only couple to base_coupler classes or "
            "subclasses. Found int") in str(excinfo.value)

def test_call_language_function_coupled_system():
    '''Test the coupled system functionality'''
    backend = Cabana()
    coupler = coupler_test()
    backend.add_coupler(coupler)
    assert "\"test.h\"" in backend._includes
    rval = backend.call_language_function("a_function")
    assert rval == "test_string"
    rval = backend.call_language_function("b_function", "thing->thing2")
    assert rval == "thing->thing2"
    rval = backend.call_language_function("unknown_func")
    assert rval == "unknown_func(  );\n"

def test_get_particle_position():
    '''Test the get_particle_position function of Cabana'''
    backend = Cabana()
    a = "core_part_position[0]"
    b = "core_part_position[1]"
    c = "core_part_position[2]"
    assert backend.get_particle_position("x") == a
    assert backend.get_particle_position("y") == b
    assert backend.get_particle_position("z") == c
    with pytest.raises(InvalidNameError) as excinfo:
        backend.get_particle_position("abc")
    assert "The dimension argument should be x, y, or z" in str(excinfo.value)

def test_create_variable():
    '''Test the create_variable function of Cabana'''
    backend = Cabana()

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
    assert out == "int64_t d;\n"

    # Test int32
    out = backend.create_variable(c_int32_t, "e")
    assert out == "int32_t e;\n"

    # Test int8
    out = backend.create_variable(c_int8_t, "foo")
    assert out == "int8_t foo;\n"

    #Test bool
    out = backend.create_variable(c_bool, "g")
    assert out == "bool g;\n"

    #Test some illegal names
    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "space variable")
    assert ("Cabana does not support \"space variable\" as a name"
           " for variables.") in str(excinfo.value)

    with pytest.raises(InvalidNameError) as excinfo:
        out = backend.create_variable(c_int, "0abc")
    assert ("Cabana does not support \"0abc\" as a name"
           " for variables.") in str(excinfo.value)

    #Test illegal type
    with pytest.raises(UnsupportedTypeError) as excinfo:
        out = backend.create_variable("not a type", "a")
    assert ("Cabana does not support type \"not a type\" in"
            " created variables.") in str(excinfo.value)

def test_get_pointer():
    '''Test the get_pointer method of the Cabana backend'''
    backend = Cabana()
    strin = backend.get_pointer("variable")
    assert strin == "&(variable)"

def test_access_to_string():
    '''Test the access_to_string method of the C_AOS backend'''
    backend = Cabana()
    HartreeParticleDSL.set_backend(backend)
    backend.variable_scope.add_variable("var1", "c_double", True)
    backend.variable_scope.add_variable("var2", "c_double", False)

    var1_access = variable_access(backend.variable_scope.get_variable("var1"))
    var2_access = variable_access(backend.variable_scope.get_variable("var2"))
    var1_access.child = var2_access
    var2_access.add_array_index("i")

    z = backend.access_to_string(var1_access, True)
    assert z == "var1->var2[i]"

    var1_access.add_array_index("z")
    z = str(var1_access)
    assert z == "var1[z].var2[i]"

    backend.variable_scope.add_variable("var3", "UNKNOWN", False)
    var3_access = variable_access(backend.variable_scope.get_variable("var3"))
    with pytest.raises(UnsupportedTypeError) as excinfo:
       z = backend.access_to_string(var3_access, True)
    assert ("Accessing a variable of type UNKNOWN which is not supported "
            "by Cabana backend.") in str(excinfo.value)

def test_add_type():
    '''Test the add_type method of the Cabana backend'''
    backend = Cabana()
    backend.add_type('test_string', 'test_type')
    assert Cabana._type_map['test_string'] == 'test_type'

def test_get_particle_access():
    '''Test the get_particle_access method of the Cabana backend'''
    backend = Cabana()
    rval = backend.get_particle_access("part1", "field")
    assert rval == "_field.access(i, a)"
    rval = backend.get_particle_access("part1", "\"field\"")
    assert rval == "_field.access(i, a)"

def test_per_particle_loop_start():
    backend = Cabana()
    rval = backend.per_particle_loop_start("i")
    assert rval is None

def test_write_output():
    backend = Cabana()
    class temp_module(IO_Module, Cabana_IO_Mixin):
        def __init__(self):
            pass
        def call_output_cabana(self, num_parts, filename, variable):
            return "Success"
    a = temp_module()
    mod = Random_Particles()
    backend.set_io_modules(mod, a)
    assert backend.write_output("a") == "Success\n"
