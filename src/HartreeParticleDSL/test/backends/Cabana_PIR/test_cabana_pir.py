import ast
import inspect
import pytest

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config, _HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import part, config

from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR

from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import *

from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *

from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType, INT_TYPE, type_mapping_str

from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
import HartreeParticleDSL.kernel_types.kernels as kernels

from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError, \
                                                            IRGenerationError

def test_cabana_pir_init():
    a = Cabana_PIR()
    assert "<Cabana_Core.hpp>" in a._includes
    assert "<Kokkos_Core.hpp>" in a._includes
    assert "<Kokkos_Random.hpp>" in a._includes
    assert "<iostream>" in a._includes
    assert "<cmath>" in a._includes
    assert "<cstdio>" in a._includes
    assert "\"part.h\"" in a._includes

def test_cabana_pir_register_kernel():
    a = Cabana_PIR()
    Kern1 = Kern()
    Kern2 = Kern()
    a.register_kernel("x", Kern1)
    assert a._kernels["x"] == Kern1
    with pytest.raises(InvalidNameError) as excinfo:
        a.register_kernel("x", Kern2)
    assert "Kernel with name x already exists." in str(excinfo.value)

    a.register_kernel("y", Kern2)
    assert a._kernels["y"] == Kern2

def test_cabana_pir_add_structure():
    a = Cabana_PIR()
    x = StructureType()
    a.add_structure(x, "x")
    assert a._structures["x"] == x

def test_cabana_pir_add_kernel_slice():
    a = Cabana_PIR()
    a.add_kernel_slices("mykern", "slice1")
    assert "slice1" in  a._kernel_slices["mykern"]

    a.add_kernel_slices("mykern", "slice1")
    assert len(a._kernel_slices.keys()) == 1
    assert len(a._kernel_slices["mykern"]) == 1

    a.add_kernel_slices("mykern", "slice2")
    assert len(a._kernel_slices["mykern"]) == 2

def test_cabana_pir_add_type():
    backend = Cabana_PIR()
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_type("t", 123)
    assert ("Attempting to add a new type but got <class 'int'> "
            "but expected a Particle IR DataType instead." in str(excinfo.value))
    backend.add_type("t", INT_TYPE)
    assert type_mapping_str["t"] == INT_TYPE
    with pytest.raises(InvalidNameError) as excinfo:
        backend.add_type("t", INT_TYPE)
    assert ("Attempting to add new type t but a type with that name already "
            "exists." in str(excinfo.value))
    del type_mapping_str["t"]

def test_cabana_pir_create_global_variable():
    backend = Cabana_PIR()
    backend.create_global_variable(INT_TYPE, "hello", "0")
    assert backend._global_values["hello"] == "0"
    assert HartreeParticleDSL.global_symbol_table().lookup("hello").datatype == INT_TYPE
    with pytest.raises(IRGenerationError) as excinfo:
        backend.create_global_variable(INT_TYPE, "hello", "0")
    assert "Tried to add a new symbol hello but it was already present in the symbol table." in str(excinfo.value)

def test_cabana_pir_set_io_modules():
    mod = Random_Particles()
    backend = Cabana_PIR()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)

def test_cabana_pir_add_include():
    mod = Random_Particles()
    backend = Cabana_PIR()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)

def test_cabana_pir_generate_includes():
    backend = Cabana_PIR()
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
    backend = Cabana_PIR()
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

class coupler_test2(base_coupler):
    def __init__(self):
        pass

    def a_function(self):
        return "test_string"

    def b_function(self, arg):
        return arg

    def get_includes(self):
        return ["a"]

    def get_includes_header(self):
        return ["a"]

class dummy_module(IO_Module, Cabana_IO_Mixin):

    def __init__(self):
        pass

    def gen_code_cabana(self, part):
        return "hello"


def test_cabana_pir_gen_headers(capsys):
    backend = Cabana_PIR()
    part = Particle()
    config = Config()
    coupler = coupler_test2()
    backend.add_coupler(coupler)
    mod = dummy_module()
    backend.set_io_modules(mod, mod)
    _HartreeParticleDSL.the_instance = None
    backend.create_global_variable(INT_TYPE, "hello", "0")
    backend.gen_headers(config, part)
    assert "hello\n\n\nhello" in capsys.readouterr().out
    f_str = ""
    with open('part.h', 'r') as f:
        f_str = f.readlines()
    f_str = "".join(f_str)
    correct = '''#ifndef PART_H
#define PART_H
#include <Kokkos_Core.hpp>
#include <Cabana_Core.hpp>
#include a
/*using MemorySpace = Kokkos::CudaSpace;*/
using MemorySpace = Kokkos::HostSpace;
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using DeviceType = Kokkos::Device<Kokkos::DefaultExecutionSpace, MemorySpace>;
using HostType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
const int VectorLength = 16;

int hello = 0;
struct boundary{
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
};

using config_struct_type = Kokkos::View<struct config_view_type*, MemorySpace>;
using config_struct_host = config_struct_type::HostMirror;
struct config_type{
    config_struct_type config;
    config_struct_host config_host;
};
struct core_part_type{
    double position[3];
    double velocity[3];
};

struct neighbour_part_type{
    double cutoff;
};

enum FieldNames{core_part_velocity = 0,
                 core_part_position,
                 neighbour_part_cutoff
               };
using DataTypes = Cabana::MemberTypes<double[3],
    double[3],
    double>;
#endif'''
    print(f_str)
    assert f_str == correct

def kern2(part: part, config: config):
    create_variable(c_double, a, 2.0)
    part.core_part.position[0] = part.core_part.position[0] + 2.0

def test_cabana_pir_gen_kernel():
    _HartreeParticleDSL.the_instance = None
    backend = Cabana_PIR()
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    assert backend._kernels_pir["kern2"] is not None

def main():
    create_variable(c_int, a, 0)
    a = a + 1

def test_print_main():
    '''Test the print_main function of Cabana'''
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    m = ast.parse(inspect.getsource(main))
    config = Config()
    part = Particle()
    backend.gen_headers(config, part)
    backend.print_main(m)
    with open('code.cpp') as f:
        out = f.read()
        correct = '''int main( int argc, char* argv[] ){
    struct config config;
    int a;
    a = 0;
    a = (a + 1);
}
'''
        assert correct in out

def test_cabana_pir_set_cutoff():
    a = Cabana_PIR()
    with pytest.raises(NotImplementedError) as excinfo:
        a.set_cutoff(0)
    assert ("Cabana PIR backend doesn't yet support pairwise interactions" in
            str(excinfo.value))

def test_cabana_pir_initialisation_code():
    mod = Random_Particles()
    backend = Cabana_PIR()
    backend.set_io_modules(mod, mod)
    val = backend.initialisation_code(1000, None)
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 1000);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 1000);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert val == correct

def test_cabana_pir_gen_particle():
    backend = Cabana_PIR()
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

def test_cabana_pir_gen_config():
    backend = Cabana_PIR()
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

def test_cabana_pir_cleanup():
    backend = Cabana_PIR()
    st = backend.cleanup(current_indent=0)
    assert st == "}\n"

def test_cabana_pir_println():
    backend = Cabana_PIR()
    statement = backend.println("\"string 123\"", "total", current_indent=0)
    assert statement == "std::cout << \"string 123\" << total << \"\\n\";\n"

def test_cabana_pir_initialise():
    mod = Random_Particles()
    backend = Cabana_PIR()
    part = Particle()
    config = Config()
    backend.set_io_modules(mod, mod)
    backend.gen_particle(part)
    backend.gen_config(config)
    backend.add_structure("c_int", "mystruct")
    backend._kernel_slices["slice1"] = ["s1", "s2"]
    out = backend.initialise(123, "myfile", 4)
    correct = '''    Kokkos::ScopeGuard scope_guard(argc, argv);
{
    config_type config;
    config.config = config_struct_type("config", 1);
    config.config_host = Kokkos::create_mirror_view(config.config);
    Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 123);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 123);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);

    Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());
    c_int mystruct;
    auto core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    auto neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    slice1_functor<decltype(s1_slice), decltype(s2_slice)> slice1(s1_slice, s2_slice, config.config,mystruct);'''
    print(out)
    assert correct in out

def test_cabana_pir_add_coupler():
    backedn = Cabana_PIR()
    coupler = coupler_test()
    backend.add_coupler(coupler)
    assert coupler in backend._coupled_systems
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_coupler(32)
    assert ("Can only couple to base_coupler classes or "
            "subclasses. Found int") in str(excinfo.value)

def test_cabana_pir_write_output():
    a = Cabana_PIR()
    assert False

def test_cabana_pir_call_language_function():
    backend = Cabana_PIR()
    with pytest.raises(AttributeError):
        rval1 = backend.call_language_function("a_c_call", "*part", "20")
    rval2 = backend.call_language_function("cleanup", current_indent=2, indent=1)
    assert rval2 == "}\n"
    with pytest.raises(AttributeError):
        rval3 = backend.call_language_function("a_c_call", "*part", "20", current_indent=4, indent=1)
    with pytest.raises(NotImplementedError) as excinfo:
        val4 = backend.call_language_function("set_cutoff", "0.5", "C_AOS.CONSTANT")
    assert "Cabana PIR backend doesn't yet support pairwise interactions" in str(excinfo.value)

class coupler_test(base_coupler):
    def __init__(self):
        pass

    def a_function(self):
        return "test_string"

    def b_function(self, arg):
        return arg

    def get_includes(self):
        return ["\"test.h\""]

def test_cabana_pir_add_coupler():
    backend = Cabana_PIR()
    coupler = coupler_test()
    backend.add_coupler(coupler)
    assert coupler in backend._coupled_systems
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_coupler(32)
    assert ("Can only couple to base_coupler classes or "
            "subclasses. Found int") in str(excinfo.value)

def test_cabana_pir_write_output():
    backend = Cabana_PIR()
    class temp_module(IO_Module, Cabana_IO_Mixin):
        def __init__(self):
            pass
        def call_output_cabana(self, num_parts, filename, variable):
            return "Success"
    a = temp_module()
    mod = Random_Particles()
    backend.set_io_modules(mod, a)
    assert backend.write_output("a") == "Success\n"
