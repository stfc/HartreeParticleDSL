from HartreeParticleDSL.IO_modules.random_IO.random_IO import *

def test_Random_Particles_gencodec():
    mod = Random_Particles()
    assert mod.gen_code_c(None) is ""

def test_Random_Particles_call_input_c():
    mod = Random_Particles()
    assert mod.call_input_c(1000, "filename") == "random_io(1000, config);"

def test_Random_Particles_call_output_c():
    mod = Random_Particles()
    assert mod.call_output_c(0, None) is ""

def test_getincludes_c():
    mod = Random_Particles()
    includes = mod.get_includes_c()
    assert "<stdlib.h>" in includes
    assert "\"random_io.h\"" in includes

def test_getincludes_fdps():
    mod = Random_Particles()
    includes = mod.get_includes_fdps()
    assert "\"random_io.hpp\"" in includes

def test_getincludes_cabana():
    mod = Random_Particles()
    includes = mod.get_includes_cabana()
    assert "\"random_io_cabana.hpp\"" in includes

def test_gencode_fdps():
    mod = Random_Particles()
    assert mod.gen_code_fdps(None) == ""

def test_call_input_fdps():
    mod = Random_Particles()
    x = mod.call_input_fdps(123, "test.file")
    correct = '''PS::ParticleSystem<FullParticle> particle_system;
    particle_system.initialize();
    particle_system.setNumberOfParticleLocal(123);
    random_io( particle_system, config);
'''
    assert x == correct

def test_gencode_cabana():
    mod = Random_Particles()
    assert mod.gen_code_cabana(None) == ""

def test_call_input_cabana():
    mod = Random_Particles()
    x = mod.call_input_cabana(part_count=123, filename="blah")
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", 123);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", 123);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert x == correct

    x = mod.call_input_cabana(part_count="part_count=123", filename="blah")
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", 123);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", 123);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert x == correct

def test_input_cabana_pir():
    mod = Random_Particles()
    x = mod.call_input_cabana_pir(part_count="part_count=123", filename="blah")
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( \"particle_list\", 123);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( \"particle_list_host\", 123);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert x == correct

def test_gen_code_cabana_pir():
    mod = Random_Particles()
    x = mod.gen_code_cabana_pir(None)
    assert x == ""

def test_get_includes_cabana_pir():
    mod = Random_Particles()
    x = mod.get_includes_cabana_pir()
    assert "\"random_io_cabana.hpp\"" in x

def test_get_header_includes_and_box_size_cabana_pir():
    mod = Random_Particles()
    assert mod.call_get_box_size_pir(None, None) == ""
    assert len(mod.get_header_includes_cabana_pir()) == 0

