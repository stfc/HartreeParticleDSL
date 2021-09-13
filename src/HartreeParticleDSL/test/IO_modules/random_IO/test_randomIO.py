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
