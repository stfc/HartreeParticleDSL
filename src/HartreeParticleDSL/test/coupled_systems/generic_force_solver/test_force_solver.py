from HartreeParticleDSL.coupled_systems.generic_force_solver.force_solver import *

def test_init():
    a = force_solver()
    assert len(a.get_includes()) == 0
    assert len(a.get_includes_header()) == 0
