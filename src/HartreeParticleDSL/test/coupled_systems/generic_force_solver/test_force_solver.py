import pytest

from HartreeParticleDSL.coupled_systems.generic_force_solver.force_solver import *

def test_init():
    a = force_solver()
    assert len(a.get_includes()) == 0
    assert len(a.get_includes_header()) == 0
    assert a.has_preferred_decomposition() == False
    with pytest.raises(NotImplementedError):
        a.get_preferred_decomposition(0, 0, 0)
    assert len(a.get_extra_symbols(0)) == 0
    with pytest.raises(NotImplementedError):
        a.copy_files()
    assert len(a.compilation_files()) == 0
    assert len(a.get_required_packages()) == 0
