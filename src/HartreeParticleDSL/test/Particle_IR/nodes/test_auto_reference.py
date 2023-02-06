import pytest

from HartreeParticleDSL.Particle_IR.nodes.auto_reference import AutoReference
from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol

def test_auto_reference():
    sym = AutoSymbol("x", "0")

    ref = AutoReference(sym)

    correct = "AutoReference[x]"
    assert correct == ref.node_str()
    assert ref.symbol is sym

    assert ref.is_array == False
