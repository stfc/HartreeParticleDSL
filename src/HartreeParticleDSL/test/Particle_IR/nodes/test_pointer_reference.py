import pytest

from HartreeParticleDSL.Particle_IR.nodes.pointer_reference import PointerReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, PointerType
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol

def test_pr():
    structure = PointerSymbol("structure1", PointerType(INT_TYPE))

    with pytest.raises(TypeError) as excinfo:
        PointerReference("x")
    assert ("Attempted to make a PointerReference to a "
            "non-PointerSymbol. Got <class 'str'> as input."
            in str(excinfo.value))

    a = PointerReference(structure)
    assert a.symbol is structure

    correct = "PointerReference[name:'structure1']"
    assert correct == a.node_str()
