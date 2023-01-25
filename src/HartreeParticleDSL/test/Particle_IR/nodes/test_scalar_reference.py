import pytest

from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

def test_scalar_reference():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    
    ref = ScalarReference(sym)

    correct = "ScalarReference[x]"
    assert correct == ref.node_str()
    assert ref.symbol is sym

    assert ref.is_array == False

def test_scalar_ref_illegalsym():
    with pytest.raises(TypeError) as excinfo:
        ScalarReference("abc")
    assert ("Attempted to make a ScalarReference to a non-ScalarTypeSymbol. "
            "Got <class 'str'> as input.") in str(excinfo.value)
