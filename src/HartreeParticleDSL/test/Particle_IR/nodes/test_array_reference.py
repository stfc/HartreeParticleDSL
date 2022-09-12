import pytest

from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.array_reference import ArrayReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, ArrayType
from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol

def test_array_reference_invalid_symbol():
    with pytest.raises(TypeError) as excinfo:
        ArrayReference("notvalid", [])

    assert ("Attempted to make an ArrayReference to a non-ArraySymbol. "
            "Got <class 'str'> as input.") in str(excinfo.value)

def test_array_reference_init():
    one = Literal("1", INT_TYPE)
    two = Literal("2", INT_TYPE)

    sym = ArraySymbol("x", ArrayType(INT_TYPE, [5, 5]))

    ar = ArrayReference(sym, [one, two])

    assert ar.symbol is sym
    assert ar.indices[0] is one
    assert ar.indices[1] is two

    ar2 = ArrayReference(sym)
    assert ar2.symbol is sym
    assert len(ar2.indices) == 0

def test_array_reference_node_str():
    one = Literal("1", INT_TYPE)
    two = Literal("2", INT_TYPE)


    sym = ArraySymbol("x", ArrayType(INT_TYPE, [5, 5]))
    ar = ArrayReference(sym, [one, two])

    correct = ("ArrayReference[x: (Literal['1', Scalar<INTEGER, SINGLE>], Literal['2', "
               "Scalar<INTEGER, SINGLE>])]")
    assert correct == ar.node_str()
