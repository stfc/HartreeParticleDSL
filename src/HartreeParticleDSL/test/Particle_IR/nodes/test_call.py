import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

def test_call_children():
    mycall = Call("mycall")

    with pytest.raises(IRGenerationError) as excinfo:
        mycall.addchild("123")
    assert ("Item '<class 'str'>' can't be child 0 of '<class "
            "'HartreeParticleDSL.Particle_IR.nodes.call.Call'>'." 
            in str(excinfo.value))

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    mycall.addchild(ref1)
    assert mycall.children[0] is ref1
    assert mycall.func_name == "mycall"


def test_call_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("y", INT_TYPE)
    ref2 = ScalarReference(sym2)
    mycall = Call.create("mycall", [ref1, ref2])

    assert mycall.func_name == "mycall"
    assert mycall.children[0] is ref1
    assert mycall.children[1] is ref2


def test_call_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("y", INT_TYPE)
    ref2 = ScalarReference(sym2)
    mycall = Call.create("mycall", [ref1, ref2])

    correct = "Call[mycall: (ScalarReference[x], ScalarReference[y])]"
    assert correct == mycall.node_str()
