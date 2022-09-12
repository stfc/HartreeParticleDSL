import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.while_loop import While
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

def test_while_init():
    x = While()
    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert "While node should have 2 children but found 0" in str(excinfo.value)

def test_while_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.BinaryOp.LESS_THAN,[ref1, ref2])

    assert While._validate_child(0, a) == True
    assert While._validate_child(1, a) == False

    bdy = Body()
    assert While._validate_child(1, bdy) == True

def test_while_create():
    with pytest.raises(IRGenerationError) as excinfo:
        While.create("a", "b")
    assert ("While condition needs to be a DataNode but got <class 'str'>." 
            in str(excinfo.value))

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.BinaryOp.LESS_THAN,[ref1, ref2])

    call = Call("mycall")

    x = While.create(a, [call])

    assert x.condition is a
    assert isinstance(x.body, Body)
    assert x.body.children[0] is call

def test_loop_node_str():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.BinaryOp.LESS_THAN,[ref1, ref2])
    call = Call("mycall")

    x = While.create(a, [call])
    correct = '''While[BinaryOperation[BinaryOp.LESS_THAN: (ScalarReference[x], ScalarReference[x2])]: Body[
    Call[mycall: ()]
] End Body
]'''
    assert correct == x.node_str()
