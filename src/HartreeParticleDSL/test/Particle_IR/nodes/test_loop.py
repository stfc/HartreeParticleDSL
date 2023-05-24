import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.nodes.loop import Loop
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

def test_loop_init():
    with pytest.raises(TypeError) as excinfo:
       x = Loop("not symbol")
    assert ("Loop variable needs to be a symbol but got <class 'str'>." 
            in str(excinfo.value))
    sym = ScalarTypeSymbol("x", INT_TYPE)
    x = Loop(sym)
    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert "Loop node should have 4 children but found 0" in str(excinfo.value)

def test_loop_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref = ScalarReference(sym)

    assert Loop._validate_child(0, ref) == True
    assert Loop._validate_child(1, ref) == True
    assert Loop._validate_child(2, ref) == True
    assert Loop._validate_child(3, ref) == False

    bdy = Body()
    assert Loop._validate_child(3, bdy) == True

def test_loop_create():
    with pytest.raises(IRGenerationError) as excinfo:
        Loop.create("a", "b", "c", "d", [])
    assert ("Loop variable needs to be a symbol but got <class 'str'>." 
            in str(excinfo.value))

    sym = ScalarTypeSymbol("i", INT_TYPE)

    with pytest.raises(IRGenerationError) as excinfo:
        Loop.create(sym, "b", "c", "d", [])
    assert ("Loop start needs to be a DataNode but got <class 'str'>."
            in str(excinfo.value))

    ref1 = ScalarReference(ScalarTypeSymbol("start", INT_TYPE))
    with pytest.raises(IRGenerationError) as excinfo:
        Loop.create(sym, ref1, "c", "d", [])
    assert ("Loop stop needs to be a DataNode but got <class 'str'>."
            in str(excinfo.value))

    ref2 = ScalarReference(ScalarTypeSymbol("stop", INT_TYPE))
    with pytest.raises(IRGenerationError) as excinfo:
        Loop.create(sym, ref1, ref2, "d", [])
    assert ("Loop step needs to be a DataNode but got <class 'str'>."
            in str(excinfo.value))

    ref3 = ScalarReference(ScalarTypeSymbol("step", INT_TYPE))

    call = Call("mycall")

    x = Loop.create(sym, ref1, ref2, ref3, [call])

    assert x.variable is sym
    assert x.start_expr is ref1
    assert x.stop_expr is ref2
    assert x.step_expr is ref3
    assert isinstance(x.body, Body)
    assert x.body.children[0] is call

def test_loop_node_str():
    sym = ScalarTypeSymbol("i", INT_TYPE)
    ref1 = ScalarReference(ScalarTypeSymbol("start", INT_TYPE))
    ref2 = ScalarReference(ScalarTypeSymbol("stop", INT_TYPE))
    ref3 = ScalarReference(ScalarTypeSymbol("step", INT_TYPE))
    call = Call("mycall")

    x = Loop.create(sym, ref1, ref2, ref3, [call])
    correct = '''Loop[ (ScalarReference[name:'start'], ScalarReference[name:'stop'] ScalarReference[name:'step']): Body[
    Call[mycall: ()]
] End Body
]'''
    assert correct == x.node_str()
