import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, UnaryOperation

def test_binaryop_create_badargs():
    with pytest.raises(TypeError) as excinfo:
        a = BinaryOperation.create(123, [])
    assert ("BinaryOperation operator argument must be of type BinaryOperation.Operator but found int." in str(excinfo.value))

    with pytest.raises(IRGenerationError) as excinfo:
        a = BinaryOperation.create(BinaryOperation.BinaryOp.ADDITION, [])
    assert ("Attempting to create a BinaryOperation with wrong number of "
            "children. Was provided 0 but can only accept 2." 
            in str(excinfo.value))

    with pytest.raises(IRGenerationError) as excinfo:
        a = BinaryOperation.create(BinaryOperation.BinaryOp.ADDITION,["123", "123"])
    assert ("Attempting to create a BinaryOperation but first provided child "
            "is <class 'str'> instead of a DataNode." in str(excinfo.value))

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    with pytest.raises(IRGenerationError) as excinfo:
        a = BinaryOperation.create(BinaryOperation.BinaryOp.ADDITION,[ref1, "123"])
    assert ("Attempting to create a BinaryOperation but second provided child "
            "is <class 'str'> instead of a DataNode." in str(excinfo.value))


def test_binaryop_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    assert BinaryOperation._validate_child(0, "123") is False
    assert BinaryOperation._validate_child(0, ref1) is True
    assert BinaryOperation._validate_child(1, "123") is False
    assert BinaryOperation._validate_child(1, ref1) is True
    assert BinaryOperation._validate_child(2, ref1) is False

def test_binaryop_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.BinaryOp.ADDITION,[ref1, ref2])

    assert a.operator is BinaryOperation.BinaryOp.ADDITION
    assert a.children[0] is ref1
    assert a.children[1] is ref2

def test_binaryop_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.BinaryOp.ADDITION,[ref1, ref2])

    assert ("BinaryOperation[BinaryOp.ADDITION: (ScalarReference[x], "
            "ScalarReference[x2])]" in a.node_str())

def test_unaryop_create_badargs():
    with pytest.raises(TypeError) as excinfo:
        a = UnaryOperation.create(123, [])
    assert("UnaryOperation operator argument must be of type UnaryOperation.Operator but found int." in str(excinfo.value))

    with pytest.raises(IRGenerationError) as excinfo:
        a = UnaryOperation.create(UnaryOperation.UnaryOp.UNARYSUB, 212)
    assert ("Attempting to create a UnaryOperation but provided child is "
            "<class 'int'> instead of a DataNode." in str(excinfo.value))

def test_unaryop_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    assert UnaryOperation._validate_child(0, "123") is False
    assert UnaryOperation._validate_child(0, ref1) is True
    assert UnaryOperation._validate_child(1, ref1) is False

def test_unaryop_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    a = UnaryOperation.create(UnaryOperation.UnaryOp.UNARYSUB, ref1)
    assert a.operator is UnaryOperation.UnaryOp.UNARYSUB
    assert a.children[0] is ref1

def test_unaryop_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    a = UnaryOperation.create(UnaryOperation.UnaryOp.UNARYSUB, ref1)

    assert "UnaryOperation[UnaryOp.UNARYSUB: ScalarReference[x]]" in a.node_str()
