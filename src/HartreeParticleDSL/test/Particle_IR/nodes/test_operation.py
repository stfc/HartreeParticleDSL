import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from psyclone.psyir.nodes import BinaryOperation, UnaryOperation
from psyclone.errors import GenerationError

def test_binaryop_create_badargs():
    with pytest.raises(GenerationError) as excinfo:
        a = BinaryOperation.create(123, "", "")
    assert ("operator argument in create method of BinaryOperation class "
            "should be a PSyIR BinaryOperation Operator but found 'int'." in str(excinfo.value))

    with pytest.raises(GenerationError) as excinfo:
        a = BinaryOperation.create(BinaryOperation.Operator.ADD,"123", "123")
    assert ("Item 'str' can't be child 0 of 'BinaryOperation'. "
            "The valid format is: 'DataNode, DataNode'." in str(excinfo.value))

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    with pytest.raises(GenerationError) as excinfo:
        a = BinaryOperation.create(BinaryOperation.Operator.ADD,ref1, "123")
    assert ("Item 'str' can't be child 1 of 'BinaryOperation'. "
            "The valid format is: 'DataNode, DataNode'." in str(excinfo.value))


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
    a = BinaryOperation.create(BinaryOperation.Operator.ADD,ref1, ref2)

    assert a.operator is BinaryOperation.Operator.ADD
    assert a.children[0] is ref1
    assert a.children[1] is ref2

def test_binaryop_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    sym2 = ScalarTypeSymbol("x2", INT_TYPE)
    ref2 = ScalarReference(sym2)
    a = BinaryOperation.create(BinaryOperation.Operator.ADD,ref1, ref2)

    assert ("BinaryOperation[operator:'ADD']"
            in a.node_str())

def test_unaryop_create_badargs():
    with pytest.raises(GenerationError) as excinfo:
        a = UnaryOperation.create(123, [])
    assert (" operator argument in create method of UnaryOperation class should "
            "be a PSyIR UnaryOperation Operator but found 'int'." in str(excinfo.value))

    with pytest.raises(GenerationError) as excinfo:
        a = UnaryOperation.create(UnaryOperation.Operator.MINUS, 212)
    assert (" Item 'int' can't be child 0 of 'UnaryOperation'. "
            "The valid format is: 'DataNode'." in str(excinfo.value))

def test_unaryop_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)
    assert UnaryOperation._validate_child(0, "123") is False
    assert UnaryOperation._validate_child(0, ref1) is True
    assert UnaryOperation._validate_child(1, ref1) is False

def test_unaryop_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    a = UnaryOperation.create(UnaryOperation.Operator.MINUS, ref1)
    assert a.operator is UnaryOperation.Operator.MINUS
    assert a.children[0] is ref1

def test_unaryop_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref1 = ScalarReference(sym)

    a = UnaryOperation.create(UnaryOperation.Operator.MINUS, ref1)

    assert "UnaryOperation[operator:'MINUS']" in a.node_str()
