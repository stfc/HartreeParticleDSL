import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.ifelse import IfElseBlock
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.call import Call


def test_ifelse_check_completeness():
    x = IfElseBlock()

    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert ("IfElseBlock must have 3 children but found only 0."
            in str(excinfo.value))

    # Override the children to meet the failure conditions that the Children
    # list implementation should usually prevent.
    x._children = []
    x.addchild("str")
    x.addchild("str2")
    x.addchild("str3")
    
    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert ("IfElseBlock first child must be a Node but found <class 'str'>."
            in str(excinfo.value))

    x.children[0] = BinaryOperation(BinaryOperation.BinaryOp.ADDITION)
    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert ("IfElseBlock second child must be a Body but found <class 'str'>."
            in str(excinfo.value))

    x.children[1] = Body()
    with pytest.raises(IRGenerationError) as excinfo:
        x._check_completeness()
    assert ("IfElseBlock third child must be a Body but found <class 'str'>."
            in str(excinfo.value))
    x.children[2] = Body()

    x._check_completeness()



def test_ifelse_validate_child():
    assert IfElseBlock._validate_child(0, "123") == False
    assert IfElseBlock._validate_child(1, "123") == False
    assert IfElseBlock._validate_child(2, "123") == False
    assert IfElseBlock._validate_child(0, BinaryOperation(BinaryOperation.BinaryOp.ADDITION)) == True
    assert IfElseBlock._validate_child(1, Body()) == True
    assert IfElseBlock._validate_child(2, Body()) == True

def test_ifelse_create():
    with pytest.raises(TypeError) as excinfo:
        IfElseBlock.create("str", [], [])
    assert ("The condition input needs to be a Node, but found <class 'str'>."
            in str(excinfo.value))

    body1 = Call("hello")
    body2 = Call("hello2")
    bo = BinaryOperation(BinaryOperation.BinaryOp.LESS_THAN_EQUAL)
    bo.addchild(Literal("1", INT_TYPE))
    bo.addchild(Literal("2", INT_TYPE))
    ifelse = IfElseBlock.create(bo, [body1], [body2])

    assert ifelse.condition is bo
    assert ifelse.ifbody.children[0] is body1
    assert ifelse.elsebody.children[0] is body2
    assert ifelse.is_else_if() is False

    bo2 = BinaryOperation(BinaryOperation.BinaryOp.LESS_THAN_EQUAL)
    ifelse2 = IfElseBlock.create(bo2, [Call("hello3")], [ifelse])

    assert ifelse2.elsebody.children[0] is ifelse
    assert len(ifelse2.elsebody.children) > 0
    assert isinstance(ifelse2.elsebody.children[0], IfElseBlock)
    assert ifelse2.is_else_if() is True

def test_ifelse_nodestr():
    body1 = Call("hello")
    body2 = Call("hello2")
    bo = BinaryOperation(BinaryOperation.BinaryOp.LESS_THAN_EQUAL)
    bo.addchild(Literal("1", INT_TYPE))
    bo.addchild(Literal("2", INT_TYPE))
    ifelse = IfElseBlock.create(bo, [body1], [body2])

    correct = '''IfElseBlock[BinaryOperation[BinaryOp.LESS_THAN_EQUAL: (Literal['1', Scalar<INTEGER, SINGLE>], Literal['2', Scalar<INTEGER, SINGLE>])]:
Body[
    Call[hello: ()]
] End Body
Body[
    Call[hello2: ()]
] End Body
]'''
    print(ifelse.node_str())
    assert correct == ifelse.node_str()
