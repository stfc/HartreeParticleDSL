import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.nodes.funcdef import FuncDef
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable


def test_funcdef_validate_child():
    assert FuncDef._validate_child(0, Body()) == True
    assert FuncDef._validate_child(0, "asd") == False
    assert FuncDef._validate_child(1, Body()) == False

def test_funcdef_init():
    pk1 = FuncDef("kernel")

    assert pk1.name == "kernel"
    pk1.name = "kernel_rename"
    assert pk1.name == "kernel_rename"

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert ("Expected FuncDef name to be a str but got <class 'int'>."
            in str(excinfo.value))

def test_funcdef_arguments_setter():
    pk1 = FuncDef("kernel")

    with pytest.raises(TypeError) as excinfo:
        pk1.arguments = ["s", "a", "b"]
    assert ("Each argument must be a DataNode, but "
            "found <class 'str'>" in str(excinfo.value))

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))
    pk1.arguments = [arg1, arg2, arg3]
    assert pk1.arguments[0] is arg1
    assert pk1.arguments[1] is arg2
    assert pk1.arguments[2] is arg3

def test_funcdef_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs)

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))

    pk1 = FuncDef.create("Kernel", [arg1, arg2, arg3], [assign])

    assert pk1.name == "Kernel"
    assert pk1.arguments[0] is arg1
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)

def test_funcdef_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs)

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))

    pk1 = FuncDef.create("Kernel", [arg1, arg2, arg3], [assign])
    correct = '''FuncDef[ScalarReference[name:'y'], ScalarReference[name:'z'], ScalarReference[name:'a']: Body[
    Assignment[ScalarReference[name:'x'], Literal[value:'25', Scalar<INTEGER, SINGLE>]]
] End Body]'''
    assert correct == pk1.node_str()
