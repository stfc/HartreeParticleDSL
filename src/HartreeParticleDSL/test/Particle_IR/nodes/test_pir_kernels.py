import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.kernels import PairwiseKernel, PerPartKernel, \
                                                         MainKernel
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable

def test_kern_validate_child():
    assert PairwiseKernel._validate_child(0, Body()) == True
    assert PairwiseKernel._validate_child(0, "asd") == False
    assert PairwiseKernel._validate_child(1, Body()) == False

def test_pairwise_init():
    pk1 = PairwiseKernel("kernel")

    assert pk1.name == "kernel"
    pk1.name = "kernel_rename"
    assert pk1.name == "kernel_rename"

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert ("Expected PairwiseKernel name to be a str but got <class 'int'>." 
            in str(excinfo.value))

def test_pairwise_arguments_setter():
    pk1 = PairwiseKernel("kernel")

    with pytest.raises(IRGenerationError) as excinfo:
        pk1.arguments = []
    assert ("Pairwise kernel requires at least three arguments, "
            "but only got 0." == str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        pk1.arguments = ["s", "a", "b"]
    assert ("Each argument must be a Reference, but "
            "found <class 'str'>" in str(excinfo.value))

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))
    pk1.arguments = [arg1, arg2, arg3]
    assert pk1.arguments[0] is arg1
    assert pk1.arguments[1] is arg2
    assert pk1.arguments[2] is arg3


def test_pairwise_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))

    pk1 = PairwiseKernel.create("Kernel", [arg1, arg2, arg3], [assign])

    assert pk1.name == "Kernel"
    assert pk1.arguments[0] is arg1
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)

def test_pairwise_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    arg3 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))

    pk1 = PairwiseKernel.create("Kernel", [arg1, arg2, arg3], [assign])
    correct = '''PairwiseKernel[ScalarReference[y], ScalarReference[z], ScalarReference[a]: Body[
    Assignment[ScalarReference[x], Literal['25', Scalar<INTEGER, SINGLE>]]
] End Body]'''
    assert correct == pk1.node_str()

def test_perpart_init():
    pk1 = PerPartKernel("kernel")

    assert pk1.name == "kernel"
    pk1.name = "kernel_rename"
    assert pk1.name == "kernel_rename"

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert ("Expected PerPartKernel name to be a str but got <class 'int'>." 
            in str(excinfo.value))

def test_perpart_arguments_setter():
    pk1 = PerPartKernel("kernel")

    with pytest.raises(IRGenerationError) as excinfo:
        pk1.arguments = []
    assert ("Perpart kernel requires at least two arguments, "
            "but only got 0." == str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        pk1.arguments = ["s", "a"]
    assert ("Each argument must be a Reference, but "
            "found <class 'str'>" in str(excinfo.value))

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))
    pk1.arguments = [arg1, arg2]
    assert pk1.arguments[0] is arg1
    assert pk1.arguments[1] is arg2

def test_perpart_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = PerPartKernel.create("Kernel", [arg1, arg2], [assign])

    assert pk1.name == "Kernel"
    assert pk1.arguments[0] is arg1
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)

def test_perpart_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = PerPartKernel.create("Kernel", [arg1, arg2], [assign])
    correct = '''PerPartKernel[ScalarReference[y], ScalarReference[z]: Body[
    Assignment[ScalarReference[x], Literal['25', Scalar<INTEGER, SINGLE>]]
] End Body]'''
    assert correct == pk1.node_str()

def test_main_init():
    pk1 = MainKernel("kernel")

    assert pk1.name == "kernel"
    pk1.name = "kernel_rename"
    assert pk1.name == "kernel_rename"

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert ("Expected MainKernel name to be a str but got <class 'int'>." 
            in str(excinfo.value))

def test_main_arguments_setter():
    pk1 = MainKernel("kernel")

    with pytest.raises(IRGenerationError) as excinfo:
        pk1.arguments = ["s", "a"]
    assert ("Main kernel requires zero arguments, but got 2."
            in str(excinfo.value))
    assert len(pk1.arguments) == 0

def test_main_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 

    pk1 = MainKernel.create("Kernel", [assign])

    assert pk1.name == "Kernel"
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)
