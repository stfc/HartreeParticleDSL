import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.kernels import PairwiseKernel, PerPartKernel, \
                                                         MainKernel, SourceBoundaryKernel, \
                                                         SinkBoundaryKernel
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.nodes.particle_position_reference import ParticlePositionReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, StructureType
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable

def test_kern_validate_child():
    assert PairwiseKernel._validate_child(0, Body()) == True
    assert PairwiseKernel._validate_child(0, "asd") == False
    assert PairwiseKernel._validate_child(1, Body()) == False

def test_kernel_update_position():
    pk1 = PairwiseKernel("kernel")
    struct_type = StructureType()
    structure = StructureSymbol(name="structure1", datatype=struct_type)

    a = ParticlePositionReference(structure, 0)

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    assign = Assignment.create(ref_lhs, a) 
    pk1.body.addchild(assign)
    assert pk1.does_update_position() == False

    b = ParticlePositionReference(structure, 1)
    sym2 = ScalarTypeSymbol("y", INT_TYPE)
    ref_rhs = ScalarReference(sym2)
    assign2 = Assignment.create(b, ref_rhs)
    pk1.body.addchild(assign2)
    assert pk1.does_update_position() == True

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
    correct = '''PairwiseKernel[ScalarReference[name:'y'], ScalarReference[name:'z'], ScalarReference[name:'a']: Body[
    Assignment[ScalarReference[name:'x'], Literal['25', Scalar<INTEGER, SINGLE>]]
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
    assert pk1.arguments == [arg1, arg2]

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
    correct = '''PerPartKernel[ScalarReference[name:'y'], ScalarReference[name:'z']: Body[
    Assignment[ScalarReference[name:'x'], Literal['25', Scalar<INTEGER, SINGLE>]]
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

def test_source_boundary():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 
    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = SourceBoundaryKernel.create("myname", 100000, [arg1, arg2], [assign])
    
    assert pk1.name == "myname"
    assert len(pk1.arguments) == 2
    assert pk1.arguments[0] is arg1
    assert pk1.arguments[1] is arg2
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)

    pk1.name = "newname"
    assert pk1.name == "newname"

    with pytest.raises(IRGenerationError) as excinfo:
        pk1.arguments = []
    assert "Source boundary kernel requires at least two arguments, but only got 0." in str(excinfo.value)

    with pytest.raises(TypeError) as excinfo:
        pk1.arguments = ["a", 1]
    assert "Each argument must be a Reference, but found <class 'str'>." in str(excinfo.value)

    assert pk1.source_count == 100000
    pk1.source_count = 123
    assert pk1.source_count == 123

    with pytest.raises(TypeError) as excinfo:
        pk1.source_count = "a"
    assert "Expected SourceBoundaryKernel source_count to be an int but got <class 'str'>." in str(excinfo.value)

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert "Expected SourceBoundaryKernel name to be a str but got <class 'int'>." in str(excinfo.value)

def test_sink_boundary():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs) 
    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = SinkBoundaryKernel.create("myname", [arg1, arg2], [assign])
    
    assert pk1.name == "myname"
    assert len(pk1.arguments) == 2
    assert pk1.arguments[0] is arg1
    assert pk1.arguments[1] is arg2
    assert pk1.children[0].children[0] is assign
    assert pk1.body.children[0] is assign
    assert isinstance(pk1.symbol_table, SymbolTable)

    pk1.name = "newname"
    assert pk1.name == "newname"

    with pytest.raises(IRGenerationError) as excinfo:
        pk1.arguments = []
    assert "Sink boundary kernel requires at least two arguments, but only got 0." in str(excinfo.value)

    with pytest.raises(TypeError) as excinfo:
        pk1.arguments = ["a", 1]
    assert "Each argument must be a Reference, but found <class 'str'>." in str(excinfo.value)

    with pytest.raises(TypeError) as excinfo:
        pk1.name = 123
    assert "Expected SinkBoundaryKernel name to be a str but got <class 'int'>." in str(excinfo.value)
