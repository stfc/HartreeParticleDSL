import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from psyclone.errors import GenerationError

def test_assignment_construction():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)

    lit_rhs = Literal("25", INT_TYPE)

    assign = Assignment()

    with pytest.raises(GenerationError) as excinfo:
        assign.addchild("123")
    assert ("Item 'str' can't be child 0 of 'Assignment'. "
            "The valid format is: 'DataNode, DataNode'." in str(excinfo.value))

    assign.addchild(ref_lhs)

    with pytest.raises(GenerationError) as excinfo:
        assign.addchild("123")
    assert ("Item 'str' can't be child 1 of 'Assignment'. "
            "The valid format is: 'DataNode, DataNode'." in str(excinfo.value))

    assign.addchild(lit_rhs)

    assert assign.lhs is ref_lhs
    assert assign.rhs is lit_rhs


def test_assignment_create():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)

    lit_rhs = Literal("25", INT_TYPE)

    assign = Assignment.create(ref_lhs, lit_rhs)
    assert assign.lhs is ref_lhs
    assert assign.rhs is lit_rhs

def test_assignment_nodestr():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)

    lit_rhs = Literal("25", INT_TYPE)

    assign = Assignment.create(ref_lhs, lit_rhs)

    correct = "Assignment[ScalarReference[name:'x'], Literal[value:'25', Scalar<INTEGER, SINGLE>]]"
    assert correct == assign.node_str()
