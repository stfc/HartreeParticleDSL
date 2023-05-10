import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

from psyclone.errors import GenerationError

def test_body():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)

    lit_rhs = Literal("25", INT_TYPE)

    assign = Assignment.create(ref_lhs, lit_rhs)

    body = Body()

    with pytest.raises(GenerationError) as excinfo:
        body.addchild("123")
    assert ("Item 'str' can't be child 0 of 'Schedule'. The valid format is: '[Statement]*'" in str(excinfo.value))

    body.addchild(assign)

    correct = '''Body[
    Assignment[ScalarReference[x], Literal['25', Scalar<INTEGER, SINGLE>]]
] End Body'''
    assert correct == body.node_str()
