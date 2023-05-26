import pytest

from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.nodes.array_member import ArrayMember

def test_array_member():
    one = Literal("1", INT_TYPE)
    two = Literal("2", INT_TYPE)

    am = ArrayMember("x", [one, two])

    assert am.children[0] is one
    assert am.children[1] is two

    correct = "ArrayMember[x: (Literal[value:'1', Scalar<INTEGER, SINGLE>], Literal[value:'2', Scalar<INTEGER, SINGLE>])]"
    assert correct == am.node_str()

    assert am.is_array

    assert am.indices[0] is one
    assert am.indices[1] is two

