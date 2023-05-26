from HartreeParticleDSL.Particle_IR.nodes.statement import Return
from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE


def test_return_validate_child():
    l = Literal("1", INT_TYPE)
    assert Return._validate_child(0, l) == True
    assert Return._validate_child(1, l) == False

