import pytest

from HartreeParticleDSL.Particle_IR.nodes.particle_position_reference import ParticlePositionReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, StructureType
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

def test_ppr():
    struct_type = StructureType()
    structure = StructureSymbol(name="structure1", datatype=struct_type)

    with pytest.raises(TypeError) as excinfo:
        ParticlePositionReference("x", 0)
    assert ("Attempted to make a ParticlePositionReference to a "
            "non-StructureSymbol. Got <class 'str'> as input."
            in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        ParticlePositionReference(structure, "abc")
    assert ("Attempted to make a ParticlePositionReference with a "
            "non-int dimension. Got <class 'str'> as input."
            in str(excinfo.value))

    a = ParticlePositionReference(structure, 0)
    assert a.symbol is structure
    assert a.dimension == 0

    correct = "ParticlePositionReference[structure1: 0]"
    assert correct == a.node_str()
