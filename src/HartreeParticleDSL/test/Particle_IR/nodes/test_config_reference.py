import pytest

from HartreeParticleDSL.Particle_IR.nodes.config_reference import ConfigReference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, StructureType
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

def test_config_reference():
    struct_type = StructureType()
    structure = StructureSymbol(name="structure1", datatype=struct_type)

    with pytest.raises(TypeError) as excinfo:
        ConfigReference("x", 0)
    assert ("Attempted to make a ConfigReference to a "
            "non-StructureSymbol. Got <class 'str'> as input."
            in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        ConfigReference(structure, "abc")
    assert ("Attempted to make a ConfigReference with a "
            "non-Member access. Got <class 'str'> as input."
            in str(excinfo.value))
    mem = Member("x")
    a = ConfigReference(structure, mem)
    assert a.symbol is structure
    assert a.member is mem

    correct = "ConfigReference[structure1: Member[x]]"
    assert correct == a.node_str()
