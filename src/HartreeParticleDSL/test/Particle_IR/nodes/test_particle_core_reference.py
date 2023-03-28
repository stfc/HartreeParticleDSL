import pytest

from HartreeParticleDSL.Particle_IR.nodes.particle_core_reference import ParticleCoreReference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, StructureType
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

def test_pcr():
    struct_type = StructureType()
    structure = StructureSymbol(name="structure1", datatype=struct_type)

    with pytest.raises(TypeError) as excinfo:
        ParticleCoreReference("x", 0)
    assert ("Attempted to make a ParticleCoreReference to a "
            "non-StructureSymbol. Got <class 'str'> as input."
            in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        ParticleCoreReference(structure, "abc")
    assert ("Attempted to make a ParticleCoreReference with a "
            "non-Member access. Got <class 'str'> as input."
            in str(excinfo.value))
    mem = Member("x")
    a = ParticleCoreReference(structure, mem)
    assert a.symbol is structure
    assert a.member is mem

    correct = "ParticleCoreReference[structure1: Member[x]]"
    assert correct == a.node_str()
