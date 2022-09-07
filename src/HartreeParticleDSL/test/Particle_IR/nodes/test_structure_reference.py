import pytest

from HartreeParticleDSL.Particle_IR.nodes.structure_reference import StructureReference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE, StructureType
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

def test_pr():
    structure = StructureSymbol("structure1", StructureType())
    mem = Member("x")

    with pytest.raises(TypeError) as excinfo:
        StructureReference("x", "a")
    assert ("Attempted to make a StructureReference to a "
            "non-StructureSymbol. Got <class 'str'> as input."
            in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        StructureReference(structure, "b")
    assert ("Attempted to make a StructureReference with a non-Member access. "
            "Got <class 'str'> as input." in str(excinfo.value))

    a = StructureReference(structure, mem)
    assert a.symbol is structure
    assert a.member is mem

    correct = "StructureReference[structure1: Member[x]]"
    assert correct == a.node_str()
