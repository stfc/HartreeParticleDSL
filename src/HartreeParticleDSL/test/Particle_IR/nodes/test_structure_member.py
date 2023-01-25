import pytest

from HartreeParticleDSL.Particle_IR.nodes.structure_member import StructureMember
from HartreeParticleDSL.Particle_IR.nodes.member import Member

def test_array_member():
    mem = Member("a")

    sm = StructureMember("x", mem)

    assert sm.children[0] is mem

    correct = "StructureMember[x: Member[a]]"
    assert correct == sm.node_str()

    assert sm.is_array == False

    assert sm.member is mem

    assert StructureMember._validate_child(0, "123") is False
