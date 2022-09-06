import pytest

from HartreeParticleDSL.Particle_IR.nodes.member import Member

def test_member():
    x = Member("myname")
    assert x.name == "myname"
    assert x.is_array == False
