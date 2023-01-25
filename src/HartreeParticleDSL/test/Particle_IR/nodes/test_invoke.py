import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.kernels import PerPartKernel, MainKernel
from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import STRING_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable

def test_invoke_validate_child():
    pk1 = Literal("mykern", STRING_TYPE)
    assert Invoke._validate_child(0, pk1) == True
    assert Invoke._validate_child(1, pk1) == True

    pk2 = MainKernel.create("main", [])
    assert Invoke._validate_child(0, pk2) == False

def test_invoke_creates():
    pk1 = Literal("mykern", STRING_TYPE)
    pk2 = Literal("mykern2", STRING_TYPE)

    invoke1 = Invoke.create(pk1)

    assert invoke1.children[0] is pk1

    pk1.detach()

    invoke2 = Invoke.create([pk1, pk2])

    assert invoke2.children[0] is pk1
    assert invoke2.children[1] is pk2
