import pytest

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.kernels import PerPartKernel, MainKernel
from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.datatypes.datatype import INT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable

def test_invoke_validate_child():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs)

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = PerPartKernel.create("Kernel", [arg1, arg2], [assign])
    assert Invoke._validate_child(0, pk1) == True
    assert Invoke._validate_child(1, pk1) == True

    pk2 = MainKernel.create("main", [])
    assert Invoke._validate_child(0, pk2) == False

def test_invoke_creates():
    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs)

    arg1 = ScalarReference(ScalarTypeSymbol("y", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("z", INT_TYPE))

    pk1 = PerPartKernel.create("Kernel", [arg1, arg2], [assign])

    sym = ScalarTypeSymbol("x", INT_TYPE)
    ref_lhs = ScalarReference(sym)
    lit_rhs = Literal("25", INT_TYPE)
    assign = Assignment.create(ref_lhs, lit_rhs)

    arg1 = ScalarReference(ScalarTypeSymbol("a", INT_TYPE))
    arg2 = ScalarReference(ScalarTypeSymbol("b", INT_TYPE))

    pk2 = PerPartKernel.create("Kernel", [arg1, arg2], [assign])

    invoke1 = Invoke.create(pk1)

    assert invoke1.children[0] is pk1

    pk1.detach()

    invoke2 = Invoke.create([pk1, pk2])

    assert invoke2.children[0] is pk1
    assert invoke2.children[1] is pk2
