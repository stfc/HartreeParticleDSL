from __future__ import annotations

from HartreeParticleDSL.HartreeParticleDSL import set_mpi
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.backends.Fortran_base.pir_to_fortran_visitor import *
from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import *

from HartreeParticleDSL.inbuilt_kernels.boundaries.periodic_boundaries import periodic_boundaries

from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol
from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, ScalarType, \
        StructureType, ArrayType, PointerType, \
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, reset_part_and_config

from psyclone.psyir.symbols import Symbol

import ast
import inspect
import textwrap
import pytest

def test_pir_fortran_visit_break():
    def a():
        while(True):
            break
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    a = ast_to_pir_visitor()
    pir = a.visit(c)

    b = Fortran_PIR_Writer()
    out = b(pir)
    correct = '''Subroutine a()
    while(.True.) do
        EXIT
    enddo
End Subroutine a'''
    assert correct == out

@pytest.mark.xfail
def test_pir_fortran_visit_return():
    # Not yet implemented.
    assert False

def test_pir_fortran_visit_ifdef():
    def a():
        create_variable(c_int, b)
        if True:
            b = 1
        else:
            b = 2

    astpir = ast_to_pir_visitor()
    cpir = Fortran_PIR_Writer()

    correct = '''Subroutine a()
    Integer :: b
    if(.True.) then
        b = 1
    else
        b = 2
    endif
End Subroutine a'''
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    pir = astpir.visit(c)
    out = cpir(pir)
    assert correct == out

    def b():
        create_variable(c_int, b)
        if b < 1:
            b = 1
        elif b < 2:
            b = 2
        elif b < 3:
            b = 3
        else:
            b = 4
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine b()
    Integer :: b
    if(b < 1) then
        b = 1
    else if(b < 2) then
        b = 2
    else if(b < 3) then
        b = 3
    else
        b = 4
    endif
End Subroutine b'''
    assert correct == out

def test_pir_fortran_visit_loop():
    astpir = ast_to_pir_visitor()
    cpir = Fortran_PIR_Writer()
    def a():
        create_variable(c_int, b)
        for b in range(0, 3):
            b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine a()
    Integer :: b
    do b = 0, 3, 1
        b = b + 1
    enddo
End Subroutine a'''
    assert correct == out

def test_pir_fortran_visit_operations():
    astpir = ast_to_pir_visitor()
    cpir = Fortran_PIR_Writer()
    def a(b: c_int, c: c_double, d: c_bool):
        b = b + 1
        c = c - 1.0
        b = b * 2
        c = c / 2.0

        d = b < 2
        d = b <= 3
        d = b > 4
        d = b >= 1
        d = b == 2
        d = d and d
        d = d or d

        d = not d
        b = -1


    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine a(b, c, d)
    Integer :: b
    Real*8 :: c
    Logical :: d
    b = b + 1
    c = c - 1.0
    b = b * 2
    c = c / 2.0
    d = b < 2
    d = b <= 3
    d = b > 4
    d = b >= 1
    d = b == 2
    d = d .and. d
    d = d .or. d
    d = .not.d
    b = -1
End Subroutine a'''

    assert correct == out

def test_pir_fortran_visit_call():
    astpir = ast_to_pir_visitor()
    cpir = Fortran_PIR_Writer()

    def d():
        create_variable(c_int, x)
        mycall(x)
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine d()
    Integer :: x
    call mycall(x)
End Subroutine d'''
    assert correct == out

def test_pir_fortran_visit_members():
    astpir = ast_to_pir_visitor()
    cpir = Fortran_PIR_Writer()
    mystruc1 = StructureType()
    mystruc1.add("d", INT_TYPE, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc1"]=mystruc1
    def a():
        create_variable(mystruc1, d)
        d.d = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine a()
    Type(mystruc1) :: d
    d%d = 1
End Subroutine a'''
    assert correct == out
    del(type_mapping_str["mystruc1"])

    mystruc2 = StructureType()
    substruc = StructureType()
    substruc.add("e", INT_TYPE, Symbol.Visibility.PUBLIC)
    mystruc2.add("d", substruc, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2
    def b():
        create_variable(mystruc2, b)
        b.d.e = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    pir = astpir.visit(c)
    out = cpir(pir)

    correct = '''Subroutine b()
    Type(mystruc2) :: b
    b%d%e = 1
End Subroutine b'''
    assert correct == out
    del(type_mapping_str["mystruc2"])

    mystruc2 = StructureType()
    substruc1 = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC])
    mystruc2.add("c", substruc1, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2
    def d():
        create_variable(mystruc2, b)
        b.c[0] = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''Subroutine d()
    Type(mystruc2) :: b
    b%c(0) = 1
End Subroutine d'''
    assert correct == out
    del(type_mapping_str["mystruc2"])

@pytest.mark.xfail
def test_pir_fortran_pointer_reference():
    a = PointerReference(PointerSymbol("a", PointerType(INT_TYPE)))
    cpir = Fortran_PIR_Writer()
    out = cpir(a)

@pytest.mark.xfail
def test_pir_fortran_pointer_symbol():
    a = PointerSymbol("a", PointerType(INT_TYPE))
    cpir = Fortran_PIR_Writer()
    out = cpir.visit_pointersymbol_node(a)

@pytest.mark.xfail
def test_pir_fortran_auto_symbol():
    a = AutoSymbol("a", "fake_call()")
    cpir = Fortran_PIR_Writer()
    out = cpir.visit_autosymbol_node(a)

def test_pir_fortran_arrayreference():
    cpir = Fortran_PIR_Writer()
    ar = ArrayReference(ArraySymbol("x", ArrayType(INT_TYPE, [2,2])), [Literal("2", INT_TYPE), Literal("2", INT_TYPE)])
    out = cpir(ar)
    assert out == "x(2,2)"

    ar = ArrayReference(ArraySymbol("ea", ArrayType(DOUBLE_TYPE, [24])), [Literal("2", INT_TYPE)])
    out = cpir(ar)
    assert out == "ea(2)"

@pytest.mark.xfail
def test_pir_fortran_config_reference():

    config = type_mapping_str["config"]
    cpir = Fortran_PIR_Writer()

    def x(arg: part, conf: config):
        conf.space = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    v = pir_perpart_visitor()
    pir = v.visit(c)
    out = cpir(pir)
    correct = '''struct x_functor{
    config_struct_type _conf;

    KOKKOS_INLINE_FUNCTION
     x_functor( config_struct_type conf):
    _conf(conf){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int a) const{
        _conf(0).space = 1;
    }
};
'''
    assert out == correct

    def main():
        config.space = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    v = pir_main_visitor()
    pir = v.visit(c)
    out = cpir(pir)

    correct = "config.config_host(0).space = 1;\n"
    assert correct in out

@pytest.mark.xfail
def test_pir_fortran_particle_references_and_perpart():
    assert False

@pytest.mark.xfail
def test_pir_fortran_pairwise():
    assert False

@pytest.mark.xfail
def test_pir_fortran_sink_boundary():
    assert False

@pytest.mark.xfail
def test_pir_fortran_source_boundary():
    assert False

@pytest.mark.xfail
def test_pir_fortran_mainkernel():
    assert False


