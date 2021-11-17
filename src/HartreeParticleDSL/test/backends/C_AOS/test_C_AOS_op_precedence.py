from HartreeParticleDSL.backends.C_AOS.visitors import *
from HartreeParticleDSL.backends.C_AOS.C_AOS import *
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL

def test_plus_mul():
    '''Test the plus_mul order for C_AOS'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 + 2 * 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( 1 + ( 2 * 3 ) )" in out

def test_plus_div():
    '''Test the plus_div order for C_AOS'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 + 2 / 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( 1 + ( 2 / 3 ) )" in out

def test_mul_div():
    '''Test the mul_div order for C_AOS'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 * 2 / 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( ( 1 * 2 ) / 3 )" in out
    def a():
        b = 1 / 2 * 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( ( 1 / 2 ) * 3 )" in out

def test_bracket_plus_mul():
    '''Test how brackets affect the plus_mul order for C_AOS'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = (1 + 2) * 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( ( 1 + 2 ) * 3 )" in out

def test_bracket_plus_div():
    '''Test how brackets affect the plus_div order for C_AOS'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = (1 + 2) / 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( ( 1 + 2 ) / 3 )" in out

def test_land():
    '''Test how l_and works'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = z and y
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "z && y" in out

