from HartreeParticleDSL.backends.FDPS_backend.FDPS_visitors import *
from HartreeParticleDSL.backends.FDPS_backend.FDPS import *
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
from HartreeParticleDSL.HartreeParticleDSL import set_backend

def test_fdps_visitor_visit_Attribute():
    '''Test the visit_Attribute function in fdps_visitor'''
    aos = FDPS()
    aos.disable_variable_checks()
    set_backend(aos)
    v = fdps_visitor(aos)
    def a():
        d.b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "d.b" in out
    def aa():
        d.e.f = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(aa)))
    out = v.visit(c)
    assert "d.e.f" in out
    def ab():
        d.e[1] = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ab)))
    out = v.visit(c)
    assert "d.e[1]" in out
    def a():
        d.b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "d.b" in out
    def aa():
        d.e.f = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(aa)))
    out = v.visit(c)
    assert "d.e.f" in out
    def ab():
        d.e[1] = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ab)))
    out = v.visit(c)
    assert "d.e[1]" in out
    def ac():
        d.e[1][2][3][4] = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ac)))
    out = v.visit(c)
    assert "d.e[1][2][3][4]" in out
    def ad():
        d[1].e[2] = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ad)))
    out = v.visit(c)
    assert "d[1].e[2]" in out
    def ae():
        d[1][3].e = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ae)))
    out = v.visit(c)
    assert "d[1][3].e" in out
    def af():
        d.e[2].f = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(af)))
    out = v.visit(c)
    assert "d.e[2].f" in out

    def ag():
        d.core_part.position.x = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(ag)))
    out = v.visit(c)
    assert "d.core_part.position.x = 1" in out


def test_fdps_main_visitor_visit_Expr():
    '''Test the visit_Expr function in the fdps_main_visitor'''
    aos = FDPS()
    aos.disable_variable_checks()
    set_backend(aos)
    v = fdps_main_visitor(aos)
    def main():
        b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    out = v.visit(c)
    assert "b = ( b + 1 );" in out

def test_fdps_perpart_visitor_visit_arguments():
    aos = FDPS()
    aos.disable_variable_checks()
    set_backend(aos)
    v = fdps_perpart_visitor(aos)
    def func(p1):
        a = a - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(func)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        out = v.visit(c)
    assert "Per part function must have 2 arguments for FDPS backend" in str(excinfo.value)
