from HartreeParticleDSL.backends.C_AOS.visitors import *
from HartreeParticleDSL.backends.C_AOS.C_AOS import *
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL


def test_c_visitor_visit_Str():
    '''Test the visit_Str function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = "string"
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "\"string\"" in out

def test_c_visitor_visit_str():
    '''Test the visit_str function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    aos.variable_scope.add_variable("b", "c_double", False)
    def a():
        b = part.str
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = part.str;" in out
    aos.variable_scope.add_variable("yz", "c_int", False)
    def z():
        b = x.yz
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    out = v.visit(c)
    assert "b = x.yz" in out

def test_c_visitor_visit_int():
    '''Test the visit_int function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 13
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = 13;" in out

def test_c_visitor_visit_Add():
    '''Test the visit_Add function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 + 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "+" in out

def test_c_visitor_visit_Mult():
    '''Test the visit_Mult function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 * 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "*" in out
    
def test_c_visitor_visit_Sub():
    '''Test the visit_Sub function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "-" in out

def test_c_visitor_visit_Div():
    '''Test the visit_Div function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 / 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "/" in out


def test_c_visitor_visit_LtE():
    '''Test the visit_LtE function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 <= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "<=" in out

def test_c_visitor_visit_GtE():
    '''Test the visit_GtE function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 >= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert ">=" in out

def test_c_visitor_visit_Lt():
    '''Test the visit_Lt function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 < 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "<" in out

def test_c_visitor_visit_Gt():
    '''Test the visit_Gt function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 > 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert ">" in out


def test_c_visitor_visit_USub():
    '''Test the visit_USub function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 + -a
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "1 + -a" in out

def test_c_visitor_visit_UnaryOp():
    '''Test the visit_UnaryOp function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = -3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = -3" in out

def test_c_visitor_visit_Compare():
    '''Test the visit_Compare function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = (2 <= 3)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "2 <= 3" in out

def test_c_visitor_visit_BinOp():
    '''Test the visit_BinOp function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = ( 1 - 2 )" in out

def test_c_visitor_visit_And():
    '''Test the visit_And function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = z and y
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "&&" in out
    
def test_c_visitor_visit_Or():
    '''Test the visit_Or function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = z or y
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "||" in out

def test_c_visitor_visit_Not():
    '''Test the visit_Not function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = not y
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "!y" in out
    

def test_c_visitor_visit_BoolOp():
    '''Test the visit_BoolOp function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = z and y
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = ( z && y )" in out


def test_c_visitor_visit_Name():
    '''Test the visit_Name function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        d.b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "d" in out


def test_c_visitor_visit_Attribute():
    '''Test the visit_Attribute function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
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
    def abc():
        a.core_part.position.x = 3.0
        a.core_part.position.y = 4.0
        a.core_part.position.z = 5.0
    c = ast.parse(textwrap.dedent(inspect.getsource(abc)))
    out = v.visit(c)
    assert "a.core_part.position[0] = 3.0" in out
    assert "a.core_part.position[1] = 4.0" in out
    assert "a.core_part.position[2] = 5.0" in out

    def illegal():
        a.core_part.position.xd = 6.0
    c = ast.parse(textwrap.dedent(inspect.getsource(illegal)))
    with pytest.raises(InvalidNameError) as excinfo:
        out = v.visit(c)
    assert "The dimension argument should be x, y, or z" in str(excinfo.value)


def test_c_visitor_visit_Num():
    '''Test the visit_Num function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "1.2" in out

def test_c_visitor_visit_Assign():
    '''Test the visit_Assign function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = 1.2;" in out

def test_c_visitor_visit_arg():
    '''Test the visit_arg function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a(arg):
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "arg" in out


def test_c_visitor_visit_arguments():
    '''Test the visit_arguments function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a(arg, arg2):
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "( struct part *arg, struct part *arg2 )" in out

def test_c_visitor_visit_FunctionDef():
    '''Test the visit_FunctionDef function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "void a(  )\n{\n" in out
    assert "}\n" in out


def test_c_visitor_visit_If():
    '''Test the visit_If function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = 0
        if( c > 1 ):
            b = 1
        elif( c <= -1):
            b = 2
        else:
            b = 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "    if(  ( c > 1 )  ){\n" in out
    assert "}else if(  ( c <= -1 )  ){\n" in out
    assert "}else{\n" in out

def test_c_visitor_visit_Call():
    '''Test the visit_Call function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        func(arg1)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "func( arg1 );" in out
    def b():
        func(arg1, arg2=32)
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert "func( arg1, arg2=32 );" in out

    def b():
        func(arg1, arg2)
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert "func( arg1, arg2 );" in out

    def f():
        mod1.mod2.func()
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    out = v.visit(c)
    assert "mod1.mod2.func( )" in out

def test_c_visitor_visit_Module():
    '''Test the visit_Module function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        asd.fgh = 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "asd.fgh" in out

def test_c_visitor_visit_Index():
    '''Test the visit_Index function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        d = s[2]
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "s[2]" in out

def test_c_visitor_visit_Subscript():
    '''Test the visit_Subscript function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        s = d[1]
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "[1]" in out

def test_c_visitor_visit_For():
    '''Test the visit_For function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        for i in range(3):
            d[i] = 0.2
        for j in range(3,5):
            d[j] = 0.2
        for k in range(3,0,-1):
            d[k] = 0.2
        for s in range(3, 9, 2):
            d[s] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "for( int i = 0; i < 3; i = i + 1)\n    {" in out
    assert "for( int j = 3; j < 5; j = j + 1)\n    {" in out
    assert "for( int k = 3; k >= 0; k = k + -1)\n    {" in out
    assert "for( int s = 3; s < 9; s = s + 2)\n    {" in out
    def b():
        for i in range(3):
            d[i] = 0.2
        else:
            x = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    with pytest.raises(IllegalLoopError) as excinfo:
        out = v.visit(c)
    assert "Else clauses on Loops are not supported in C_AOS" in str(excinfo.value)
    def f():
        for i in a:
            d[i] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    with pytest.raises(IllegalLoopError) as excinfo:
        out = v.visit(c)
    assert "Only range loops are supported in C_AOS" in str(excinfo.value)

    def g():
        for i in func(a):
            d[i] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(g)))
    with pytest.raises(IllegalLoopError) as excinfo:
        out = v.visit(c)
    assert "Only range loops are supported in C_AOS" in str(excinfo.value)


def test_c_visitor_visit_While():
    '''Test the visit_While function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        while b < 3:
            b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "while( ( b < 3 ) ){\n" in out

def test_c_visitor_visit_Expr():
    '''Test the visit_Expr function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    def a():
        b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert "b = ( b + 1 );" in out

def test_c_visitor_generic_visit():
    '''Test the generic_visit function in c_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_visitor(aos)
    with pytest.raises(UnsupportedCodeError) as excinfo:
        out = v.visit(v)
    assert f"Found unsupported node of type {type(v)}" in str(excinfo.value)

def test_c_pairwise_visitor_visit_arguments():
    '''Test the visit_arguments function in the c_pairwise_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_pairwise_visitor(aos)
    def a():
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        out = v.visit(c)
    assert "Pairwise function must have 4 arguments for C_AOS backend" in str(excinfo.value)
    def f(part1, part2, r2, config):
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    out = v.visit(c)
    assert "struct part *part1, struct part *part2, double r2, struct config_type *config" in out

def test_c_main_visitor_visit_Expr():
    '''Test the visit_Expr function in c_main_visitor'''
    backend = C_AOS()
    HartreeParticleDSL.set_backend(backend)
    v = c_main_visitor(backend)
    def main():
        b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    out = v.visit(c)
    assert "b = ( b + 1 );" in out

def kern4(part1, part2, r2, config):
    part1.a = part1.a + 2.0

def test_c_main_visit_Call():
    '''Test the visit_Call function in c_main_visitor'''
    backend = C_AOS()
    HartreeParticleDSL.set_backend(backend)
    kernel = kernels.pairwise_interaction(kern4)
    v = c_main_visitor(backend)
    def main():
        invoke(kern4)
        cleanup()
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    out = v.visit(c)
    assert "/* INVOKE generated for kern4 */" in out
    assert "/* End of INVOKE generated for kern4 */" in out
    assert "    free(config);\n    free(parts);" in out


def test_c_perpart_visitor_visit_arguments():
    '''Test the visit_arguments function in the c_perpart_visitor'''
    aos = C_AOS()
    HartreeParticleDSL.set_backend(aos)
    v = c_perpart_visitor(aos)
    def a():
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        out = v.visit(c)
    assert "Per part function must have 2 arguments for C_AOS backend" in str(excinfo.value)
    def f(part1, config):
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    out = v.visit(c)
    assert "struct part *part1, struct config_type *config" in out
