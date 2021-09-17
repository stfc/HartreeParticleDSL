from HartreeParticleDSL.backends.C_AOS.visitors import *
from HartreeParticleDSL.backends.C_AOS.C_AOS import *
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError

def test_c_visitor_visit_Str(capsys):
    '''Test the visit_Str function in c_visitor'''
    v = c_visitor()
    def a():
        b = "string"
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    captured = capsys.readouterr()
    v.visit(c)
    captured = capsys.readouterr()
    assert "\"string\"" in captured.out

def test_c_visitor_visit_str(capsys):
    '''Test the visit_str function in c_visitor'''
    v = c_visitor()
    def a():
        b = part.str
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = part->str;" in captured.out

def test_c_visitor_visit_int(capsys):
    '''Test the visit_int function in c_visitor'''
    v = c_visitor()
    def a():
        b = 13
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = 13;" in captured.out

def test_c_visitor_visit_Add(capsys):
    '''Test the visit_Add function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 + 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "+" in captured.out

def test_c_visitor_visit_Mult(capsys):
    '''Test the visit_Mult function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 * 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "*" in captured.out
    
def test_c_visitor_visit_Sub(capsys):
    '''Test the visit_Sub function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "-" in captured.out

def test_c_visitor_visit_LtE(capsys):
    '''Test the visit_LtE function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 <= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "<=" in captured.out

def test_c_visitor_visit_GtE(capsys):
    '''Test the visit_GtE function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 >= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert ">=" in captured.out

def test_c_visitor_visit_Lt(capsys):
    '''Test the visit_Lt function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 < 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "<" in captured.out

def test_c_visitor_visit_Gt(capsys):
    '''Test the visit_Gt function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 > 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert ">" in captured.out


def test_c_visitor_visit_USub(capsys):
    '''Test the visit_USub function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 + -a
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "1 +  -a" in captured.out

def test_c_visitor_visit_UnaryOp(capsys):
    '''Test the visit_UnaryOp function in c_visitor'''
    v = c_visitor()
    def a():
        b = -3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b =  -3" in captured.out

def test_c_visitor_visit_Compare(capsys):
    '''Test the visit_Compare function in c_visitor'''
    v = c_visitor()
    def a():
        b = (2 <= 3)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "2 <= 3" in captured.out

def test_c_visitor_visit_BinOp(capsys):
    '''Test the visit_BinOp function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = ( 1 - 2 )" in captured.out

def test_c_visitor_visit_Name(capsys):
    '''Test the visit_Name function in c_visitor'''
    v = c_visitor()
    def a():
        d.b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "d" in captured.out


def test_c_visitor_visit_Attribute(capsys):
    '''Test the visit_Attribute function in c_visitor'''
    v = c_visitor()
    def a():
        d.b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "d->b" in captured.out


def test_c_visitor_visit_Num(capsys):
    '''Test the visit_Num function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "1.2" in captured.out

def test_c_visitor_visit_Assign(capsys):
    '''Test the visit_Assign function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = 1.2;" in captured.out

def test_c_visitor_visit_arg(capsys):
    '''Test the visit_arg function in c_visitor'''
    v = c_visitor()
    def a(arg):
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "arg" in captured.out


def test_c_visitor_visit_arguments(capsys):
    '''Test the visit_arguments function in c_visitor'''
    v = c_visitor()
    def a(arg, arg2):
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "( struct part *arg, struct part *arg2 )" in captured.out

def test_c_visitor_visit_FunctionDef(capsys):
    '''Test the visit_FunctionDef function in c_visitor'''
    v = c_visitor()
    def a():
        b = 1.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "void a(  )\n{\n" in captured.out
    assert "}\n" in captured.out


def test_c_visitor_visit_If(capsys):
    '''Test the visit_If function in c_visitor'''
    v = c_visitor()
    def a():
        b = 0
        if( c > 1 ):
            b = 1
        elif( c <= -1):
            b = 2
        else:
            b = 3
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "    if(  ( c > 1 )  ){\n" in captured.out
    assert "}else if(  ( c <=  -1 )  ){\n" in captured.out
    assert "}else{\n" in captured.out

#Fails atm
def test_c_visitor_visit_Call(capsys):
    '''Test the visit_Call function in c_visitor'''
    v = c_visitor()
    def a():
        func(arg1)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    print(captured.out)
    assert "func( arg1 );" in captured.out

def test_c_visitor_visit_Module(capsys):
    '''Test the visit_Module function in c_visitor'''
    v = c_visitor()
    def a():
        asd.fgh = 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "asd->fgh" in captured.out

def test_c_visitor_visit_Index(capsys):
    '''Test the visit_Index function in c_visitor'''
    v = c_visitor()
    def a():
        d = s[2]
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "s[2]" in captured.out

def test_c_visitor_visit_Subscript(capsys):
    '''Test the visit_Subscript function in c_visitor'''
    v = c_visitor()
    def a():
        s = d[1]
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "[1]" in captured.out

def test_c_visitor_visit_For(capsys):
    '''Test the visit_For function in c_visitor'''
    v = c_visitor()
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
    v.visit(c)
    captured = capsys.readouterr()
    print(captured.out)
    assert "for( int i = 0; i < 3; i = i + 1)\n    {" in captured.out
    assert "for( int j = 3; j < 5; j = j + 1)\n    {" in captured.out
    assert "for( int k = 3; k >= 0; k = k + -1)\n    {" in captured.out
    assert "for( int s = 3; s < 9; s = s + 2)\n    {" in captured.out
    def b():
        for i in range(3):
            d[i] = 0.2
        else:
            x = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    with pytest.raises(IllegalLoopError) as excinfo:
        v.visit(c)
    assert "Else clauses on Loops are not supported in C_AOS" in str(excinfo.value)
    def f():
        for i in a:
            d[i] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    with pytest.raises(IllegalLoopError) as excinfo:
        v.visit(c)
    assert "Only range loops are supported in C_AOS" in str(excinfo.value)

    def g():
        for i in func(a):
            d[i] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(g)))
    with pytest.raises(IllegalLoopError) as excinfo:
        v.visit(c)
    assert "Only range loops are supported in C_AOS" in str(excinfo.value)


def test_c_visitor_visit_While(capsys):
    '''Test the visit_While function in c_visitor'''
    v = c_visitor()
    def a():
        while b < 3:
            b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "while( ( b < 3 ) ){\n" in captured.out

def test_c_visitor_visit_Expr(capsys):
    '''Test the visit_Expr function in c_visitor'''
    v = c_visitor()
    def a():
        b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = ( b + 1 );" in captured.out

def test_c_visitor_generic_visit(capsys):
    '''Test the generic_visit function in c_visitor'''
    v = c_visitor()
    with pytest.raises(UnsupportedCodeError) as excinfo:
        v.visit(v)
    assert f"Found unsupported node of type {type(v)}" in str(excinfo.value)

def test_c_pairwise_visitor_visit_arguments(capsys):
    '''Test the visit_arguments function in the c_pairwise_visitor'''
    v = c_pairwise_visitor()
    def a():
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        v.visit(c)
    assert "Pairwise function must have 4 arguments for C_AOS backend" in str(excinfo.value)
    capsys.readouterr()
    def f(part1, part2, r2, config):
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "struct part *part1, struct part *part2, double r2, struct config_type *config" in captured.out

# No tests yet for c_main_visitor, definitely needs a rethink anyway.
def test_c_main_visitor_visit_Expr(capsys):
    '''Test the visit_Expr function in c_main_visitor'''
    backend = C_AOS()
    v = c_main_visitor(backend)
    def main():
        b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "b = ( b + 1 );" in captured.out

def kern4(part1, part2, r2, config):
    part1.a = part1.a + 2.0

def test_c_main_visit_Call(capsys):
    '''Test the visit_Call function in c_main_visitor'''
    backend = C_AOS()
    kernel = kernels.pairwise_interaction(kern4)
    v = c_main_visitor(backend)
    def main():
        invoke(kern4)
        cleanup()
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "/* INVOKE generated for kern4 */" in captured.out
    assert "/* End of INVOKE generated for kern4 */" in captured.out
    assert "    free(config);\n    free(parts);" in captured.out


def test_c_perpart_visitor_visit_arguments(capsys):
    '''Test the visit_arguments function in the c_perpart_visitor'''
    v = c_perpart_visitor()
    def a():
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        v.visit(c)
    assert "Per part function must have 2 arguments for C_AOS backend" in str(excinfo.value)
    capsys.readouterr()
    def f(part1, config):
        b = 0
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    v.visit(c)
    captured = capsys.readouterr()
    assert "struct part *part1, struct config_type *config" in captured.out
