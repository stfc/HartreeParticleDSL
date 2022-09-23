from __future__ import annotations

from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import *
from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import *

import ast
import inspect
import textwrap
import pytest

def test_pir_cabana_visit_break():
    def a():
        while(True):
            break
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    a = ast_to_pir_visitor()
    pir = a.visit(c)

    b = Cabana_PIR_Visitor(None)
    out = b(pir)
    correct = '''void a(){
    while(True){
        break;
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_return():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        return
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    return;
}
'''
    assert correct == out

    def b():
        return 1
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''int b(){
    return 1;
}
'''
    assert correct == out

    def d(a: c_float):
        return a
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''float d(float a){
    return a;
}
'''
    assert correct == out


def test_pir_cabana_visit_ifdef():
    def a():
        create_variable(c_int, b)
        if True:
            b = 1
        else:
            b = 2

    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)

    correct = '''void a(){
    int b;
    if(True){
        b = 1;
    }else{
        b = 2;
    }
}
'''
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    pir = astpir.visit(c)
    out = cpir(pir)
    assert correct == out

    #TODO More tests here for incompleted functionality

def test_pir_cabana_visit_loop():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        create_variable(c_int, b)
        for b in range(0, 3):
            b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    int b;
    for(b = 0; b < 3; b += 1){
        b = (b + 1);
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_while():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        while True:
            break
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    while(True){
        break;
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_operations():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
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
    correct = '''void a(int b, double c, bool d){
    b = (b + 1);
    c = (c - 1.0);
    b = (b * 2);
    c = (c / 2.0);
    d = (b < 2);
    d = (b <= 3);
    d = (b > 4);
    d = (b >= 1);
    d = (b == 2);
    d = (d && d);
    d = (d || d);
    d = (!d);
    b = -1;
}
'''

    assert correct == out

