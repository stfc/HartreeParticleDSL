from HartreeParticleDSL.backends.Cabana_backend.Cabana_visitors import *
from HartreeParticleDSL.backends.Cabana_backend.Cabana import *
import ast
import inspect
import textwrap
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
from HartreeParticleDSL.HartreeParticleDSL import set_backend

def test_cabana_visitor_visit_Attribute():
    '''Test the visit_Attribute function in cabana_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_visitor(aos)
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
    def ag():
        core_part.position.x = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(ag)))
    out = v.visit(c)
    assert "core_part_position[0] = 1" in out

def test_cabana_pairwise_visitor_check_position():
    '''Test the check_position function in cabana_pairwise_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_pairwise_visitor(aos)
    assert v.check_position(None) is None
    v.addSlice(None)

def test_cabana_pairwise_visitor_visit_arguments():
    '''Test the visit_arguments function in cabana_pairwise_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_pairwise_visitor(aos)
    def a(arg1, arg2):
        d = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(UnsupportedCodeError) as excinfo:
        v.visit(c)
    assert "Cabana doesn't yet support pairwise operations" in str(excinfo.value)

def test_cabana_main_visitor_addSlice():
    '''Test the addSlice function in cabana_main_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_main_visitor(aos)
    assert len(v._slices) == 0
    v.addSlice("pot")
    assert v._slices[0] == "pot"
    v.addSlice("pot")
    assert len(v._slices) == 1

def test_cabana_main_visitor_visit_FunctionDef():
    '''Test the visit_FunctionDef function in cabana_main_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_main_visitor(aos)
    def main():
        a = 3
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    r = v.visit(c)
    assert "int main( int argc, char* argv[] )\n{\n" in r

def test_cabana_perpart_visitor_check_position():
    ''' Test the check_position function in cabana_perpart_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_perpart_visitor(aos)
    assert v.check_position(None) == None

def test_cabana_perpart_visitor_slices():
    ''' Test the addSlice and resetSlices functions in cabana_perpart_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_perpart_visitor(aos)
    assert len(v._slices) == 0
    v.addSlice("pot")
    assert v._slices[0] == "pot"
    v.addSlice("pot")
    assert len(v._slices) == 1
    v.resetSlices()
    assert len(v._slices) == 0

def test_cabana_perpart_visitor_visit_arguments():
    ''' Test the visit_arguments function in cabana_perpart_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_perpart_visitor(aos)
    def f1(arg1):
        a = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(f1)))
    with pytest.raises(IllegalArgumentCountError) as excinfo:
        v.visit(c)
    assert "Per part function must have 2 arguments for Cabana backend" in str(excinfo.value)
    def f2(arg1, arg2):
        a = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(f2)))
    x = v.visit(c.body[0].args)
    assert v._part1 == "arg1"
    assert v._config == "arg2"
    assert "const int arg1" in x

def test_cabana_perpart_visitor_visit_FunctionDef():
    ''' Test the visit_FucntionDef function in cabana_perpart_visitor'''
    aos = Cabana()
    aos.disable_variable_checks()
    set_backend(aos)
    v = cabana_perpart_visitor(aos)
    aos._per_part_visitor = v

    def f1(part1, arg2):
        part1.v = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(f1)))
    x = v.visit(c)
    assert v._part1 == ""
    assert v._config == ""
    correct = '''template < class V >
struct f1_functor{
    config_struct_type _arg2;
    V _v;

    KOKKOS_INLINE_FUNCTION
     f1_functor( V v, config_struct_type arg2) :
    _v(v), _arg2(arg2) {}

    void operator()(const int i, const int a) const{
        _v.access(i, a) = 1;

    }
};
'''
    assert correct == x
