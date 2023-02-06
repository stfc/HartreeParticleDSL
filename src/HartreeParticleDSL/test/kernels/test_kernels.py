import ast
import inspect
from HartreeParticleDSL.HartreeParticleDSL import reset_for_tests
from HartreeParticleDSL.kernel_types.kernels import *

def test_kernel():
    kern = kernel()
    assert kern.get_kernel_tree() is None

def test_main_function_wrapper():
    a = ast.parse("a = a + 1")
    mfw = main_function_wrapper(a)
    assert mfw.get_kernel_tree() is a

def main():
    create_variable(c_int, "a")
    a = a + 1

def test_main_declaration_function():
    reset_for_tests()
    parser = main_declaration(main)
    assert isinstance(parser, main_function_wrapper)

def source_kernel(arg1, arg2):
    create_variable(c_int, "a")
    a = a + 1

def test_source_boundary_kernel_wrapper():
    reset_for_tests()
    source = source_boundary(source_kernel)
    assert isinstance(source, source_boundary_kernel_wrapper)
    assert source.get_source_count() == 0
    source.set_source_count(1000)
    assert source.get_source_count() == 1000

    assert isinstance(source.get_kernel_tree(), ast.Module)
   
def sink_kernel(arg1, arg2):
    create_variable(c_int, "a")
    a = a + 2

def test_sink_boundary_kernel_wrapper():
    reset_for_tests()
    sink = sink_boundary(sink_kernel)
    assert isinstance(sink, sink_boundary_kernel_wrapper)
    assert isinstance(sink.get_kernel_tree(), ast.Module)
