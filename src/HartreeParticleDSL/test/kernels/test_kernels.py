import ast
import inspect
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
    parser = main_declaration(main)
    assert isinstance(parser, main_function_wrapper)
