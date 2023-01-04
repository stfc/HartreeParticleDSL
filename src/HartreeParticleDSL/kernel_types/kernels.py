import ast
import inspect
from functools import wraps

class kernel():
    def __init__(self, kernel_tree=None):
        pass
    def get_kernel_tree(self):
        '''
        Returns the kernel tree held by this object.

        :returns: The Kernel tree held by this object
        :rtype: :py:class:`ast.Module`
        '''
        pass


class pairwise_kernel_wrapper(kernel):
    '''
    Pairwise Kernel wrapper class used by the DSL to store
    user-written code as an AST that can be used by backends
    to produce output.

    :param kernel_tree: The AST parsed tree produced for a pairwise kernel.
    :type kernel_tree: :py:class:`ast.Module`
    '''
    def __init__(self, kernel_tree):
        self._pairwise_kernel_function = kernel_tree

    def get_kernel_tree(self):
        '''
        Returns the kernel tree held by this object.

        :returns: The Kernel tree held by this object
        :rtype: :py:class:`ast.Module`
        '''
        return self._pairwise_kernel_function

class perpart_kernel_wrapper(kernel):
    '''
    Perpart Kernel wrapper class used by the DSL to store
    user-written code as an AST that can be used by backends
    to produce output.

    :param kernel_tree: The AST parsed tree produced for a per-part kernel.
    :type kernel_tree: :py:class:`ast.Module`
    '''
    def __init__(self, kernel_tree):
        self._perpart_kernel_function = kernel_tree

    def get_kernel_tree(self):
        '''
        Returns the kernel tree held by this object.

        :returns: The Kernel tree held by this object
        :rtype: :py:class:`ast.Module`
        '''
        return self._perpart_kernel_function

class main_function_wrapper():
    '''
    Main function wrapper class used by the DSL to store
    user-written code as an AST that can be used by backends
    to produce output.

    :param kernel_tree: The AST parsed tree produced for a per-part kernel.
    :type kernel_tree: :py:class:`ast.Module`
    '''
    def __init__(self, kernel_tree):
        self._main_function = kernel_tree

    def get_kernel_tree(self):
        '''
        Returns the function tree held by this object.

        :returns: The Function tree held by this object
        :rtype: :py:class:`ast.Module`
        '''
        return self._main_function


def pairwise_interaction(function):
    '''
    Wrapper decorator for pairwise kernel functions.
    The function takes the kernel function, and generates the
    AST, and then creates the pairwise kernel wrapper from that.
    Finally, it registers the kernel with the DSL.

    Example use:

    >>> @pairwise_interaction
    >>> def my_kernel(...):
    >>>    ....

    :returns: The pairwise kernel wrapper corresponding to the
              given kernel.
    :rtype: :py:class:`HartreeParticleDSL.kernel_types.kernels.pairwise_kernel_wrapper`
    '''
    @wraps(function)
    def wrapper():
        import HartreeParticleDSL.HartreeParticleDSL as HPDSL
        tree = ast.parse(inspect.getsource(function))
        parser = pairwise_kernel_wrapper(tree)
        HPDSL.register_kernel(tree.body[0].name, parser)
        return parser
    return wrapper()

def perpart_interaction(function):
    '''
    Wrapper decorator for perpart kernel fucntions.
    The function takes the kernel function, and generates the
    AST, and then creates the perpart kernel wrapper from that.
    Finally, it registers the kernel with the DSL.

    Example use:
    >>> @perpart_interaction
    >>> def my_kernel(...):
    >>>    ....

    :returns: The perpart kernel wrapper corresponding to the
              given kernel.
    :rtype: :py:class:`HartreeParticleDSL.kernel_types.kernels.perpart_kernel_wrapper`
    '''
    @wraps(function)
    def wrapper():
        import HartreeParticleDSL.HartreeParticleDSL as HPDSL
        tree = ast.parse(inspect.getsource(function))
        parser = perpart_kernel_wrapper(tree)
        HPDSL.register_kernel(tree.body[0].name, parser)
        return parser
    return wrapper()

def main_declaration(function):
    '''
    Wrapper decorator for the main function (Algorithm).
    The function takes the main function, generates the ast and 
    main_function_wrapper, and then calls the DSL to print the
    main function according to the chosen backend. This is currently
    output to stdout.

    Example use:
    >>> @main_declaration
    >>> def main(...):
    >>>    ....

    :returns: The main function wrapper corresponding to the 
              given function.
    :rtype: :py:class:`HartreeParticleDSL.kernel_types.kernels.main_function_wrapper`
    '''
    @wraps(function)
    def wrapper():
        import HartreeParticleDSL.HartreeParticleDSL as HPDSL
        tree = ast.parse(inspect.getsource(function))
        parser = main_function_wrapper(tree)
        HPDSL.print_main(tree)
        return parser
    return wrapper()
