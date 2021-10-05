import ast
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
import HartreeParticleDSL.backends.C_AOS.visitors as c_visitors

class fdps_visitor(c_visitors.c_visitor):
    ''' FPDS visitor class. Inherits from the base C_visitor class

    :param parent: The FDPS backend object creating this visitor
    :type parent: :py:class:`HartreeParticleDSL.backends.FDPS_backend.FDPS.FDPS`
    :param int indent: The indentation level to use when indentations are needed.
                       Default value is 4
    '''
    def __init__(self, parent, indent=4):
        super().__init__(parent, indent=indent)

    def visit_Attribute(self, node):
        rval = ""
        if type(node.value) is ast.Attribute:
            rval = rval + self.visit(node.value)
            rval = rval + f".{node.attr}"
        else:
            rval = rval + self.visit(node.value)
            rval = rval + "."
            rval = rval + self.visit(node.attr)
        return rval

class fdps_pairwise_visitor(fdps_visitor):
    ''' FDPS pairwise visitor class. Inherits from the fdps_visitor class
    '''

    def visit_arguments(self, node):
        raise UnsupportedCodeError("FDPS doesn't yet support pairwise operations")

class fdps_main_visitor(fdps_visitor):
    ''' FDPS main visitor class. '''
    def __init__(self, parent, indent=4):
        self._parent = parent
        super().__init__(parent, indent=indent)

    def visit_Expr(self, node):
        rval = ""
        for a in ast.iter_child_nodes(node):
            rval = rval + self.visit(a)
        return rval

    def visit_FunctionDef(self, node):
        assert node.name is "main"
        rval = f"int {node.name}( "
        rval = rval + self.visit(node.args)
        rval = rval + " )\n{\n"
        self.incrementIndent()
        for a in node.body:
            rval = rval + self.visit(a)
        self.decrementIndent()
        rval = rval + "}\n"
        return rval

class fdps_perpart_visitor(fdps_visitor):
    ''' FDPS per_part visitor class. '''
    def visit_arguments(self, node):
        if len(node.args) != 2:
            raise IllegalArgumentCountError("Per part function must have 2 arguments for FDPS backend")
        rval = "FullParticle &" + self.visit(node.args[0])
        # Assuming config type is a class for C++
        rval = rval + ", config &" + self.visit(node.args[1])
        return rval
