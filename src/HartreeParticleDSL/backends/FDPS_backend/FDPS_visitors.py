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
        x = self.check_position(node)
        if x is not None:
            return x
        var = self.visit(node.value)
        attr = self.visit(node.attr)
        childless = var
        while childless.child is not None:
            childless = childless.child
        childless.child = attr
        return var

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
        self._parent.variable_scope.add_variable(self.visit(node.args[0]), "FullParticle", False)
        self._parent.variable_scope.add_variable(self.visit(node.args[1]), "config_type", False)
        rval = "FullParticle& " + self.visit(node.args[0])
        # Assuming config type is a class for C++
        rval = rval + ", config_type& " + self.visit(node.args[1])
        return rval
