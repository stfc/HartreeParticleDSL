import ast
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
import HartreeParticleDSL.backends.C_AOS.visitors as c_visitors
import re

class cabana_visitor(c_visitors.c_visitor):
    ''' Cabana visitor class. Inherits from the base C_visitor class

    :param parent: The Cabana backend object creating this visitor
    :type parent: :py:class:`HartreeParticleDSL.backends.Cabana_backend.Cabana.Cabana`
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

class cabana_pairwise_visitor(cabana_visitor):
    ''' Cabana pairwise visitor class. Inherits from the cabana_visitor class
    '''
    def __init__(self, parent, indent=4):
        self._part1 = ""
        self._config = ""
        self._part2 = ""

    def check_position(self, node):
        # Disabling check_position for now
#        rval = ""
#        if type(node.value) is ast.Attribute:
#            rval = rval + str(self.visit(node.value))
#            rval = rval + f".{node.attr}"
#        # Search for the core_part.position.[...] in the string.
#        # Upon locating it, we pull out the contents of the [ ]
#        a = re.match("core_part\([a-zA-Z0-9, ]*\).position.[a-zA-Z]*$", rval[rval.find("core_part"):])
#        if a is not None:
#            end = self._parent.get_particle_position(a.group(0)[(a.group(0).find("position.") +9 ):])
#            return rval[0:rval.find(".position")] + end
        return None

    def addSlice(self, slice_name):
        #TODO When adding cabana pairwise support
        pass

    def visit_arguments(self, node):
        raise UnsupportedCodeError("Cabana doesn't yet support pairwise operations")

class cabana_main_visitor(cabana_visitor):
    ''' Cabana main visitor class. '''
    def __init__(self, parent, indent=4):
        self._parent = parent
        self._slices = []
        super().__init__(parent, indent=indent)

    def addSlice(self, slice_name):
        if slice_name not in self._slices:
            self._slices.append(slice_name)

    def visit_FunctionDef(self, node):
        assert node.name is "main"
        assert len(node.args.args) == 0
        self._parent.in_kernel_code = False
        rval = f"int {node.name}( "
        rval = rval + "int argc, char* argv[] )\n{\n"
        self.incrementIndent()
        for a in node.body:
            rval = rval + self.visit(a)
        self.decrementIndent()
        rval = rval + "}\n"
        return rval

class cabana_perpart_visitor(cabana_visitor):
    ''' Cabana per_part visitor class. '''
    def __init__(self, parent, indent=4):
        self._part1 = ""
        self._config = ""
        self._parent = parent
        self._slices = []
        super().__init__(parent, indent=indent)

    def check_position(self, node):
        # Disabling check_position for now
#        rval = ""
#        if type(node.value) is ast.Attribute:
#            rval = rval + str(self.visit(node.value))
#            rval = rval + f".{node.attr}"
#        # Search for the core_part.position.[...] in the string.
#        # Upon locating it, we pull out the contents of the [ ]
#        a = re.match("core_part\([a-zA-Z0-9, ]*\).position.[a-zA-Z]*$", rval[rval.find("core_part"):])
#        if a is not None:
#            end = self._parent.get_particle_position(a.group(0)[(a.group(0).find("position.") +9 ):])
#            return rval[0:rval.find(".position")] + end
        return None

    def resetSlices(self):
        self._slices = []

    def addSlice(self, slice_name):
        if slice_name not in self._slices:
            self._slices.append(slice_name)

    def visit_arguments(self, node):
        if len(node.args) != 2:
            raise IllegalArgumentCountError("Per part function must have 2 arguments for Cabana backend")
        self._part1 = self.visit(node.args[0])
        self._config = self.visit(node.args[1])
        self._parent.variable_scope.add_variable(self.visit(node.args[0]), "FULLPART", False)
        self._parent.variable_scope.add_variable(self.visit(node.args[1]), "FULLCONF", False)
        rval = "const int " + self._part1
        return rval


    def visit_FunctionDef(self, node):
        # Reset slices required for this run
        self._parent.in_kernel_code = True
        self.resetSlices()
        for structure in self._parent._structures:
            self._parent.variable_scope.add_variable(structure, self._parent._structures[structure], False)
        arglist = self.visit(node.args)

        # The main body will be double incremented
        self.incrementIndent()
        self.incrementIndent()
        body = ""
        for a in node.body:
            body = body + self.visit(a)
        self.decrementIndent()
        self.decrementIndent()

        # Find the node name
        nodename = f"{node.name}"

        # Update the slices required for this kernel.
        for slices in self._slices:
            self._parent.add_kernel_slice(f"{nodename}", slices)

        # We template some slices
        rval = "template < "
        all_slices = []
        for slices in self._slices:
            all_slices.append("class " + slices.upper())
        classes = ", ".join(all_slices)
        rval = rval + classes + " >\n"
        rval = rval + f"struct {nodename}_functor" + "{\n"
        self.incrementIndent()
        #FIXME variables
        rval = rval + self.addIndent() + f"config_struct_type _{self._config};" + "\n"
        all_slices = []
        for slices in self._slices:
            rval = rval + self.addIndent() + slices.upper() + " _" + slices + ";\n"
        for structure in self._parent._structures:
            rval = rval + self.addIndent() + self._parent._structures[structure] + " " + structure + ";\n"

        #FIXME Constructor
        rval = rval + "\n"
        rval = rval + self.addIndent() + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + self.addIndent() + f" {nodename}_functor( "
        all_slices = []
        for slices in self._slices:
            all_slices.append(slices.upper() + " " + slices)
        all_slices.append("config_struct_type " + self._config)
        for structure in self._parent._structures:
            all_slices.append(self._parent._structures[structure] + " " + structure.upper())
        classes = ", ".join(all_slices)
        rval = rval + classes + ") :\n"
        all_slices = []
        for slices in self._slices:
            all_slices.append("_" + slices + "(" + slices + ")")
        for structure in self._parent._structures:
            all_slices.append(structure + "(" + structure.upper() + ")")
        classes = ", ".join(all_slices)
        rval = rval + self.addIndent() + classes
        rval = rval + ", _" + f"{self._config}({self._config}) " + "{}\n"

        # If we should update the structures defined in the parents we need a function to
        # enable this
        if len(self._parent._structures) > 0:
            rval = rval + self.addIndent() + "void update_structs("
            all_structs = []
            for struct in self._parent._structures:
                all_structs.append(self._parent._structures[struct] + " " + struct.upper())
            rval = rval + ", ".join(all_structs) + "){\n"
            self.incrementIndent()
            for struct in self._parent._structures:
                rval = rval + self.addIndent() + struct + " = " + struct.upper() + ";\n"
            self.decrementIndent()
            rval = rval + self.addIndent() + "}\n"

        rval = rval + "\n"
        rval = rval + self.addIndent() + "void operator()(const int i, const int a) const{\n"
        rval = rval + body
        rval = rval + "\n" + self.addIndent() + "}\n"

        self.decrementIndent()
        #End of structure
        rval = rval + self.addIndent() + "};\n"

        #Reset part1, config names
        self._part1 = ""
        self._config = ""

        return rval

