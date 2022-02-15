import ast
from HartreeParticleDSL.backends.base_backend.visitors import baseVisitor
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError
from HartreeParticleDSL.language_utils.variable_scope import variable_scope, variable_access, variable
from HartreeParticleDSL.c_types import c_double, c_int

class c_visitor(baseVisitor):

    _currentIndent = 0
    _else_if = False

    def __init__(self, parent, indent=4):
        self._currentIndent = 0
        self._parent = parent
        self._indent = indent
        super().__init__(parent, indent=indent)

    def incrementIndent(self):
        self._currentIndent = self._currentIndent + self._indent

    def decrementIndent(self):
        self._currentIndent = self._currentIndent - self._indent

    def addIndent(self):
        rval = " " * self._currentIndent
        return rval

    def visit_Str(self, node):
        return f"\"{node.s}\""

    def visit_str(self, node):
        v = self._parent.variable_scope.get_variable(node)
        if v is not None:
            var = variable_access(v)
            return var
        else:
            v = variable(node, "UNKNOWN", False)
            var = variable_access(v)
            return var

    def visit_int(self, node):
        return f"{node}"

    def visit_Add(self, node):
        return " + "

    def visit_Mult(self, node):
        return " * "

    def visit_Sub(self, node):
        return " - "

    def visit_Div(self, node):
        return " / "

    def visit_LtE(self, node):
        return " <= "

    def visit_GtE(self, node):
        return " >= "

    def visit_Lt(self, node):
        return " < "

    def visit_Gt(self, node):
        return " > "

    def visit_USub(self, node):
        return "-"

    def visit_And(self, node):
        return " && "

    def visit_Or(self, node):
        return " || "

    def visit_Not(self, node):
        return " !"

    def visit_BoolOp(self, node):
        code = "( " + str(self.visit(node.values[0]))
        code = code + self.visit(node.op)
        code = code + str(self.visit(node.values[1])) + " )"
        return code

    def visit_UnaryOp(self, node):
        rval = self.visit(node.op) + str(self.visit(node.operand))
        return rval

    def visit_Compare(self, node):
        assert len(node.ops) == 1
        assert len(node.comparators) == 1
        rval = " ( "
        rval = rval + str(self.visit(node.left))
        rval = rval + self.visit(node.ops[0])
        rval = rval + str(self.visit(node.comparators[0]))
        rval = rval + " ) "
        return rval

    def visit_BinOp(self, node):
        rval = "( "
        rval = rval + str(self.visit(node.left))
        rval = rval + self.visit(node.op)
        rval = rval + str(self.visit(node.right))
        rval = rval + " )"
        return rval

    def visit_Name(self, node):
        v = self._parent.variable_scope.get_variable(f"{node.id}")
        if v is not None:
            var = variable_access(v)
            return var
        else:
            v = variable(f"{node.id}", "UNKNOWN", False)
            var = variable_access(v)
            return var

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

    # For pre-python 3.8
    def visit_Num(self, node):
        return f"{node.n}"

    def visit_NameConstant(self, node):
        if str(node.value) == "True":
            return "true"
        elif str(node.value) == "False":
            return "false"

    def visit_Break(self, node):
        return self.addIndent() + "break;\n"

    def visit_Constant(self, node):
        rval = ""
        if type(node.value) is str:
            rval = rval + f"\"{node.value}\""
        else:
            rval = rval + f"{node.value}"
        return rval

    def visit_Assign(self, node):
        assert len(node.targets) == 1
        rval = self.addIndent()
        rval = rval + str(self.visit(node.targets[0]))
        rval = rval + " = "
        rval = rval + str(self.visit(node.value)).lstrip()
        rval = rval + ";\n"
        return rval

    def visit_arg(self, node):
        return f"{node.arg}"

    def visit_arguments(self, node):
        rval = ""
        for child in node.args:
            rval = rval + "struct part *"
            rval = rval + self.visit(child)
            if child is not node.args[-1]:
                rval = rval + ", "
        return rval

    def visit_FunctionDef(self, node):
        rval = f"void {node.name}( "
        rval = rval + self.visit(node.args)
        rval = rval + " )\n{\n"
        self.incrementIndent()
        for a in node.body:
            rval = rval + self.visit(a)
        self.decrementIndent()
        rval = rval + "}\n"
        return rval

    def visit_If(self, node):
        rval = ""
        if not self._else_if:
            rval = rval + self.addIndent()
        self.elseIf = False
        rval = rval + "if( "
        rval = rval + str(self.visit(node.test))
        rval = rval + " ){\n"
        self.incrementIndent()
        for child in node.body:
            rval = rval + self.visit(child)
        self.decrementIndent()
        rval = rval + self.addIndent()
        rval = rval + "}"
        for child in node.orelse:
            if type(child) == ast.If:
                rval = rval + "else "
                self._else_if = True
                rval = rval + self.visit(child)
            else:
                rval = rval + "else{\n"
                self.incrementIndent()
                rval = rval + self.visit(child)
                self.decrementIndent()
                rval = rval + self.addIndent()
                rval = rval + "}"
        rval = rval + "\n"
        return rval

    def visit_Call(self, node):
        from HartreeParticleDSL.HartreeParticleDSL import gen_invoke
        # Disable variable checking during visit call as they might be python variables
        self._parent.disable_variable_checks()
        # Recursively build up the function name if it is of style mod1.mod2.func
        function_name = ""
        if isinstance(node.func, ast.Name):
            function_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            temp_node = node.func
            while isinstance(temp_node, ast.Attribute):
                function_name = temp_node.attr + "." + function_name
                temp_node = temp_node.value
            function_name = temp_node.id + "." + function_name
            #Remove the . at the end
            function_name = function_name[:-1]
        func_string = f"{function_name}( "
        arguments = []
        # TODO #2L: Temporary version until visitors return strings instead of printing
        for child in node.args:
            vis = str(self.visit(child))
            arguments.append(vis)
        for child in node.keywords:
            string = f"{child.arg}={str(self.visit(child.value))}"
            arguments.append(string)
        rval = ""
        if function_name != "invoke":
            rval = self._parent.call_language_function(function_name, *arguments, current_indent=self._currentIndent,
                                                       indent=self._indent)
        elif function_name == "invoke":
            for child in node.args:
                rval = gen_invoke(child.id, self._currentIndent, self._indent)
        self._parent.enable_variable_checks()
        return rval

    def visit_Module(self, node):
        rval = ""
        for a in ast.iter_child_nodes(node):
            rval = rval + self.visit(a)
        return rval

    def visit_Index(self, node):
        return self.visit(node.value)

    def visit_Subscript(self, node):
        rval = self.visit(node.value)
        # node.slice must be an ast.Index right now
        # Add the array index to the lowest child
        last_child = rval
        while last_child.child is not None:
            last_child = last_child.child
        last_child.add_array_index(self.visit(node.slice))
#        rval = rval + "["
#        rval = rval + self.visit(node.slice)
#        rval = rval + "]"
        return rval

    def visit_For(self, node):
        if len(node.orelse) != 0:
            raise IllegalLoopError(f"Else clauses on Loops are not supported in {self._parent.__class__.__name__}")
        assert len(node.orelse) == 0
        if type(node.iter) != ast.Call:
            raise IllegalLoopError(f"Only range loops are supported in {self._parent.__class__.__name__}")
        if type(node.iter) == ast.Call:
            if node.iter.func.id != "range":
                raise IllegalLoopError(f"Only range loops are supported in {self._parent.__class__.__name__}")
        startval = 0
        endval = node.iter.args[0]
        increment = 1
        if len(node.iter.args) == 2:
            startval = node.iter.args[0]
            endval = node.iter.args[1]
        if len(node.iter.args) == 3:
            startval = node.iter.args[0]
            endval = node.iter.args[1]
            increment = node.iter.args[2]
            if type(increment) is ast.UnaryOp:
                if type(increment.op) is ast.USub:
                    increment = -1 *  increment.operand.n
            else:
                increment = increment.n
       
        # Add the loop variable to the scope
        self._parent.variable_scope.add_variable(node.target.id, c_int, False)

        rval = ""
        rval = rval + self.addIndent()
        rval = rval + "for( int "
        rval = rval + str(self.visit(node.target))
        rval = rval + " = "
        rval = rval + self.visit(startval)
        rval = rval + "; "
        rval = rval + str(self.visit(node.target))
        if increment > 0:
            rval = rval + " < "
        else:
            rval = rval + " >= "
        rval = rval + self.visit(endval)
        rval = rval + "; "
        rval = rval + str(self.visit(node.target))
        rval = rval + " = "
        rval = rval + str(self.visit(node.target))
        rval = rval + " + "
        rval = rval + self.visit(increment)
        rval = rval + ")\n"
        rval = rval + self.addIndent()
        rval = rval + "{\n"
        self.incrementIndent()
        for child in node.body:
            rval = rval + self.visit(child)
        self.decrementIndent()
        rval = rval + self.addIndent()
        rval = rval + "}\n"
        self._parent.variable_scope.remove_variable(node.target.id)
        return rval

    def visit_While(self, node):
        rval = self.addIndent()
        rval = rval + "while("
        rval = rval + self.visit(node.test)
        rval = rval + "){\n"
        self.incrementIndent()
        for child in node.body:
            rval = rval + self.visit(child)
        self.decrementIndent()
        rval = rval + self.addIndent()
        rval = rval + "}\n"
        return rval

    def visit_Expr(self, node):
        rval = ""
        for a in ast.iter_child_nodes(node):
            rval = rval + str(self.visit(a))
#            if type(a) is ast.Call:
#                rval = rval + ";\n"
        return rval

    def generic_visit(self, node):
        raise UnsupportedCodeError(f"Found unsupported node of type {type(node)}")

class c_pairwise_visitor(c_visitor):
    def visit_arguments(self, node):
        if len(node.args) != 4:
            raise IllegalArgumentCountError("Pairwise function must have 4 arguments for C_AOS backend")
        self._parent.variable_scope.add_variable(self.visit(node.args[0]), "struct part", True)
        self._parent.variable_scope.add_variable(self.visit(node.args[1]), "struct part", True)
        self._parent.variable_scope.add_variable(self.visit(node.args[2]), c_double, False)
        self._parent.variable_scope.add_variable(self.visit(node.args[3]), "struct config_type", True)
        rval = "struct part *"
        rval = rval + self.visit(node.args[0])
        rval = rval + ", struct part *"
        rval = rval + self.visit(node.args[1])
        rval = rval + ", double "
        rval = rval + self.visit(node.args[2])
        rval = rval + ", struct config_type *"
        rval = rval + self.visit(node.args[3])
        return rval

class c_main_visitor(c_visitor):
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

class c_perpart_visitor(c_visitor):
    def visit_arguments(self, node):
        if len(node.args) != 2:
            raise IllegalArgumentCountError("Per part function must have 2 arguments for C_AOS backend")
        self._parent.variable_scope.add_variable(self.visit(node.args[0]), "struct part", True)
        self._parent.variable_scope.add_variable(self.visit(node.args[1]), "struct config_type", True)

        rval = "struct part *"
        rval = rval + self.visit(node.args[0])
        rval = rval + ", struct config_type *"
        rval = rval + self.visit(node.args[1])  
        return rval
