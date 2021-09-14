import ast
from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError, \
                                                            IllegalArgumentCountError

class c_visitor(ast.NodeVisitor):

    _currentIndent = 0
    _else_if = False

    def __init__(self, indent=4):
        self._currentIndent = 0
        self._indent = indent

    def incrementIndent(self):
        self._currentIndent = self._currentIndent + self._indent

    def decrementIndent(self):
        self._currentIndent = self._currentIndent - self._indent

    def addIndent(self):
        for i in range(self._currentIndent):
            print(" ", end="")

    def visit_Str(self, node):
        print(f"\"{node.s}\"", end="")

    def visit_str(self, node):
        print(node, end="")

    def visit_int(self, node):
        print(node, end="")

    def visit_Add(self, node):
        print(" + ", end = "")

    def visit_Mult(self, node):
        print(" * ", end = "")

    def visit_Sub(self, node):
        print(" - ", end = "")

    def visit_LtE(self, node):
        print(" <= ", end="")

    def visit_GtE(self, node):
        print(" >= ", end="")

    def visit_Lt(self, node):
        print(" < ", end="")

    def visit_Gt(self, node):
        print(" > ", end="")

    def visit_USub(self, node):
        print(" -", end="")

    def visit_UnaryOp(self, node):
        self.visit(node.op)
        self.visit(node.operand)

    def visit_Compare(self, node):
        assert len(node.ops) == 1
        assert len(node.comparators) == 1
        print(" ( ", end="")
        self.visit(node.left)
        self.visit(node.ops[0])
        self.visit(node.comparators[0])
        print(" ) ", end="")

    def visit_BinOp(self, node):
        print("( ", end = "")
        self.visit(node.left)
        self.visit(node.op)
        self.visit(node.right)
        print(" )", end = "")

    def visit_Name(self, node):
        print(node.id, end="")

    def visit_Attribute(self, node):
        if type(node.value) is ast.Attribute:
            self.visit(node.value)
            print(f".{node.attr}", end="")
        elif type(node.value) is ast.Subscript:
            self.visit(node.value)
        elif node.value.id.startswith("var_"):
            print(node.value.id.replace("var_", ""), end="")
            print(" ", end="")
            self.visit(node.attr)
            pass
        else:
            self.visit(node.value)
            print("->", end="")
            self.visit(node.attr)

    # For pre-python 3.8
    def visit_Num(self, node):
        print(node.n, end = "")

    def visit_Constant(self, node):
        print(node.value, end = "")

    def visit_Assign(self, node):
        assert len(node.targets) == 1
        self.addIndent()
        self.visit(node.targets[0])
        print(" = ", end="")
        self.visit(node.value)
        print(";")

    def visit_arg(self, node):
        print(node.arg, end="")

    def visit_arguments(self, node):
        for child in node.args:
            print("struct part *", end="")
            self.visit(child)
            if child is not node.args[-1]:
                print(", ", end="")

    def visit_FunctionDef(self, node):
        print(f"void {node.name}( ", end="")
        self.visit(node.args)
        print(" )\n{\n",end="")
        self.incrementIndent()
        for a in node.body:
            self.visit(a)
        self.decrementIndent()
        print("}\n")


    def visit_If(self, node):
        if not self._else_if:
            self.addIndent()
        self.elseIf = False
        print("if( ", end="")
        self.visit(node.test)
        print(" ){\n", end="")
        self.incrementIndent()
        for child in node.body:
            self.visit(child)
        self.decrementIndent()
        self.addIndent()
        print("}", end="")
        for child in node.orelse:
            if type(child) == ast.If:
                print("else ", end="")
                self._else_if = True
                self.visit(child)
            else:
                print("else{\n", end="")
                self.incrementIndent()
                self.visit(child)
                self.decrementIndent()
                self.addIndent()
                print("}", end="")
        print("\n", end="")

    def visit_Call(self, node):
        self.visit(node.func)
        print("( ", end="")
        for child in node.args:
            self.visit(child)
        print(" )", end="")
        for child in node.keywords:
            self.visit(child)

    def visit_Module(self, node):
        for a in ast.iter_child_nodes(node):
            self.visit(a)

    def visit_Index(self, node):
#        print("[", end="")
        self.visit(node.value)
#        print("]", end="")

    def visit_Subscript(self, node):
        self.visit(node.value)
        # node.slice must be an ast.Index right now
        print("[", end="")
        self.visit(node.slice)
        print("]", end="")

    def visit_For(self, node):
        if len(node.orelse) != 0:
            raise IllegalLoopError("Else clauses on Loops are not supported in C_AOS")
        assert len(node.orelse) == 0
        if type(node.iter) != ast.Call:
            raise IllegalLoopError("Only range loops are supported in C_AOS")
        if type(node.iter) == ast.Call:
            if node.iter.func.id != "range":
                raise IllegalLoopError("Only range loops are supported in C_AOS")
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
            
        self.addIndent()
        print("for( int ", end="")
        self.visit(node.target)
        print(" = ", end="")
        self.visit(startval)
        print("; ", end="")
        self.visit(node.target)
        if increment > 0:
            print(" < ", end = "")
        else:
            print(" >= ", end = "")
        self.visit(endval)
        print("; ", end="")
        self.visit(node.target)
        print(" = ", end="")
        self.visit(node.target)
        print(" + ", end="")
        self.visit(increment)
        print(")\n", end="")
        self.addIndent()
        print("{\n", end="")
        self.incrementIndent()
        for child in node.body:
            self.visit(child)
        self.decrementIndent()
        self.addIndent()
        print("}\n", end="")

    def visit_While(self, node):
        self.addIndent()
        print("while(", end="")
        self.visit(node.test)
        print("){")
        self.incrementIndent()
        for child in node.body:
            self.visit(child)
        self.decrementIndent()
        self.addIndent()
        print("}\n", end="")

    def visit_Expr(self, node):
        for a in ast.iter_child_nodes(node):
            self.visit(a)
        print(";")

    def generic_visit(self, node):
        raise UnsupportedCodeError(f"Found unsupported node of type {type(node)}")

class c_pairwise_visitor(c_visitor):
    def visit_arguments(self, node):
        if len(node.args) != 4:
            raise IllegalArgumentCountError("Pairwise function must have 4 arguments for C_AOS backend")
        print("struct part *", end = "")
        self.visit(node.args[0])
        print(", struct part *", end = "")
        self.visit(node.args[1])
        print(", double ", end = "")
        self.visit(node.args[2])
        print(", struct config_type *", end="")
        self.visit(node.args[3])

class c_main_visitor(c_visitor):

    def visit_Call(self, node):
        from HartreeParticleDSL.HartreeParticleDSL import initialise, gen_invoke, println
        ##Function calls must have a Name
        if not isinstance(node.func, ast.Name):
            print(f"The node has fields {node._fields}")
            print(type(node.func))
            print(f"func fields {node.func._fields}")
        assert isinstance(node.func, ast.Name)
        function_name = node.func.id
        if function_name == "initialise":
            particle_count = 0
            filename = ""
            for child in node.keywords:
                if child.arg == "particle_count":
                    particle_count = child.value.n
                elif child.arg == "filename":
                    filename = child.value.s
            string = initialise(particle_count=particle_count, filename=filename)
            self.addIndent()
            print(f"struct config_type* config = malloc(sizeof(struct config_type));")
            self.addIndent()
            print(f"struct part* parts = {string}")
        elif function_name == "cleanup":
            self.addIndent()
            print("free(config);");
            self.addIndent()
            print("free(parts);")
        elif function_name == "invoke":
            for child in node.args:
                gen_invoke(child.id, self._currentIndent, self._indent)
        elif function_name == "println":
            s = ast.Expression()
            s.body = node
            code = compile(s, '<string>', 'eval')
            string = eval(code)
            self.addIndent()
            print(string, end="")
        else:
            self.visit(node.func)
            print("(", end="")
            for child in node.args:
                self.visit(child)
            print(")", end="")
            for child in node.keywords:
                self.visit(child)

    def visit_FunctionDef(self, node):
        assert node.name is "main"
        print(f"int {node.name}( ", end="")
        self.visit(node.args)
        print(" )\n{\n",end="")
        self.incrementIndent()
        for a in node.body:
            self.visit(a)
        self.decrementIndent()
        print("}\n")

class c_perpart_visitor(c_visitor):
    def visit_arguments(self, node):
        if len(node.args) != 2:
            raise IllegalArgumentCountError("Per part function must have 2 arguments for C_AOS backend")
        print("struct part *", end="")
        self.visit(node.args[0])
        print(", struct config_type *", end="")
        self.visit(node.args[1])  
