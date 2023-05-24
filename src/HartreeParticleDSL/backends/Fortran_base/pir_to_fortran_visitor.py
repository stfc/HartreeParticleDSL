from __future__ import annotations

import inspect

from HartreeParticleDSL.backends.base_backend.pir_visitor import PIR_Visitor, VisitorError
from psyclone.psyir.backend.fortran import FortranWriter
from HartreeParticleDSL.Particle_IR.nodes.statement import Return, EmptyStatement
from psyclone.psyir.nodes import BinaryOperation, UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.pointer_reference import PointerReference

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str,\
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, ArrayType
from psyclone.psyir.symbols.datatypes import ScalarType, StructureType

class Fortran_PIR_Writer(PIR_Visitor, FortranWriter):

    def __init__(self, indent_string: str="    ", initial_indent_depth: int=0):
        super(PIR_Visitor, self).__init__(indent_string=indent_string, initial_indent_depth=initial_indent_depth)

        self._operator_2_str = {}
        self._operator_2_str[BinaryOperation.Operator.ADD] = "+"
        self._operator_2_str[BinaryOperation.Operator.SUB] = "-"
        self._operator_2_str[BinaryOperation.Operator.MUL] = "*"
        self._operator_2_str[BinaryOperation.Operator.DIV] = "/"
        self._operator_2_str[BinaryOperation.Operator.LT] = "<"
        self._operator_2_str[BinaryOperation.Operator.LE] = "<="
        self._operator_2_str[BinaryOperation.Operator.GT] = ">"
        self._operator_2_str[BinaryOperation.Operator.GE] = ">="
        self._operator_2_str[BinaryOperation.Operator.EQ] = "=="
        self._operator_2_str[BinaryOperation.Operator.AND] = ".and."
        self._operator_2_str[BinaryOperation.Operator.OR] = ".or."

        self._operator_2_str[UnaryOperation.Operator.MINUS] = "-"
        self._operator_2_str[UnaryOperation.Operator.NOT] = ".not."

    @property
    def _nindent(self) -> str:
        '''
        :returns: the current indentation string.
        :rtype: str

        '''
        return self._depth * self._indent

    def indent(self) -> None:
        '''
        Increase the current indentation.
        '''
        self._depth = self._depth + 1

    def dedent(self) -> None:
        '''
        Decrease the current indentation.
        '''
        self._depth = self._depth - 1

    def _visit(self, node : Node) -> str:
        # Make a list of the node's ancestor classes (including
        # itself) in method resolution order (mro), apart from the
        # base "object" class.
        possible_method_names = ["visit_"+curr_class.__name__.lower()+"_node"
                                 for curr_class in inspect.getmro(type(node))]
        possible_method_names.remove("visit_object_node")
        more_possible_method_names = [curr_class.__name__.lower()+"_node" for
                                      curr_class in inspect.getmro(type(node))]
        more_possible_method_names.remove("object_node")
        final_method_names = []
        for i, value in enumerate(possible_method_names):
            final_method_names.append(more_possible_method_names[i])
            final_method_names.append(value)

        # Try to call methods using the class names in the order of
        # the class hierarchy (starting from the current class name).
        for method_name in final_method_names:
            try:
                # pylint: disable=eval-used
                node_result = eval(f"self.{method_name}(node)")

                return node_result

            except AttributeError as excinfo:
                if f"attribute '{method_name}'" in str(excinfo):
                    # This attribute error is because the method that
                    # was tried does not exist.
                    pass
                else:
                    # The method does exist but it has raised an
                    # attribute error so re-raise it here.
                    raise AttributeError(excinfo) from excinfo

        # We haven't found a handler for this node.
        raise VisitorError(
            f"Unsupported node '{type(node).__name__}' found: method names "
            f"attempted were {possible_method_names}.")

    @classmethod
    def get_fortran_datatype(cls, datatype: DataType):
        type_map = {ScalarType.Intrinsic.INTEGER: {ScalarType.Precision.SINGLE: "Integer",
                                                   32: "Integer*4", 64: "Integer*8"},
                    ScalarType.Intrinsic.REAL: {ScalarType.Precision.SINGLE: "Real*4", ScalarType.Precision.DOUBLE: "Real*8"},
                    ScalarType.Intrinsic.CHARACTER: {ScalarType.Precision.UNDEFINED: "Character(len = *)"},
                    ScalarType.Intrinsic.BOOLEAN: {ScalarType.Precision.UNDEFINED: "Logical"}
                   }

        intrinsic = None
        precision = None
        if isinstance(datatype, ScalarType):
            intrinsic = datatype.intrinsic
            precision = datatype.precision
        elif isinstance(datatype, ArrayType):
            intrinsic = datatype._datatype.intrinsic
            precision = datatype._datatype.precision
        type_dict = type_map.get(intrinsic)
        type_str = None
        if type_dict is not None:
            type_str = type_dict.get(precision)
        if type_str is not None:
            return type_str

        for key in type_mapping_str.keys():
            if type_mapping_str[key] == datatype:
                if isinstance(datatype, StructureType):
                    return "Type(" + key + ")"
                else:
                    return key

    def visit_symbol_node(self, symbol: Symbol) -> str:
        # Find the name in the type_mapping_str
        type_str = ""
        for key in type_mapping_str.keys():
            if type_mapping_str[key] == symbol.datatype:
                type_str = key
        return Fortran_PIR_Writer.get_fortran_datatype(symbol.datatype)

    def visit_funcdef_node(self, node: FuncDef) -> str:
        # Don't expect this to appear much in real code, however
        # its used extensively in testing the rest of this code
        # so the implementation here exists.

        return_type = "void"
        returns = node.walk(Return)
        if len(returns) != 0:
            if len(returns[0].children) != 0:
                if isinstance(returns[0].children[0], Literal):
                    datatype = returns[0].children[0].datatype
                    if datatype.intrinsic == ScalarType.Intrinsic.INTEGER:
                        return_type = "int"
                    elif datatype.intrinsic == ScalarType.Intrinsic.REAL:
                        return_type = "double"
                    elif datatype.intrinsic == ScalarType.Intrinsic.BOOLEAN:
                        return_type = "bool"
                    elif datatype.intrinsic == ScalarType.Intrinsic.CHARACTER:
                        return_type = "char*"
                elif isinstance(returns[0].children[0], Reference):
                        assert False
                        return_type = Cabana_PIR_Visitor.get_cpp_datatype(returns[0].children[0].symbol.datatype)
                else:
                    raise NotImplementedError()

        body = self._visit(node.body)

        argument_names = []
        arg_base_names = []
        for arg in node.arguments:
            name = arg.symbol.name
            arg_base_names.append(name)
            if isinstance(arg, PointerReference):
                assert False
            else:
                argument_names.append(f"{name}")

        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        self.indent()
        for arg in node.arguments:
            sym = self._visit(arg.symbol)
            symbol_list = symbol_list + self._nindent + sym + " :: " + arg.symbol.name + "\n"
            
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in (arg_base_names):
                sym = self._visit(symbol)
                symbol_list = symbol_list + self._nindent + sym + " :: " + pairs + "\n"
        self.dedent()

        if return_type == "void":
            rval = f"{self._nindent}Subroutine {node.name}("
            rval = rval + ", ".join(argument_names) + ")\n"
            rval = rval + symbol_list
            rval = rval + body
            rval = rval + self._nindent + f"End Subroutine {node.name}"
        else:
            assert False

        return rval

    def visit_body_node(self, node: Body) -> str:
        self.indent()
        string = ""
        for child in node.children:
            if not isinstance(child, EmptyStatement):
                string = string + self._visit(child)
        self.dedent()
        return string

    def visit_while_node(self, node: While) -> str:
        cond = self._visit(node.condition)
        body = self._visit(node.body)
        string = self._nindent + f"while({cond})" + " do\n"
        string = string + body
        string = string + self._nindent + "enddo\n"
        return string

    def visit_break_node(self, node: Break) -> str:
        return self._nindent + "EXIT\n"

    def visit_ifelseblock_node(self, node: IfElseBlock) -> str:
        cond = self._visit(node.condition)
        ifbody = self._visit(node.ifbody)
        string = self._nindent + f"if({cond}) " + "then\n"
        string = string + ifbody
        if len(node.elsebody.children) > 0:
            string = string + self._nindent + "else"
            if node.is_else_if():
                self.dedent()
                child = self._visit(node.elsebody)
                child = child.lstrip()
                child = child.rstrip()
                self.indent()
                string = string + " " + child
            else:
                string = string + "\n"
                string = string + self._visit(node.elsebody)
        if not node.is_else_if():
            string = string + self._nindent + "endif"
        string = string + "\n"
        return string

    def visit_pointerreference_node(self, node: PointerReference) -> str:
        assert False
