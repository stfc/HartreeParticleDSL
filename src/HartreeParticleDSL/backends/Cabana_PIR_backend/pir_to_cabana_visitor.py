from __future__ import annotations

from HartreeParticleDSL.HartreeParticleDSL import get_mpi

from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError

from HartreeParticleDSL.Particle_IR.nodes.array_mixin import ArrayMixin
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.funcdef import FuncDef
from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.node import Node
from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.nodes.statement import Return, EmptyStatement
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.particle_position_reference import \
        ParticlePositionReference
from HartreeParticleDSL.Particle_IR.nodes.particle_reference import ParticleReference
from HartreeParticleDSL.Particle_IR.nodes.pointer_reference import PointerReference
from HartreeParticleDSL.backends.base_backend.pir_visitor import PIR_Visitor

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, ScalarType, \
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, StructureType, ArrayType

class Cabana_PIR_Visitor(PIR_Visitor):
    '''
    Generic Cabana PIR visitor.

    :param str indent_string: The indentation to use in this Visitor. Default \
                              "    ".
    :param int initial_indent_depth: Specifies the initial indentation, default \
                                     is 0.

    '''

    def __init__(self, parent, **kwargs):
        super().__init__(**kwargs)
        self._in_kernel = False
        self._slices = []
        self._parent = parent

    @classmethod
    def get_cpp_datatype(cls, datatype: DataType):
        type_map = {INT_TYPE: "int",
                FLOAT_TYPE: "float",
                DOUBLE_TYPE: "double",
                INT64_TYPE: "int64_t",
                INT32_TYPE: "int32_t",
                BOOL_TYPE: "bool",
                STRING_TYPE: "char*",
                BASE_PARTICLE_TYPE: "part"}
        t = type_map.get(datatype)
        if isinstance(datatype, ArrayType):
            childtype = datatype._datatype
            t = type_map.get(childtype)
        if t is not None:
            return t
        for key in type_mapping_str.keys():
            if isinstance(datatype, StructureType):
                if type_mapping_str[key] == datatype:
                    return "struct " + key
            if type_mapping_str[key] == datatype:
                return key

    def addSlice(self, slice_name: str):
        if slice_name not in self._slices:
            self._slices.append(slice_name)

    def visit_symbol_node(self, symbol: Symbol) -> str:
        # Find the name in the type_mapping_str
        type_str = ""
        for key in type_mapping_str.keys():
            if type_mapping_str[key] == symbol.datatype:
                type_str = key
        return Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)

    def visit_pointersymbol_node(self, symbol: PointerSymbol) -> str:
        # Find the name in the type_mapping_str
        type_str = ""
        type_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype._datatype)
        return f"{type_str}*"

    def visit_arraysymbol_node(self, symbol: ArraySymbol) -> str:
        # Find the name in the type_mapping_str
        type_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)
        all_dynamic = True
        for index, extent in enumerate(symbol.datatype.shape):
            if isinstance(extent, int) and (index == 0 or all_dynamic == False):
                all_dynamic = False
            elif isinstance(extent, int) or not all_dynamic:
                raise UnsupportedCodeError("Got an ArraySymbol with a mixed "
                                           "dynamic and fixed typing, which is "
                                           "not supported.")
        if all_dynamic:
            type_str = type_str + ("*"*len(symbol.datatype.shape))
        else:
            for extent in symbol.datatype.shape:
                type_str = type_str + "[" + f"{extent}" + "]"
        return type_str

    # Implement visitors
    def visit_pairwisekernel_node(self, node: PairwiseKernel) -> str:
        self._in_kernel = True
        self._slices = []
        raise NotImplementedError("Pairwise kernels not yet supported in cabana.")
        # self._parent.register_kernel(node.name, node)

    def visit_perpartkernel_node(self, node: PerPartKernel) -> str:
        self._in_kernel = True
        self._slices = []
        # Check rules. Doesn't contain invokes.
        if len(node.walk(Invoke)) != 0:
            raise UnsupportedCodeError("Per part kernel cannot contain Invokes.")

        # Arg 1 is always a particle
        #TODO Save the particle 1 symbol name.
        self._parent.register_kernel(node.name, node)
        argument_names = []
        for arg in node.arguments:
            # Need the correct typing on these arguments
            argument_names.append(arg.symbol.name)

        self.indent()
        # Extra indentation for symbol table
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in argument_names:
                sym = self._visit(symbol)
                symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
        self.dedent()

        #Current indentation level isn't quite correct
        kernel_body = self._visit(node.body)
        self.dedent()

        # Start codegen.

        # Create the templated functor
        rval = self._nindent + "template < "
        all_slices = []
        for slices in self._slices:
            all_slices.append("class " + slices.upper())
        classes = ", ".join(all_slices)
        rval = rval + classes + " >\n"
        rval = rval + f"{self._nindent}struct {node.name}_functor" + "{\n"
        self.indent()
        rval = rval + f"{self._nindent}config_struct_type _{node.arguments[1].symbol.name};\n"
        all_slices = []
        for slices in self._slices:
            rval = rval + f"{self._nindent}{slices.upper()} _{slices};\n"
        for structure in self._parent.structures:
            rval = rval + f"{self._nindent}{structure} _{structure};\n"

        # Constructor
        rval = rval + "\n"
        rval = rval + f"{self._nindent}KOKKOS_INLINE_FUNCTION\n"
        rval = rval + f"{self._nindent} {node.name}_functor( "
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"{slices.upper()} {slices}")
        all_slices.append("config_struct_type " + node.arguments[1].symbol.name)
        for structure in self._parent.structures:
            all_slices.append(f"{structure} {structure.upper()}")
        classes = ", ".join(all_slices)
        rval = rval + classes + "):\n"
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"_{slices}({slices})")
        for structure in self._parent.structures:
            all_slices.append(f"_{structure}({structure.upper()})")
        classes = ", ".join(all_slices)
        rval = rval + f"{self._nindent}{classes}"
        rval = rval + ", " + f"_{node.arguments[1].symbol.name}({node.arguments[1].symbol.name})" + "{}\n"

        # Need to have an update_structs call if there are structures
        if len(self._parent.structures) > 0:
            rval = rval + f"\n{self._nindent}void update_structs("
            all_structs = []
            for struct in self._parent.structures:
                all_structs.append(struct + " " + struct.upper())
            rval = rval + ", ".join(all_structs) + "){\n"
            self.indent()
            for struct in self._parent.structures:
                rval = rval + self._nindent + "_" + struct + " = " + struct.upper() + ";\n"
            self.dedent()
            rval = rval + self._nindent + "}\n"

        rval = rval + "\n"
        rval = rval + self._nindent + "void operator()(const int i, const int a) const{\n"
        rval = rval + symbol_list
        rval = rval + kernel_body
        rval = rval + self._nindent + "}\n"

        self.dedent()
        rval = rval + self._nindent + "};\n"
        for slices in self._slices:
            self._parent.add_kernel_slices(node.name, slices)
        return rval

    def visit_mainkernel_node(self, node: MainKernel) -> str:
        self._in_kernel = False

        #Current indentation level isn't quite correct
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if symbol.name == "config":
                continue
            sym = self._visit(symbol)
            symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"

        self.dedent()
        body = self._visit(node.body)
        rval = self._nindent + f"int {node.name}( int argc, char* argv[] )" + "{\n"
        rval = rval + symbol_list
        # Append the code body
        rval = rval + body
        rval = rval + self._nindent + "}\n"

        return rval

    def visit_body_node(self, node: Body) -> str:
        self.indent()
        string = ""
        for child in node.children:
            if not isinstance(child, EmptyStatement):
                string = string + self._nindent + self._visit(child) + "\n"
        self.dedent()
        return string

    def visit_ifelseblock_node(self, node: IfElseBlock) -> str:
        cond = self._visit(node.condition)
        ifbody = self._visit(node.ifbody)
#        if node.is_else_if():
#            raise NotImplementedError()
        string = f"if({cond})" + "{\n"
        string = string + ifbody
        string = string + self._nindent + "}"
        if len(node.elsebody.children) > 0:
            string = string + "else"
            if node.is_else_if():
                self.dedent()
                child = self._visit(node.elsebody)
                child = child.lstrip()
                child = child.rstrip()
                self.indent()
                string = string + " " + child
            else:
                string = string + "{\n"
                string = string + self._visit(node.elsebody)
                string = string + self._nindent + "}"
        return string

    def visit_invoke_node(self, node: Invoke) -> str:
        rval = ""
        for invoke in node.children:
            # TODO Check if the call updates particle position.
            name = invoke.value
            pir_kernel = self._parent._kernels[name]
            assigns = pir_kernel.walk(Assignment)
            updates_part_pos = False
            for assign in assigns:
                if isinstance(assign.lhs, ParticlePositionReference):
                    updates_part_pos = True
                # TODO If it contains a "Call" node to an unknown call maybe
                # this should just be set to True as well?
            rval = rval + f"Kokkos::deep_copy(config.config, config.config_host);\n"
            if len(self._parent.structures) > 0:
                rval = rval + self._nindent
                rval = rval + f"{invoke.value}.update_structs("
                struct_list = []
                for struct in self._parent.structures:
                    struct_list.append(struct)
                rval = rval + ", ".join(struct_list) + ");\n"
            rval = rval + f"{self._nindent}Cabana::simd_parallel_for(simd_policy, {invoke.value}, "
            rval = rval + "\"" + invoke.value + "\");\n"
            # TODO Check if we need to block (i.e. only if kernel dependencies or final kernel)
            rval = rval + self._nindent + "Kokkos::fence();"
            if updates_part_pos and (self._parent.boundary_condition is not None):
                bound = self._parent.boundary_condition
                if len(self._parent.structures) > 0:
                    rval = rval + self._nindent
                    rval = rval + f"{invoke.value}.update_structs("
                    struct_list = []
                    for struct in self._parent.structures:
                        struct_list.append(struct)
                    rval = rval + ", ".join(struct_list) + ");\n"
                rval = rval + f"{self._nindent}Cabana::simd_parallel_for(simd_policy, {bound.name}"
                rval = rval + "\"" + bound.name + "\");\n"
                rval = rval + self._nindent + "Kokkos::fence();"
                # What if MPI?
                if get_mpi():
                    pass
                    assert False
                    # TODO  Fix this for general purpose
                    #rval = rval + "migrator.exchange_data( particle_aosoa, neighbors, myrank, particle_aosoa.size(), dx*10.0, movement_since_remesh, host_field(0).x_min_local, host_field(0).x_max_local);"
        return rval

    def visit_arraymember_node(self, node: ArrayMember) -> str:
        indices = []
        for index in node.indices:
            indices.append(f"[{self._visit(index)}]")
        return node.name + "".join(indices)

    def visit_arrayreference_node(self, node: ArrayReference) -> str:
        indices = []
        for index in node.indices:
            indices.append(f"[{self._visit(index)}]")
        return node.symbol.name + "".join(indices)

    def visit_assignment_node(self, node: Assignment) -> str:
        lhs = self._visit(node.lhs)
        rhs = self._visit(node.rhs)
        return f"{lhs} = {rhs};"

    def visit_call_node(self, node: Call) -> str:
        # Need to handle special function calls
        # Use the call in the parent
        func_name = node.func_name
        args = []
        for arg in node.children:
            args.append(self._visit(arg))
        try:
            current_indent = len(self._nindent)
            indent = len(self._indent)
            return self._parent.call_language_function(func_name, *args,
                    current_indent=current_indent, indent=indent)
        except:
            pass
        arg_string = ", ".join(args)
        end = ";"
        return f"{func_name}({arg_string}){end}"

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
                    elif datatype.intrinsic == ScalarType.Intrinsic.FLOAT:
                        return_type = "double"
                    elif datatype.intrinsic == ScalarType.Intrinsic.BOOLEAN:
                        return_type = "bool"
                    elif datatype.intrinsic == ScalarType.Intrinsic.CHARACTER:
                        return_type = "char*"
                elif isinstance(returns[0].children[0], Reference):
                        return_type = Cabana_PIR_Visitor.get_cpp_datatype(returns[0].children[0].symbol.datatype)
                else:
                    raise NotImplementedError()

        body = self._visit(node.body)

        argument_names = []
        arg_base_names = []
        for arg in node.arguments:
            name = arg.symbol.name
            arg_base_names.append(name)
            typ = Cabana_PIR_Visitor.get_cpp_datatype(arg.symbol.datatype)
            if isinstance(arg, PointerReference):
                argument_names.append(f"*{typ} {name}")
            else:
                argument_names.append(f"{typ} {name}")

        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        self.indent()
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in arg_base_names:
                sym = self._visit(symbol)
                symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
        self.dedent()

        rval = f"{self._nindent}{return_type} {node.name}("
        rval = rval + ", ".join(argument_names) + "){\n"
        rval = rval + symbol_list
        rval = rval + body
        rval = rval + self._nindent + "}\n"
        return rval


    def visit_literal_node(self, node: Literal) -> str:
        if node.datatype.intrinsic == ScalarType.Intrinsic.CHARACTER:
            return "\"" + node.value + "\""
        return node.value

    def visit_loop_node(self, node: Loop) -> str:
        start = self._visit(node.start_expr)
        stop = self._visit(node.stop_expr)
        step = self._visit(node.step_expr)
        body = self._visit(node.body)
        string = f"for({node.variable.name} = {start}; {node.variable.name} < {stop}; {node.variable.name} += {step})" + "{\n"
        string = string + body
        string = string + self._nindent + "}"
        return string

    def visit_member_node(self, node: Member) -> str:
        return node.name

    def visit_binaryoperation_node(self, node: BinaryOperation) -> str:
        lhs = self._visit(node.children[0])
        rhs = self._visit(node.children[1])
        operator = ""
        # This should be a map
        if node.operator == BinaryOperation.BinaryOp.ADDITION:
            operator = "+"
        elif node.operator == BinaryOperation.BinaryOp.SUBTRACTION:
            operator = "-"
        elif node.operator == BinaryOperation.BinaryOp.MULTIPLY:
            operator = "*"
        elif node.operator == BinaryOperation.BinaryOp.DIVISION:
            operator = "/"
        elif node.operator == BinaryOperation.BinaryOp.LESS_THAN:
            operator = "<"
        elif node.operator == BinaryOperation.BinaryOp.LESS_THAN_EQUAL:
            operator = "<="
        elif node.operator == BinaryOperation.BinaryOp.GREATER_THAN:
            operator = ">"
        elif node.operator == BinaryOperation.BinaryOp.GREATER_THAN_EQUAL:
            operator = ">="
        elif node.operator == BinaryOperation.BinaryOp.EQUALITY:
            operator = "=="
        elif node.operator == BinaryOperation.BinaryOp.LOG_AND:
            operator = "&&"
        elif node.operator == BinaryOperation.BinaryOp.LOG_OR:
            operator = "||"

        return f"({lhs} {operator} {rhs})"

    def visit_unaryoperation_node(self, node: UnaryOperation) -> str:
        if node.operator == UnaryOperation.UnaryOp.UNARYSUB:
            return f"-{self._visit(node.children[0])}"
        elif node.operator == UnaryOperation.UnaryOp.LOG_NOT:
            return f"(!{self._visit(node.children[0])})"

    def visit_particlepositionreference_node(self, node: ParticlePositionReference) -> str:
        if not self._in_kernel:
            raise NotImplementedError()
        self.addSlice("core_part_position")
        # TODO We should check which particle this accesses when pariwise is
        # supported.
        return f"_core_part_position.access(i, a, {node.dimension})"

    def visit_particlereference_node(self, node: ParticleReference) -> str:
        # TODO We should check which particle this accesses when pariwise is
        # supported.
        slice_name = node.member.name
        self.addSlice(slice_name)

        rstr = f"_{slice_name}.access(i, a"
        if isinstance(node.member, ArrayMixin):
            indices = []
            for index in node.member.indices:
                indices.append(self._visit(index))
            rstr = rstr + ", " + ", ".join(indices) + ")"
        else:
            rstr = rstr + ")"
        return rstr

    def visit_pointerreference_node(self, node: PointerReference) -> str:
        return f"*{node.symbol.name}"

    def visit_scalarreference_node(self, node: ScalarReference) -> str:
        return f"{node.symbol.name}"

    def visit_break_node(self, node: Break) -> str:
        return "break;"

    def visit_return_node(self, node: Return) -> str:
        if len(node.children) == 1:
            return f"return {self._visit(node.children[0])};"
        return "return;"

    def visit_structuremember_node(self, node: StructureMember) -> str:
        return f"{node.name}.{self._visit(node.member)}"

    def visit_structurereference_node(self, node: StructureReference) -> str:
        return f"{node.symbol.name}.{self._visit(node.member)}"

    def visit_while_node(self, node: While) -> str:
        cond = self._visit(node.condition)
        body = self._visit(node.body)
        string = f"while({cond})" + "{\n"
        string = string + body
        string = string + self._nindent + "}"
        return string
