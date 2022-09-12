from __future__ import annotations

from typing import Union

import ast
import re

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType, BOOL_TYPE,\
                                                     type_mapping_str, INT_TYPE

from HartreeParticleDSL.Particle_IR.nodes.symbol_to_reference import symbol_to_reference

from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.array_member import ArrayMember
from HartreeParticleDSL.Particle_IR.nodes.array_mixin import ArrayMixin
from HartreeParticleDSL.Particle_IR.nodes.array_reference import ArrayReference
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.funcdef import FuncDef
from HartreeParticleDSL.Particle_IR.nodes.ifelse import IfElseBlock
from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.loop import Loop
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, \
                                                           UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.nodes.structure_reference import StructureReference
from HartreeParticleDSL.Particle_IR.nodes.structure_member import StructureMember
from HartreeParticleDSL.Particle_IR.nodes.statement import EmptyStatement, Break
from HartreeParticleDSL.Particle_IR.nodes.while_loop import While


from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol_map import datatype_to_symbol

class ast_to_pir_visitor(ast.NodeVisitor):

    def __init__(self):
        super().__init__()
        self._symbol_table = None

    def visit_Add(self, node: ast.Add) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.ADDITION

    def visit_Mult(self, node: ast.Mult) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.MULTIPLY

    def visit_Sub(self, node: ast.Sub) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.SUBTRACTION

    def visit_Div(self, node: ast.Div) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.DIVISION

    def visit_LtE(self, node: ast.LtE) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.LESS_THAN_EQUAL

    def visit_Lt(self, node: ast.Lt) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.LESS_THAN

    def visit_GtE(self, node: ast.GtE) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.GREATER_THAN_EQUAL

    def visit_Gt(self, node: ast.Gt) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.GREATER_THAN

    def visit_USub(self, node: ast.USub) -> UnaryOperation.UnaryOp:
        return UnaryOperation.UnaryOp.UNARYSUB

    def visit_And(self, node: ast.And) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.LOG_AND

    def visit_Or(self, node: ast.Or) -> BinaryOperation.BinaryOp:
        return BinaryOperation.BinaryOp.LOG_OR

    def visit_Not(self, node: ast.Not) -> UnaryOperation.UnaryOp:
        return UnaryOperation.UnaryOp.LOG_NOT

    def visit_BoolOp(self, node: ast.BoolOp) -> BinaryOperation:
        # Check our input is ok, else we assert at the moment.
        assert len(node.values) >= 2
        if len(node.values) == 2:
            operation = self.visit(node.op)
            child1 = self.visit(node.values[0])
            child2 = self.visit(node.values[1])
            return BinaryOperation.create(operation, [child1, child2])
        else:
            last = len(node.values) - 2
            operation = self.visit(node.op)
            child2 = self.visit(node.values[-1])
            child1 = self.visit(node.values[last])
            op = BinaryOperation.create(operation, [child1, child2])
            last = last - 1
            while last >= 0:
                child2 = op
                child1 = self.visit(node.values[last])
                op = BinaryOperation.create(operation, [child1, child2])
                last = last - 1
            return op

    def visit_UnaryOp(self, node: ast.UnaryOp) -> UnaryOperation:
        operation = self.visit(node.op)
        child = self.visit(node.operand)
        return UnaryOperation.create(operation, child)

    def visit_Compare(self, node: ast.Compare) -> BinaryOperation:
        # Check our input is ok, else we assert at the moment.
        assert len(node.ops) == 1
        assert len(node.comparators) == 1
        op = self.visit(node.ops[0])
        child1 = self.visit(node.left)
        child2 = self.visit(node.comparators[0])
        return BinaryOperation.create(op, [child1, child2])

    def visit_BinOp(self, node: ast.BinOp) -> BinaryOperation:
        operation = self.visit(node.op)
        child1 = self.visit(node.left)
        child2 = self.visit(node.right)
        return BinaryOperation.create(operation, [child1, child2])

    def visit_Name(self, node: ast.Name) -> Reference:
        sym = self._symbol_table.lookup(f"{node.id}")
        if sym is None:
            raise IRGenerationError("Attempted to access a symbol that has "
                                    "not been defined in this scope. Symbol "
                                    f"name was {node.id}")
#        if isinstance(sym, ScalarTypeSymbol):
#            return ScalarReference(sym)
        ref = symbol_to_reference[type(sym)](sym)
        return ref

    def visit_Attribute(self, node: ast.Attribute) -> StructureReference:
        # TODO What if Particle symbol is found. Need to create something else.
        attribute_names = []
        attribute_names.append(f"{node.attr}")
        temp_node = node.value
        while isinstance(temp_node, ast.Attribute):
            attribute_names.append(f"{temp_node.attr}")
            temp_node = temp_node.value
        attribute_names.reverse()
        if isinstance(temp_node, ast.Subscript):
            raise IRGenerationError("Array inside a reference with children "
                                    "is not currently supported in ParticleIR.")
        sym_name = f"{temp_node.id}"
        sym = self._symbol_table.lookup(sym_name)
        if len(attribute_names) == 1:
            # Check validity
            if attribute_names[0] not in sym.datatype.components:
                raise IRGenerationError("Attempted to access member "
                                        f"{attribute_names[0]} of structure "
                                        f"type {sym.datatype}.")
            mem = Member(attribute_names[0])
        else:
            inner_member = Member(attribute_names[-1])
            for i in range(len(attribute_names)-2, 0, -1):
                    inner_member = StructureMember(attribute_names[i], inner_member)
            mem = StructureMember(attribute_names[0], inner_member)

            #TODO Check validity.
            dtype = sym.datatype
            current_mem = mem
            while isinstance(current_mem, Member):
                if current_mem.name not in dtype.components:
                    raise IRGenerationError("Attempted to access member "
                                            f"{current_mem.name} of structure "
                                            f"type {dtype}.")
                dtype = dtype.components[current_mem.name]
                if hasattr(current_mem, "member"):
                    current_mem = current_mem.member
                else:
                    current_mem = None
#            raise NotImplementedError()

        return StructureReference(sym, mem)


    def visit_Break(self, node: ast.Break) -> Break:
        return Break()

    def visit_Constant(self, node: ast.Constant) -> Literal:
        if type(node.value) is str:
            dtype = ScalarType(ScalarType.Intrinsic.CHARACTER,
                               ScalarType.Precision.UNDEFINED)
            return Literal(node.value, dtype)
        elif str(node.value) == "True":
            return Literal("True", BOOL_TYPE)
        elif str(node.value) == "False":
            return Literal("False", BOOL_TYPE)
        elif isinstance(node.value, int):
            dtype = ScalarType(ScalarType.Intrinsic.INTEGER,
                               ScalarType.Precision.UNDEFINED)
            return Literal(f"{node.value}", dtype)
        elif isinstance(node.value, float):
            dtype = ScalarType(ScalarType.Intrinsic.FLOAT,
                               ScalarType.Precision.UNDEFINED)
            return Literal(f"{node.value}", dtype)
        else:
            raise NotImplementedError()

    def visit_Assign(self, node: ast.Assign) -> Assignment:
        assert len(node.targets) == 1
        lhs = self.visit(node.targets[0])
        rhs = self.visit(node.value)
        return Assignment.create(lhs, rhs)

    def visit_arg(self, node: ast.arg) -> Reference:
        if node.annotation is not None:
            if f"{node.annotation.id}" not in type_mapping_str:
                raise IRGenerationError("Argument has unexpected type. Got "
                                        f"{node.annotation.id}.")
            types = f"{node.annotation.id}"
            sym = self._symbol_table.new_symbol(f"{node.arg}", type_mapping_str[types],
                                          datatype_to_symbol[type(type_mapping_str[types])])
            ref = symbol_to_reference[type(sym)](sym)

        else:
            print("Got argument without corresponding type. Assuming int")
            sym = self._symbol_table.new_symbol(f"{node.arg}", type_mapping_str["c_int"],
                                          datatype_to_symbol[type(type_mapping_str["c_int"])])
            ref = symbol_to_reference[type(sym)](sym)
        return ref

    def visit_arguments(self, node: ast.arguments) -> List[DataNode]:
        args = []
        for arg in node.args:
            args.append(self.visit(arg))
        return args

    def visit_FunctionDef(self, node: ast.FunctionDef) -> FuncDef:
        arglist = []
        body_nodes = []
        fd = FuncDef.create(f"{node.name}", arglist, body_nodes)
        self._symbol_table = fd.symbol_table
        arglist = self.visit(node.args)
        fd.arguments = arglist
        for child in node.body:
            fd.body.children.append(self.visit(child))
            

        return fd

    def visit_If(self, node: ast.If) -> IfElseBlock:
        elseIf = False
        cond = self.visit(node.test)

        ifbody = []
        for child in node.body:
            ifbody.append(self.visit(child))

        elsebody = []
        for child in node.orelse:
            elsebody.append(self.visit(child))

        ifblock = IfElseBlock.create(cond, ifbody, elsebody)
        # Temp check for now
        return ifblock

    def visit_Call(self, node: ast.Call) -> Union[Call, EmptyStatement, Invoke, Assignment]:
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
        if function_name == "invoke":
            args = []
            for arg in node.args:
                args.append(self.visit(arg))
            return Invoke.create(args)
        elif function_name == "create_variable":
            var_type = f"{node.args[0].id}"
            name = f"{node.args[1].id}"
            if type_mapping_str.get(var_type) is None:
                raise IRGenerationError(
                        f"Provided variable type {var_type} not supported "
                        "in Particle IR.")
            # Only allowing C valid names for now.
            a = re.match("[a-zA-Z_][a-zA-Z_0-9]*", name)
            if a is None or a.group(0) != name:
                raise IRGenerationError(
                        f"Particle IR does not support {name} as a variable "
                        "name.")
            # Add to symbol table
            sym = self._symbol_table.new_symbol(name, type_mapping_str[var_type], datatype_to_symbol[type(type_mapping_str[var_type])])
            if len(node.args) == 2:
                return EmptyStatement()
            else:
                initial_value = node.args[2]
                # Create a Reference to the symbol
                lhs = ScalarReference(sym)
                rhs = self.visit(node.args[2])
                return Assignment.create(lhs, rhs)
        else:
            args = []
            for arg in node.args:
                args.append(self.visit(arg))
            return Call.create(function_name, args)

    def visit_Module(self, node: ast.Module) -> Node:
        module_nodes = []
        for a in ast.iter_child_nodes(node):
            module_nodes.append(self.visit(a))
        return module_nodes[0]

    # 3.8 only
    def visit_Index(self, node: ast.Index) -> DataNode:
        return self.visit(node.value)

    def visit_Subscript(self, node: ast.Subscript) -> Union[List[Node],Node]:
        if isinstance(node.value, ast.Attribute):
            temp_ref = self.visit(node.value)
            x = temp_ref.member
            #while hasattr(x, "member"):
#                x = temp_ref.member
#                temp_ref = temp_ref.member
            assert not hasattr(x, "member")
            # x is always the innermost member, use assert to double check but
            # couldn't create cases otherwise
            indices = self.visit(node.slice)
            if not isinstance(indices, list):
                indices = [indices]
            new_x = ArrayMember(x.name, indices)
            temp_ref.member = new_x
            return temp_ref
        elif isinstance(node.value, ast.Subscript):
            index = self.visit(node.slice)
            x = self.visit(node.value)
            if isinstance(x, StructureReference):
                child = x
                # Most cases here are already handled in visit_Attribute
                # we assert to ensure but
                # child.member should always be ArrayMixin
                assert isinstance(child.member, ArrayMixin)
#                while not isinstance(child.member, ArrayMixin):
#                    child = child.member
                child.member.indices.append(index)
            elif isinstance(x, ArrayReference):
                x.indices.append(index)
#            else:
#                raise NotImplementedError()
            return x
        else:
            r = self.visit(node.value)
            index = self.visit(node.slice)
            r.indices.append(index)
            return r

    def visit_For(self, node: ast.For) -> Loop:
        if len(node.orelse) != 0:
            raise IRGenerationError(f"Else clauses on Loops are not supported in ParticleIR")
        if type(node.iter) != ast.Call:
            raise IRGenerationError(f"Only range loops are supported in ParticleIR")
        if type(node.iter) == ast.Call:
            if node.iter.func.id != "range":
                raise IRGenerationError(f"Only range loops are supported in ParticleIR")

        if len(node.iter.args) == 2:
            startval = self.visit(node.iter.args[0])
            stopval = self.visit(node.iter.args[1])
            stepval = Literal("1", INT_TYPE)
        elif len(node.iter.args) == 3:
            startval = self.visit(node.iter.args[0])
            stopval = self.visit(node.iter.args[1])
            stepval = self.visit(node.iter.args[2])
        elif len(node.iter.args) == 1:
            startval = Literal("0", INT_TYPE)
            stopval = self.visit(node.iter.args[0])
            stepval = Literal("1", INT_TYPE)

        sym = self._symbol_table.find_or_create(node.target.id, INT_TYPE, ScalarTypeSymbol)

        body_nodes = []
        for child in node.body:
            body_nodes.append(self.visit(child))

        return Loop.create(sym, startval, stopval, stepval, body_nodes)

    def visit_While(self, node: ast.While) -> While:
        test = self.visit(node.test)
        body = []
        for child in node.body:
            body.append(self.visit(child))
        return While.create(test, body)

    def visit_Expr(self, node: ast.Expr) -> Node:
        exprs = []
        for a in ast.iter_child_nodes(node):
            exprs.append(self.visit(a))
        return exprs[0]

    def generic_visit(self, node: Any):
        raise IRGenerationError(f"Found unsupported node of type {type(node)}")

