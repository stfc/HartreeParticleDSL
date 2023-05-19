from __future__ import annotations

import pytest

from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import *
from HartreeParticleDSL.backends.base_backend.backend import Backend

import inspect
import textwrap
import ast

from HartreeParticleDSL.HartreeParticleDSL import set_backend

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType, BOOL_TYPE,\
                                                     type_mapping_str, StructureType,\
                                                     INT_TYPE, ArrayType, FLOAT_TYPE

from HartreeParticleDSL.Particle_IR.nodes.array_reference import ArrayReference
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.config_reference import ConfigReference
from HartreeParticleDSL.Particle_IR.nodes.funcdef import FuncDef
from HartreeParticleDSL.Particle_IR.nodes.ifelse import IfElseBlock
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.nodes.loop import Loop
from HartreeParticleDSL.Particle_IR.nodes.kernels import MainKernel, PairwiseKernel, \
                                                         PerPartKernel
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, \
                                                           UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.statement import EmptyStatement, Break
from HartreeParticleDSL.Particle_IR.nodes.structure_reference import StructureReference
from HartreeParticleDSL.Particle_IR.nodes.while_loop import While

from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.c_types import c_double, c_int

from psyclone.psyir.symbols import Symbol


def test_find_calls_visitor():
    def a():
        my_call()
        my_call.the_call()
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    v = find_calls_visitor()
    v.visit(c)

    calls = v.calls
    assert "my_call" in calls
    assert "my_call.the_call" in calls

def test_visit_str():
    pass


def test_visit_Name():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1 + 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert assign.lhs.symbol.name == "b"

    def b():
        b = 3
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    with pytest.raises(IRGenerationError) as excinfo:
        v.visit(c)
    assert ("Attempted to access a symbol that has not been defined in this "
            "scope. Symbol name was b" in str(excinfo.value))
    
    def d():
        mycall(hello)
    # Occasionally used in calls and similar
    v._check_valid = False
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    out = v.visit(c)
    assert isinstance(out.body.children[0], Call)
    assert isinstance(out.body.children[0].children[0], ScalarReference)
    assert out.body.children[0].children[0].symbol.datatype == STRING_TYPE
    assert out.body.children[0].children[0].symbol.name == "hello"

def test_visit_Add():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1 + 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.ADDITION

def test_visit_Mult():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1 * 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.MULTIPLY

def test_visit_Sub():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1 - 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.SUBTRACTION


def test_visit_Div():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1 / 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.DIVISION

def test_visit_LtE():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = 1 <= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LESS_THAN_EQUAL

def test_visit_Lt():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = 1 < 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LESS_THAN

def test_visit_GtE():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = 1 >= 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.GREATER_THAN_EQUAL

def test_visit_Eq():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = 1 == 2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.EQUALITY


def test_visit_Gt_and_Compare():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = (1 > 2)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.GREATER_THAN

def test_visit_Usub_and_UnaryOp():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        b = -1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.rhs, UnaryOperation)
    assert assign.rhs.operator == UnaryOperation.UnaryOp.UNARYSUB

def test_visit_And():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        create_variable(c_bool, z)
        create_variable(c_bool, y)
        b = (z and y)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], EmptyStatement)
    assert isinstance(out.body.children[3], Assignment)
    assign = out.body.children[3]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LOG_AND

def test_visit_Or():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        create_variable(c_bool, z)
        create_variable(c_bool, y)
        b = (z or y)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], EmptyStatement)
    assert isinstance(out.body.children[3], Assignment)
    assign = out.body.children[3]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LOG_OR

def test_visit_Or():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        create_variable(c_bool, z)
        b = (not z)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], Assignment)
    assign = out.body.children[2]
    assert isinstance(assign.rhs, UnaryOperation)
    assert assign.rhs.operator == UnaryOperation.UnaryOp.LOG_NOT

def test_visit_BoolOp():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_bool, b)
        create_variable(c_bool, z)
        create_variable(c_bool, y)
        b = (z and y and y)
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], EmptyStatement)
    assert isinstance(out.body.children[3], Assignment)
    assign = out.body.children[3]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LOG_AND
    assert isinstance(assign.rhs.children[0], ScalarReference)
    assert assign.rhs.children[0].symbol.name == "z"
    assert isinstance(assign.rhs.children[1], BinaryOperation)
    inner_op = assign.rhs.children[1]
    assert inner_op.operator == BinaryOperation.BinaryOp.LOG_AND
    assert isinstance(inner_op.children[0], ScalarReference)
    assert inner_op.children[0].symbol.name == "y"
    assert isinstance(inner_op.children[1], ScalarReference)
    assert inner_op.children[1].symbol.name == "y"

    def b():
        create_variable(c_bool, b)
        create_variable(c_bool, z)
        create_variable(c_bool, y)
        b = (z and y or y)
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], EmptyStatement)
    assert isinstance(out.body.children[3], Assignment)
    assign = out.body.children[3]
    assert isinstance(assign.rhs, BinaryOperation)
    assert assign.rhs.operator == BinaryOperation.BinaryOp.LOG_OR
    assert isinstance(assign.rhs.children[1], ScalarReference)
    assert assign.rhs.children[1].symbol.name == "y"
    assert isinstance(assign.rhs.children[0], BinaryOperation)
    inner_op = assign.rhs.children[0]
    assert inner_op.operator == BinaryOperation.BinaryOp.LOG_AND
    assert isinstance(inner_op.children[0], ScalarReference)
    assert inner_op.children[0].symbol.name == "z"
    assert isinstance(inner_op.children[1], ScalarReference)
    assert inner_op.children[1].symbol.name == "y"

def test_visit_Attribute():
    v = ast_to_pir_visitor()
    # Create a structure type
    mystruc1 = StructureType()
    mystruc1.add("d", INT_TYPE, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc1"]=mystruc1
    def a():
        create_variable(mystruc1, b)
        b.d = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, StructureReference)
    assert isinstance(assign.lhs.member, Member)
    assert assign.lhs.symbol.name == "b"
    assert assign.lhs.member.name == "d"

    def aa():
        create_variable(mystruc1, b)
        b.e = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(aa)))

    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Attempted to access member e of structure type StructureType<>."
            in str(excinfo.value))
    del(type_mapping_str["mystruc1"])
    assert type_mapping_str.get("mystruc1") is None

    mystruc2 = StructureType()
    substruc = StructureType()
    substruc.add("e", INT_TYPE, Symbol.Visibility.PUBLIC)
    mystruc2.add("d", substruc, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2

    def ab():
        create_variable(mystruc2, b)
        b.d.e = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(ab)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, StructureReference)
    assert isinstance(assign.lhs.member, StructureMember)
    assert isinstance(assign.lhs.member.member, Member)
    assert assign.lhs.symbol.name == "b"
    assert assign.lhs.member.name == "d"
    assert assign.lhs.member.member.name == "e"

    del(type_mapping_str["mystruc2"])
    assert type_mapping_str.get("mystruc2") is None

    mystruc2 = StructureType()
    substruc1 = StructureType()
    substruc2 = StructureType()
    substruc2.add("e", INT_TYPE, Symbol.Visibility.PUBLIC)
    substruc1.add("d", substruc2, Symbol.Visibility.PUBLIC)
    mystruc2.add("c", substruc1, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2

    def ac():
        create_variable(mystruc2, b)
        b.c.d.f = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(ac)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Attempted to access member f of structure type StructureType<"
            ">." in str(excinfo.value))

    del(type_mapping_str["mystruc2"])
    assert type_mapping_str.get("mystruc2") is None

    mystruc2 = StructureType()
    substruc1 = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC])
    mystruc2.add("c", substruc1, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2

    def ad():
        create_variable(mystruc2, b)
        b.c[2] = 3

    c = ast.parse(textwrap.dedent(inspect.getsource(ad)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, StructureReference)
    assert isinstance(assign.lhs.member, ArrayMember)
    assert assign.lhs.symbol.name == "b"
    assert assign.lhs.member.name == "c"
    assert isinstance(assign.lhs.member.indices[0], Literal)
    assert assign.lhs.member.indices[0].value == "2"

    del(type_mapping_str["mystruc2"])
    assert type_mapping_str.get("mystruc2") is None

    mystruc2 = StructureType()
    substruc1 = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC, ArrayType.Extent.DYNAMIC])
    mystruc2.add("c", substruc1, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2

    def ae():
        create_variable(mystruc2, b)
        b.c[3][4] = 2

    c = ast.parse(textwrap.dedent(inspect.getsource(ae)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, StructureReference)
    assert isinstance(assign.lhs.member, ArrayMember)
    assert assign.lhs.symbol.name == "b"
    assert assign.lhs.member.name == "c"
    assert isinstance(assign.lhs.member.indices[0], Literal)
    assert assign.lhs.member.indices[0].value == "3"
    assert isinstance(assign.lhs.member.indices[1], Literal)
    assert assign.lhs.member.indices[1].value == "4"

    def ag():
        create_variable(mystruc2,b)
        b[1].c[2] = 2
    c = ast.parse(textwrap.dedent(inspect.getsource(ag)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Array inside a reference with children is not currently "
            "supported in ParticleIR." in str(excinfo.value))

    def ah():
        b.c = 3
    c = ast.parse(textwrap.dedent(inspect.getsource(ah)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Attempted to access a symbol that has not been defined in this "
            "scope. Symbol name was b" in str(excinfo.value))

    part = type_mapping_str["part"]
    config = type_mapping_str["config"]
    def aj(a: part, conf: config):
        a.core_part.position.x = 1.0
        a.core_part.position.y = 2.0
        a.core_part.position.z = 3.0
        conf.space.box_dims.x_min = 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(aj)))
    out = v.visit(c)
    assign1 = out.body.children[0]
    assert isinstance(assign1.lhs, ParticlePositionReference)
    assert assign1.lhs.dimension == 0
    assign2 = out.body.children[1]
    assert isinstance(assign2.lhs, ParticlePositionReference)
    assert assign2.lhs.dimension == 1
    assign3 = out.body.children[2]
    assert isinstance(assign3.lhs, ParticlePositionReference)
    assert assign3.lhs.dimension == 2
    assign4 = out.body.children[3]
    assert isinstance(assign4.lhs, ConfigReference)

    def ak(a: part):
        a.core_part.position.p = 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(ak)))
    with pytest.raises(IRGenerationError) as excinfo:
       out = v.visit(c)
    assert ("Attempted to access unknown element of particle position. "
            "Got p" in str(excinfo.value))



def test_visit_Assign():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, b)
        b = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, ScalarReference)
    assert assign.lhs.symbol.name == "b"

def test_visit_arg_and_arguments_and_FunctionDef(capsys):
    v = ast_to_pir_visitor()
    def a(arg: c_int, arg2: c_float):
        create_variable(c_int, b)
        b = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 2
    assert isinstance(out.arguments[0], ScalarReference)
    assert out.arguments[0].symbol.name == "arg"
    assert out.arguments[0].symbol.datatype == INT_TYPE
    assert out.arguments[1].symbol.name == "arg2"
    assert out.arguments[1].symbol.datatype == FLOAT_TYPE


    mystruc2 = StructureType()
    substruc1 = StructureType()
    substruc2 = StructureType()
    substruc2.add("e", INT_TYPE, Symbol.Visibility.PUBLIC)
    substruc1.add("d", substruc2, Symbol.Visibility.PUBLIC)
    mystruc2.add("c", substruc1, Symbol.Visibility.PUBLIC)
    type_mapping_str["mystruc2"]=mystruc2

    def b(arg: mystruc2):
        create_variable(c_int, c)

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 1
    assert isinstance(out.arguments[0], StructureReference)
    assert out.arguments[0].symbol.name == "arg"
    assert out.arguments[0].symbol.datatype == mystruc2

    del(type_mapping_str["mystruc2"])
    assert type_mapping_str.get("mystruc2") is None

    at = ArrayType(INT_TYPE, [1])
    type_mapping_str["at"] = at
    def d(arg: at):
        create_variable(c_int, c)

    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 1
    assert isinstance(out.arguments[0], ArrayReference)
    assert out.arguments[0].symbol.name == "arg"
    assert out.arguments[0].symbol.datatype == at
    del(type_mapping_str["at"])

    def e(arg):
        create_variable(c_int, c)

    c = ast.parse(textwrap.dedent(inspect.getsource(e)))
    out = v.visit(c)
    captured = capsys.readouterr()
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 1
    assert isinstance(out.arguments[0], ScalarReference)
    assert out.arguments[0].symbol.name == "arg"
    assert out.arguments[0].symbol.datatype == INT_TYPE
    assert "Got argument without corresponding type. Assuming int" in captured.out

    def f(arg: fake):
        create_variable(c_int, c)
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Argument has unexpected type. Got fake." in str(excinfo.value))

    def main():
        create_variable(c_int, c)
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    out = v.visit(c)
    symbols = out.symbol_table.get_symbols()
    assert "config" in symbols

def test_visit_If():
    v = ast_to_pir_visitor()
    def a():
        create_variable(c_int, c)
        if c < 1:
            c = c + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], IfElseBlock)

    assert isinstance(out.body.children[1].ifbody, Body)
    assert isinstance(out.body.children[1].ifbody.children[0], Assignment)
    assert len(out.body.children[1].elsebody.children) == 0

    def b():
        create_variable(c_int, c)
        if c < 1:
            c = c + 1
        else:
            c = c + 2

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], IfElseBlock)

    assert isinstance(out.body.children[1].ifbody, Body)
    assert isinstance(out.body.children[1].ifbody.children[0], Assignment)
    assert out.body.children[1].ifbody.children[0].rhs.children[1].value == "1"
    assert isinstance(out.body.children[1].elsebody.children[0], Assignment)
    assert out.body.children[1].elsebody.children[0].rhs.children[1].value == "2"

    def d():
        create_variable(c_int, c)
        if c < 1:
            c = c + 1
        elif c < 2:
            c = c + 2

    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], IfElseBlock)

    assert isinstance(out.body.children[1].ifbody, Body)
    assert out.body.children[1].is_else_if()
    assert isinstance(out.body.children[1].ifbody.children[0], Assignment)
    assert out.body.children[1].ifbody.children[0].rhs.children[1].value == "1"
    assert isinstance(out.body.children[1].elsebody.children[0], IfElseBlock)
    inner_if = out.body.children[1].elsebody.children[0]
    assert isinstance(inner_if.ifbody.children[0], Assignment)
    assert inner_if.ifbody.children[0].rhs.children[1].value == "2"
    assert len(inner_if.elsebody.children) == 0

    def e():
        create_variable(c_int, c)
        if c < 1:
            c = c + 1
        elif c < 2:
            c = c + 2
        else:
            c = c + 3

    c = ast.parse(textwrap.dedent(inspect.getsource(e)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], IfElseBlock)

    assert isinstance(out.body.children[1].ifbody, Body)
    assert out.body.children[1].is_else_if()
    assert isinstance(out.body.children[1].ifbody.children[0], Assignment)
    assert out.body.children[1].ifbody.children[0].rhs.children[1].value == "1"
    assert isinstance(out.body.children[1].elsebody.children[0], IfElseBlock)
    inner_if = out.body.children[1].elsebody.children[0]
    assert isinstance(inner_if.ifbody.children[0], Assignment)
    assert inner_if.ifbody.children[0].rhs.children[1].value == "2"
    assert isinstance(inner_if.elsebody.children[0], Assignment)
    assert inner_if.elsebody.children[0].rhs.children[1].value == "3"


def test_visit_Call():

    v = ast_to_pir_visitor()
    def a():
        generic_call()

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], Call)
    assert out.body.children[0].func_name == "generic_call"
    assert len(out.body.children[0].children) == 0

    v = ast_to_pir_visitor()
    def b(x: c_int, y: c_int):
        generic_call_2(x, y)

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], Call)
    assert out.body.children[0].func_name == "generic_call_2"
    assert len(out.body.children[0].children) == 2
    assert out.body.children[0].children[0].symbol.name == "x"
    assert out.body.children[0].children[1].symbol.name == "y"

    def f():
        generic_call.call()

    c = ast.parse(textwrap.dedent(inspect.getsource(f)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], Call)
    assert out.body.children[0].func_name == "generic_call.call"
    assert len(out.body.children[0].children) == 0

    def g():
        invoke(mykern)


    c = ast.parse(textwrap.dedent(inspect.getsource(g)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], Invoke)
    assert len(out.body.children[0].children) == 1

    def h():
        create_variable(fake, a)

    c = ast.parse(textwrap.dedent(inspect.getsource(h)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Provided variable type fake not supported in Particle IR."
            in str(excinfo.value))

    # There isn't a known way to create an invalid variable name
    # inside python, but we're checking just in case.
    class fake_arg():
        def __init__(self, the_id):
            self.id = the_id
    class fake_Call():
        def __init__(self):
            import ast
            self.func = ast.Name("create_variable")
            self.args = []
            self.args.append(fake_arg("c_int"))
            self.args.append(fake_arg("90xsdawer"))

    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit_Call(fake_Call())
    assert ("Particle IR does not support 90xsdawer as a variable name."
            in str(excinfo.value))

    def k():
        create_variable(c_int, a, 3)
    c = ast.parse(textwrap.dedent(inspect.getsource(k)))
    out = v.visit(c)

    assert isinstance(out.body.children[0], Assignment)
    assert out.body.children[0].lhs.symbol.name == "a"
    assert out.body.children[0].lhs.symbol.datatype == INT_TYPE
    assert out.body.children[0].rhs.value == "3"



def test_visit_Module():
    v = ast_to_pir_visitor()
    def a():
        generic_call()

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)

    assert isinstance(out, FuncDef)
    assert out.name == "a"

def test_visit_Index_and_Subscript():
    v = ast_to_pir_visitor()
    substruc1 = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC, ArrayType.Extent.DYNAMIC])
    type_mapping_str["array"] = substruc1

    def a():
        create_variable(array, c)
        c[0][2] = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, ArrayReference)
    assert assign.lhs.symbol.name == "c"
    assert isinstance(assign.lhs.indices[0], Literal)
    assert assign.lhs.indices[0].value == "0"
    assert isinstance(assign.lhs.indices[1], Literal)
    assert assign.lhs.indices[1].value == "2"

    struct = StructureType()
    substruct = StructureType()
    substruct.add("array", ArrayType(INT_TYPE, [3]), Symbol.Visibility.PUBLIC)
    struct.add("sub", substruct, Symbol.Visibility.PUBLIC)
    type_mapping_str["struc"] = struct

    def b():
        create_variable(struc, c)
        c.sub.array[0] = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assign = out.body.children[1]
    assert isinstance(assign.lhs, StructureReference)
    assert isinstance(assign.lhs.member, StructureMember)
    assert isinstance(assign.lhs.member.member, ArrayMixin)
    assert assign.lhs.member.member.indices[0].value == "0"

    # Cover pre-3.9 code to ensure 3.8 installations are tested.
    x = ast.Constant(1)
    class fake_Node():
        def __init__(self, x):
            self.value = x
    y = fake_Node(x)
    out = v.visit_Index(y)
    assert isinstance(out, Literal)
    assert out.value == "1"

    del(type_mapping_str["array"])

def test_visit_For():
    v = ast_to_pir_visitor()
    substruc1 = ArrayType(FLOAT_TYPE, [ArrayType.Extent.DYNAMIC])
    type_mapping_str["array"] = substruc1
    def a():
        create_variable(array, d)
        create_variable(c_int, i)
        create_variable(c_int, j)
        create_variable(c_int, k)
        create_variable(c_int, s)
        for i in range(3):
            d[i] = 0.2
        for j in range(3,5):
            d[j] = 0.2
        for k in range(3,0,-1):
            d[k] = 0.2
        for s in range(3, 9, 2):
            d[s] = 0.2
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], EmptyStatement)
    assert isinstance(out.body.children[2], EmptyStatement)
    assert isinstance(out.body.children[3], EmptyStatement)
    assert isinstance(out.body.children[4], EmptyStatement)
    loop1 = out.body.children[5]
    assert isinstance(loop1, Loop)
    assert loop1.variable.name == "i"
    assert isinstance(loop1.start_expr, Literal)
    assert loop1.start_expr.value == "0"
    assert isinstance(loop1.stop_expr, Literal)
    assert loop1.stop_expr.value == "3"
    assert isinstance(loop1.step_expr, Literal)
    assert loop1.step_expr.value == "1"
    assert len(loop1.body.children) == 1

    loop2 = out.body.children[6]
    assert isinstance(loop2, Loop)
    assert loop2.variable.name == "j"
    assert isinstance(loop2.start_expr, Literal)
    assert loop2.start_expr.value == "3"
    assert isinstance(loop2.stop_expr, Literal)
    assert loop2.stop_expr.value == "5"
    assert isinstance(loop2.step_expr, Literal)
    assert loop2.step_expr.value == "1"
    assert len(loop2.body.children) == 1

    loop3 = out.body.children[7]
    assert isinstance(loop3, Loop)
    assert loop3.variable.name == "k"
    assert isinstance(loop3.start_expr, Literal)
    assert loop3.start_expr.value == "3"
    assert isinstance(loop3.stop_expr, Literal)
    assert loop3.stop_expr.value == "0"
    assert isinstance(loop3.step_expr, UnaryOperation)
    assert loop3.step_expr.operator == UnaryOperation.UnaryOp.UNARYSUB 
    assert loop3.step_expr.children[0].value == "1"
    assert len(loop3.body.children) == 1

    loop4 = out.body.children[8]
    assert isinstance(loop4, Loop)
    assert loop4.variable.name == "s"
    assert isinstance(loop4.start_expr, Literal)
    assert loop4.start_expr.value == "3"
    assert isinstance(loop4.stop_expr, Literal)
    assert loop4.stop_expr.value == "9"
    assert isinstance(loop4.step_expr, Literal)
    assert loop4.step_expr.value == "2"
    assert len(loop4.body.children) == 1


    def b():
        create_variable(c_int, i)
        create_variable(array, d)
        create_variable(c_int, x)
        for i in range(3):
            d[i] = 0.2
        else:
            x = 0

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Else clauses on Loops are not supported in ParticleIR" in
            str(excinfo.value))

    def d():
        create_variable(c_int, i)
        create_variable(array, a)
        create_variable(array, d)
        for i in a:
            d[i] = 0.1

    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Only range loops are supported in ParticleIR" in
            str(excinfo.value))

    def e():
        create_variable(c_int, i)
        create_variable(array, a)
        create_variable(array, d)
        for i in a():
            d[i] = 0.1
    c = ast.parse(textwrap.dedent(inspect.getsource(e)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Only range loops are supported in ParticleIR" in
            str(excinfo.value))
    del(type_mapping_str["array"])

def test_visit_While():

    v = ast_to_pir_visitor()
    def a(b: c_bool):
        while b:
            b = not b

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 1
    assert isinstance(out.body.children[0], While)
    whileloop = out.body.children[0]
    assert isinstance(whileloop.condition, ScalarReference)
    assert len(whileloop.body.children) == 1

def test_visit_Break():
    v = ast_to_pir_visitor()
    def a(b: c_bool):
        while b:
            break

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 1
    assert isinstance(out.body.children[0], While)
    whileloop = out.body.children[0]
    assert isinstance(whileloop.condition, ScalarReference)
    assert len(whileloop.body.children) == 1
    assert isinstance(whileloop.body.children[0], Break)

def test_visit_Constant():
    v = ast_to_pir_visitor()
    def a(b:c_bool, c:c_bool):
        b = True
        c = False

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out, FuncDef)
    assert len(out.arguments) == 2
    assert isinstance(out.body.children[0].rhs, Literal)
    assert out.body.children[0].rhs.value == "True"
    assert isinstance(out.body.children[1].rhs, Literal)
    assert out.body.children[1].rhs.value == "False"

    class temp():
        def __init__(self):
            self.value = None
    with pytest.raises(NotImplementedError):
        v.visit_Constant(temp())

def test_visit_Return():
    v = ast_to_pir_visitor()
    def a():
        return 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out.body.children[0], Return)
    assert len(out.body.children[0].children) == 1
    assert out.body.children[0].children[0].value == "1"

    def b():
        return
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)
    assert isinstance(out.body.children[0], Return)
    assert len(out.body.children[0].children) == 0

def test_particle_visitors():
    v = ast_to_pir_visitor()
    def a(part1: part):
        part1.core_part.position[0] = 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    out = v.visit(c)
    assert isinstance(out.body.children[0].lhs, ParticlePositionReference)
    assert out.body.children[0].lhs.symbol.name == "part1"
    assert out.body.children[0].lhs.dimension == 0

    type_mapping_str["part"].components["test"] = INT_TYPE

    def b(part1: part):
        part1.test = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    out = v.visit(c)

    assert isinstance(out.body.children[0].lhs, ParticleReference)
    assert out.body.children[0].lhs.symbol.name == "part1"
    assert out.body.children[0].lhs.member.name == "test"

    type_mapping_str["part"].components["x"] = ArrayType(INT_TYPE, [3])
    def d(part1: part):
        part1.x[0] = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    out = v.visit(c)
    assert isinstance(out.body.children[0].lhs, ParticleReference)
    assert out.body.children[0].lhs.symbol.name == "part1"
    assert out.body.children[0].lhs.member.name == "x"
    assert out.body.children[0].lhs.member.indices[0].value == "0"

    type_mapping_str["part"].components["x"] = ArrayType(INT_TYPE, [3])
    def e(part1: part, part2: part):
        part2 = part1
    c = ast.parse(textwrap.dedent(inspect.getsource(e)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert("Particle IR doesn't currently support accessing a full particle "
           "type in a single statement." in str(excinfo.value))

    def e(part1: part):
        part1.core_part.a.b = 1
    c = ast.parse(textwrap.dedent(inspect.getsource(e)))
    out = v.visit(c)
    assert isinstance(out.body.children[0].lhs, ParticleCoreReference)
    assert out.body.children[0].lhs.symbol.name == "part1"
    assert isinstance(out.body.children[0].lhs.member, StructureMember)




def test_visit_generic():
    v = ast_to_pir_visitor()

    with pytest.raises(IRGenerationError) as excinfo:
        v.visit(v)
    assert ("Found unsupported node of type <class 'HartreeParticleDSL.backends."
            "AST_to_Particle_IR.ast_to_pir_visitors.ast_to_pir_visitor'>"
            in str(excinfo.value))

    def e(part1: part):
        create_variable(part, part2, part1)



def test_visit_generic():
    v = ast_to_pir_visitor()

    with pytest.raises(IRGenerationError) as excinfo:
        v.visit(v)
    assert ("Found unsupported node of type <class 'HartreeParticleDSL.backends."
            "AST_to_Particle_IR.ast_to_pir_visitors.ast_to_pir_visitor'>"
            in str(excinfo.value))

def test_main_visitor():
    v = pir_main_visitor()

    def main(arg: c_int):
        create_variable(c_int, a)

    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Particle IR expects no arguments to main function but "
            "received 1." in str(excinfo.value))

    class temp_Backend(Backend):
        def get_extra_symbols(self, arg):
            return [ ("b", ScalarTypeSymbol("b", INT_TYPE))]
    set_backend(temp_Backend())
    
    def a():
        create_variable(c_int, a)

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Attempting to create a main function but name was not "
            "'main', got 'a'." in str(excinfo.value))

    def main():
        create_variable(c_int, a)
        b = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    out = v.visit(c)

    assert isinstance(out, MainKernel)
    assert len(out.body.children) == 2
    assert isinstance(out.body.children[0], EmptyStatement)
    assert isinstance(out.body.children[1], Assignment)
    assert out.body.children[1].lhs.symbol.name == "b"

def test_pairwise_visitor():
    v = pir_pairwise_visitor()

    class temp_Backend(Backend):
        def get_extra_symbols(self, arg):
            return [ ("b", ScalarTypeSymbol("b", INT_TYPE))]
        def get_writable_arrays(self):
            return {"temp": ("c_double", 10)}
    set_backend(temp_Backend())

    def x(arg: part, arg2: part, c: c_int):
        create_variable(c_int, a)
        b = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    out = v.visit(c)

    assert isinstance(out.arguments[0], ParticleReference)
    assert isinstance(out.arguments[1], ParticleReference)
    assert len(out.body.children) == 2
    assert out.body.children[1].lhs.symbol.name == "b"

    def y(arg: part):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(y)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Particle IR expects at least 3 arguments for a pairwise kernel, "
            "but got 1." in str(excinfo.value))

    def z(arg: part, arg2: c_int, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Pairwise Kernel needs 2 particle arguments in Particle IR."
            in str(excinfo.value))
   
def test_perpart_visitor():
    v = pir_perpart_visitor()

    class temp_Backend(Backend):
        def get_extra_symbols(self, arg):
            return [ ("b", ScalarTypeSymbol("b", INT_TYPE))]
        def get_writable_arrays(self):
            return {"temp": ("c_double", 10)}
    set_backend(temp_Backend())

    def x(arg: part, c: c_int):
        create_variable(c_int, a)
        b = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    out = v.visit(c)

    assert isinstance(out.arguments[0], ParticleReference)
    assert len(out.body.children) == 2
    assert out.body.children[1].lhs.symbol.name == "b"

    def y(arg: part):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(y)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Particle IR expects at least 2 arguments for a perpart kernel, "
            "but got 1." in str(excinfo.value))

    def z(arg: c_int, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("First argument to PerPartKernel must be a particle."
            in str(excinfo.value))

def test_source_boundary_visitor():
    class temp_class():
        def get_source_count(self):
            return 1
    v = pir_source_boundary_visitor(temp_class())

    class temp_Backend(Backend):
        def get_extra_symbols(self, arg):
            return [ ("b", ScalarTypeSymbol("b", INT_TYPE))]
        def get_writable_arrays(self):
            return {"temp": ("c_double", 10)}
    set_backend(temp_Backend())

    def y(arg: part):
        create_variable(c_int, a)

    c = ast.parse(textwrap.dedent(inspect.getsource(y)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Particle IR expects at least 2 arguments for a source boundary kernel, "
            "but got 1." in str(excinfo.value))


    def z(arg: c_int, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("First argument to SourceBoundaryKernel must be a particle."
            in str(excinfo.value))

    def z2(arg: part, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z2)))
    out = v.visit(c)
    assert isinstance(out, SourceBoundaryKernel)
    assert len(out.children[0].children) == 1
    assert isinstance(out.children[0].children[0], EmptyStatement)

def test_sink_boundary_visitor():
    v = pir_sink_boundary_visitor()

    class temp_Backend(Backend):
        def get_extra_symbols(self, arg):
            return [ ("b", ScalarTypeSymbol("b", INT_TYPE))]
        def get_writable_arrays(self):
            return {"temp": ("c_double", 10)}
    set_backend(temp_Backend())

    def y(arg: part):
        create_variable(c_int, a)

    c = ast.parse(textwrap.dedent(inspect.getsource(y)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("Particle IR expects at least 2 arguments for a sink boundary kernel, "
            "but got 1." in str(excinfo.value))


    def z(arg: c_int, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    with pytest.raises(IRGenerationError) as excinfo:
        out = v.visit(c)
    assert ("First argument to SinkBoundaryKernel must be a particle."
            in str(excinfo.value))

    def z2(arg: part, c: c_int):
        create_variable(c_int, a)
    c = ast.parse(textwrap.dedent(inspect.getsource(z2)))
    out = v.visit(c)
    assert isinstance(out, SinkBoundaryKernel)
    assert len(out.children[0].children) == 1
    assert isinstance(out.children[0].children[0], EmptyStatement)
