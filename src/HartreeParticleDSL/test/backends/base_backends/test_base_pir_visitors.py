import pytest
from HartreeParticleDSL.Particle_IR.nodes.statement import EmptyStatement,\
        Break
from HartreeParticleDSL.backends.base_backend.pir_visitor import *

def test_base_pir_init():
    a = PIR_Visitor()
    assert a._nindent == ""

    a = PIR_Visitor(initial_indent_depth=2)
    assert a._nindent == "        "

    a = PIR_Visitor(indent_string="  ", initial_indent_depth=1)
    assert a._nindent == "  "

    with pytest.raises(TypeError) as excinfo:
        a = PIR_Visitor(indent_string=2)
    assert ("Expected indent_string to be a string but got "
            "<class 'int'>." in str(excinfo.value))
    with pytest.raises(TypeError) as excinfo:
        a = PIR_Visitor(initial_indent_depth="  ")
    assert ("Expected initial_indent_depth to be an int but got "
            "<class 'str'>." in str(excinfo.value))

def test_base_pir_emptystatement():
    a = PIR_Visitor()
    es = EmptyStatement()
    assert a._visit(es) == ""

def test_base_pir_errors():
    a = PIR_Visitor()
    with pytest.raises(TypeError) as excinfo:
        a(2)
    assert ("Expected a Particle IR Node as input but got <class 'int'>."
            in str(excinfo.value))

    with pytest.raises(VisitorError) as excinfo:
        a(Break())
    print(str(excinfo.value))
    assert ("Unsupported node 'Break' found: method names attempted were "
            "['visit_break_node', 'visit_statement_node', 'visit_node_node']"
            in str(excinfo.value))

    class temp_PIR_visitor(PIR_Visitor):
        def visit_break_node(self, node):
            raise AttributeError("Child Attribute Error")

    b = temp_PIR_visitor()
    with pytest.raises(AttributeError) as excinfo:
        b(Break())
    assert ("Child Attribute Error" in str(excinfo.value))
