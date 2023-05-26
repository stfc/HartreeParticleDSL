import pytest

from psyclone.psyir.nodes import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import FLOAT_TYPE, INT_TYPE, BOOL_TYPE

def test_literal_init_bad_inputs():

    with pytest.raises(TypeError) as excinfo:
        Literal(123, INT_TYPE)
    assert ("Literals must be supplied with a value encoded as a string "
            "but found 'int'" in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        Literal("123", 123)
    assert ("The datatype of a Literal must be an instance of "
            "psyir.symbols.ScalarType or psyir.symbols.ArrayType but found 'int'"
            in str(excinfo.value))


def test_literal_init_bad_int():
    with pytest.raises(ValueError) as excinfo:
        Literal("abc", INT_TYPE)
    assert ("A scalar integer literal value must conform to the supported format "
            "('(([+-]?[1-9][0-9]*|0)|(NOT_INITIALISED))') but found 'abc'."
            in str(excinfo.value))
    with pytest.raises(ValueError) as excinfo:
        Literal("012", INT_TYPE)
    assert ("A scalar integer literal value must conform to the supported format "
            "('(([+-]?[1-9][0-9]*|0)|(NOT_INITIALISED))') but found '012'."
            in str(excinfo.value))

def test_literal_init_bad_float():
    with pytest.raises(ValueError) as excinfo:
        Literal("123.45F19", FLOAT_TYPE)
    assert("A scalar real literal value must conform to the supported "
           "format" in str(excinfo.value))

def test_literal_init_bad_bool():
    with pytest.raises(ValueError) as excinfo:
        Literal("notTrue", BOOL_TYPE)
    assert("A scalar boolean literal can only be: 'true' or 'false' "
           "but found 'notTrue'." in str(excinfo.value))

def test_valid_literals():
    Literal("true", BOOL_TYPE)
    Literal("false", BOOL_TYPE)

    Literal("1", INT_TYPE)
    Literal("-1", INT_TYPE)
    Literal("12345", INT_TYPE)
    Literal("67890", INT_TYPE)

    Literal("1", FLOAT_TYPE)
    x = Literal("1E10", FLOAT_TYPE)
    assert x.value == "1e10"
    Literal("1.234", FLOAT_TYPE)
    Literal("-0.255", FLOAT_TYPE)
    Literal("123.45678e12345", FLOAT_TYPE)
    Literal("-12.985e-13", FLOAT_TYPE)

def test_literal_nodestr():
    x = Literal("1e10", FLOAT_TYPE)
    correct = "Literal[value:'1e10', Scalar<REAL, SINGLE>]"
    assert correct == x.node_str()
