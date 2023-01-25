import pytest

from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from HartreeParticleDSL.Particle_IR.datatypes.datatype import FLOAT_TYPE, INT_TYPE, BOOL_TYPE

def test_literal_init_bad_inputs():

    with pytest.raises(TypeError) as excinfo:
        Literal(123, INT_TYPE)
    assert ("Literal value must be a string but <class 'int'> supplied." in
            str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        Literal("123", 123)
    assert("Literal datatype must be a ScalarType but <class 'int'> "
           "supplied." in str(excinfo.value))

def test_literal_init_bad_int():
    with pytest.raises(ValueError) as excinfo:
        Literal("abc", INT_TYPE)
    assert("Constructing integer Literal but got a value of 'abc' instead "
           "of an integer value." in str(excinfo.value))
    with pytest.raises(ValueError) as excinfo:
        Literal("012", INT_TYPE)
    assert("Constructing integer Literal but got a value of '012' instead "
           "of an integer value." in str(excinfo.value))

def test_literal_init_bad_float():
    with pytest.raises(ValueError) as excinfo:
        Literal("123.45F19", FLOAT_TYPE)
    assert("Constructing float Literal but got a value of '123.45F19' "
           "instead of a float value." in str(excinfo.value))

def test_literal_init_bad_bool():
    with pytest.raises(ValueError) as excinfo:
        Literal("notTrue", BOOL_TYPE)
    assert("Constructing boolean Literal but got a value of 'notTrue' "
           "instead of True or False." in str(excinfo.value))

def test_valid_literals():
    Literal("True", BOOL_TYPE)
    Literal("False", BOOL_TYPE)

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
    correct = "Literal['1e10', Scalar<FLOAT, SINGLE>]"
    assert correct == x.node_str()
