from HartreeParticleDSL.language_utils.variable_scope import *
import pytest
from HartreeParticleDSL.HartreeParticleDSLExceptions import RepeatedNameError, \
                                                            InvalidNameError


def test_variable_access_init():
    array_index = ["0", "1"]
    var = variable("A", "b", False)
    var1_access = variable_access(var)
    assert var1_access.variable is var
    assert var1_access.array_indices == []
    assert var1_access.child is None
    assert var1_access.is_child == False
    var1_access.is_child = True
    assert var1_access.is_child == True
    var1_access.is_child = False

    var2 = variable("B", "c", False)
    var2_access = variable_access(var2, child=None, array_index = array_index)
    assert var2_access.child is None
    var2_access.child = var1_access
    assert var1_access.is_child == True
    assert var2_access.array_indices[0] == "0"
    assert var2_access.array_indices[1] == "1"
    var2_access.add_array_index("2")
    assert var2_access.array_indices[2] == "2"
    assert var2_access.variable is var2

    with pytest.raises(SyntaxError) as excinfo:
        var2_access.add_array_index(12345)
    assert ("The supplied array_index was neither a "
            "variable_access or str, received <class 'int'>") in str(excinfo.value)

def test_variable_init():
    var = variable("A", "b", False)
    assert var.var_name == "A"
    assert var.var_type == "b"
    assert var.is_pointer is False

def test_variable_scope():
    scope = variable_scope()
    scope.add_variable("A", "b", False)
    assert scope.get_variable("A") is not None
    with pytest.raises(RepeatedNameError) as excinfo:
        scope.add_variable("A", "b", True)
    assert "A is already defined in the current scope" in str(excinfo.value)
    assert scope.get_variable("Z") is None
    scope.remove_variable("A")
    assert scope.get_variable("A") is None
    with pytest.raises(InvalidNameError) as excinfo:
        scope.remove_variable("A")
    assert "A is not defined in the current scope" in str(excinfo.value)
