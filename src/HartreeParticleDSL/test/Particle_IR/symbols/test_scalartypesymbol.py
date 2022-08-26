import pytest

from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType
from HartreeParticleDSL.Particle_IR.datatypes.datatype import DOUBLE_TYPE, FLOAT_TYPE

def test_scalartypesymbol():
    scalartypetype = DOUBLE_TYPE
    with pytest.raises(TypeError) as excinfo:
        ScalarTypeSymbol(name=32, datatype=scalartypetype)
    assert ("<class 'HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol.ScalarTypeSymbol'>"
            " 'name' attribute should be of type str "
            "but <class 'int'> found.") in str(excinfo.value)
    
    with pytest.raises(TypeError) as excinfo:
        ScalarTypeSymbol(name="scalartype1", datatype=12)
    assert ("The datatype of a <class 'HartreeParticleDSL.Particle_IR"
            ".symbols.scalartypesymbol.ScalarTypeSymbol'> must be specified using "
            "a ScalarType but got <class 'int'>.") in str(excinfo.value)



    scalartype = ScalarTypeSymbol(name="scalartype1", datatype=scalartypetype)
    assert scalartype.datatype is scalartypetype
    scalartypetype2 =FLOAT_TYPE
    scalartype.datatype = scalartypetype2
    assert scalartype.datatype is scalartypetype2

    assert scalartype.name == "scalartype1"
    assert scalartype.visibility == Symbol.Visibility.LOCAL
