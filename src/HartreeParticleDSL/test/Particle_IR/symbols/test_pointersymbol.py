import pytest

from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import PointerType
from HartreeParticleDSL.Particle_IR.datatypes.datatype import DOUBLE_TYPE, FLOAT_TYPE

def test_pointersymbol():
    pointertype = PointerType(DOUBLE_TYPE)
    with pytest.raises(TypeError) as excinfo:
        PointerSymbol(name=32, datatype=pointertype)
    assert ("<class 'HartreeParticleDSL.Particle_IR.symbols.pointersymbol.PointerSymbol'>"
            " 'name' attribute should be of type str "
            "but <class 'int'> found.") in str(excinfo.value)
    
    with pytest.raises(TypeError) as excinfo:
        PointerSymbol(name="pointer1", datatype=12)
    assert ("The datatype of a <class 'HartreeParticleDSL.Particle_IR"
            ".symbols.pointersymbol.PointerSymbol'> must be specified using "
            "a PointerType but got <class 'int'>.") in str(excinfo.value)



    pointer = PointerSymbol(name="pointer1", datatype=pointertype)
    assert pointer.datatype is pointertype
    pointertype2 = PointerType(FLOAT_TYPE)
    pointer.datatype = pointertype2
    assert pointer.datatype is pointertype2

    assert pointer.name == "pointer1"
    assert pointer.visibility == Symbol.Visibility.LOCAL
