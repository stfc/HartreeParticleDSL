import pytest

from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import ArrayType, DOUBLE_TYPE

def test_arraysymbol():
    atype = ArrayType(DOUBLE_TYPE, [3])
    with pytest.raises(TypeError) as excinfo:
        ArraySymbol(name=32, datatype=atype)
    assert ("<class 'HartreeParticleDSL.Particle_IR.symbols.arraysymbol.ArraySymbol'>"
            " 'name' attribute should be of type str "
            "but <class 'int'> found.") in str(excinfo.value)
    
    with pytest.raises(TypeError) as excinfo:
        ArraySymbol(name="array1", datatype=12)
    assert ("The datatype of a <class 'HartreeParticleDSL.Particle_IR"
            ".symbols.arraysymbol.ArraySymbol'> must be specified using "
            "an ArrayType but got <class 'int'>.") in str(excinfo.value)



    array = ArraySymbol(name="array1", datatype=atype)
    assert array.datatype is atype
    atype2 = ArrayType(DOUBLE_TYPE, [4])
    array.datatype = atype2
    assert array.datatype is atype2

    assert array.name == "array1"
    assert array.visibility == Symbol.Visibility.LOCAL
