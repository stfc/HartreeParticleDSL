import pytest

from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType
from HartreeParticleDSL.Particle_IR.datatypes.datatype import DOUBLE_TYPE, FLOAT_TYPE

def test_structuresymbol():
    structuretype = StructureType()
    with pytest.raises(TypeError) as excinfo:
        StructureSymbol(name=32, datatype=structuretype)
    assert ("<class 'HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol'>"
            " 'name' attribute should be of type str "
            "but <class 'int'> found.") in str(excinfo.value)
    
    with pytest.raises(TypeError) as excinfo:
        StructureSymbol(name="structure1", datatype=12)
    assert ("The datatype of a <class 'HartreeParticleDSL.Particle_IR"
            ".symbols.structuresymbol.StructureSymbol'> must be specified using "
            "a StructureType but got <class 'int'>.") in str(excinfo.value)



    structure = StructureSymbol(name="structure1", datatype=structuretype)
    assert structure.datatype is structuretype
    structuretype2 = StructureType()
    structure.datatype = structuretype2
    assert structure.datatype is structuretype2

    assert structure.name == "structure1"
    assert structure.visibility == Symbol.Visibility.LOCAL
