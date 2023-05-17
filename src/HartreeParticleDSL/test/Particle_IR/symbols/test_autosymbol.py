import pytest

from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from psyclone.psyir.nodes import Node

def test_autosymbol():
    with pytest.raises(TypeError) as excinfo:
        AutoSymbol(name=32, initial_value="0")
    assert ("<class 'HartreeParticleDSL.Particle_IR.symbols.autosymbol.AutoSymbol'>"
            " 'name' attribute should be of type str "
            "but <class 'int'> found.") in str(excinfo.value)
    
    with pytest.raises(TypeError) as excinfo:
        AutoSymbol(name="scalartype1", initial_value=12)
    assert ("The initial_value of a <class 'HartreeParticleDSL.Particle_IR"
            ".symbols.autosymbol.AutoSymbol'> must be specified using "
            "a str but got <class 'int'>.") in str(excinfo.value)



    autosym = AutoSymbol(name="autosymbol1", initial_value="0")
    assert autosym.initial_value == "0"
