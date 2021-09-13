import pytest
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSLExceptions import SingletonInstanceError

def test_single_instance():
    a = HartreeParticleDSL._HartreeParticleDSL()
    with pytest.raises(SingletonInstanceError) as excinfo:
        b = HartreeParticleDSL._HartreeParticleDSL()
    assert "Only one instance of _HartreeParticleDSL is allowed" in str(excinfo.value)

def test_config_init():
    a = HartreeParticleDSL.Config()
    assert a.config_type['space'] == {'type' : 'struct space_type', 'is_array' : False }
