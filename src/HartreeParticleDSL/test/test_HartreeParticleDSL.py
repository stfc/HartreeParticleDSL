import pytest
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSLExceptions import SingletonInstanceError, \
                                                            RepeatedNameError

def test_single_instance():
    '''Test only one instance of HartreeParticleDSL is allowed'''
    a = HartreeParticleDSL._HartreeParticleDSL()
    with pytest.raises(SingletonInstanceError) as excinfo:
        b = HartreeParticleDSL._HartreeParticleDSL()
    assert "Only one instance of _HartreeParticleDSL is allowed" in str(excinfo.value)

def test_config_init():
    '''Test config initialisation'''
    conf = HartreeParticleDSL.Config()
    assert conf.config_type['space'] == {'type' : 'struct space_type', 'is_array' : False }
    assert conf.config_type['neighbour_config'] == {'type': 'struct neighbour_config_type' , 'is_array' : False}

def test_config_add_element():
    '''Test the add_element function of the config'''
    conf = HartreeParticleDSL.Config()
    conf.add_element("value", "double")
    assert conf.config_type["value"] == {'type' : 'double', 'is_array': False}
    conf.add_element("value2", "double[4]")
    assert conf.config_type["value2"] == {'type' : 'double[4]', 'is_array' : True}
    with pytest.raises(RepeatedNameError) as excinfo:
        conf.add_element("value", "int")
    assert "The variable name value is already in the config type" in str(excinfo.value)

def test_config_reset():
    conf = HartreeParticleDSL.Config()
    conf.add_element("value", "double")
    conf.reset_config()
    assert "value" not in conf.config_type.keys()
    assert conf.config_type['space'] == {'type' : 'struct space_type', 'is_array' : False }
    assert conf.config_type['neighbour_config'] == {'type': 'struct neighbour_config_type' , 'is_array' : False}

def test_particle_init():
    '''Test particle initialisation'''
    part = HartreeParticleDSL.Particle()
    part.add_element("value", "double")
    part.reset_particle()
    assert part.particle_type['core_part'] == {'type' : 'struct core_part_type', 'is_array' : False }
    assert part.particle_type['neighbour_part'] == {'type' : 'struct neighbour_part_type', 'is_array' : False}

def test_particle_add_element():
    '''Test the add_element function of the particle'''
    part = HartreeParticleDSL.Particle()
    part.add_element("value", "double")
    assert part.particle_type["value"] == {'type' : 'double', 'is_array': False}
    part.add_element("value2", "double[4]")
    assert part.particle_type["value2"] == {'type' : 'double[4]', 'is_array' : True}
    with pytest.raises(RepeatedNameError) as excinfo:
        part.add_element("value", "int")
    assert "The variable name value is already in the particle type" in str(excinfo.value)
