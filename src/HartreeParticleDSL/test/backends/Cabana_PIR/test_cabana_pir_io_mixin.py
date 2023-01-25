from HartreeParticleDSL.backends.Cabana_PIR_backend.Cabana_PIR_IO_Mixin import Cabana_PIR_IO_Mixin
import pytest

def test_Cabana_PIR_IO_Mixin():
    '''Tests the abstract methods required for the mixin class throw
       exceptions if not defined'''
    class mixin(Cabana_PIR_IO_Mixin):
        pass

    a = mixin()
    with pytest.raises(NotImplementedError) as excinfo:
        a.gen_code_cabana_pir(None)
    assert "mixin does not implement required function gen_code_cabana_pir" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_input_cabana_pir(None, None)
    assert "mixin does not implement required function call_input_cabana_pir" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_output_cabana_pir(None, None, None)
    assert "mixin does not implement required function call_output_cabana_pir" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.get_includes_cabana_pir()
    assert "mixin does not implement required function get_includes_cabana_pir" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_get_box_size_pir(None, None)
    assert "mixin does not implement required function call_get_box_size_pir" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.get_header_includes_cabana_pir()
    assert "mixin does not implement required function get_header_includes_cabana_pir" in str(excinfo.value)
