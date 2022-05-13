from HartreeParticleDSL.backends.Cabana_backend.Cabana_IO_Mixin import Cabana_IO_Mixin
import pytest

def test_Cabana_IO_Mixin():
    '''Tests the abstract methods required for the mixin class throw
       exceptions if not defined'''
    class mixin(Cabana_IO_Mixin):
        pass

    a = mixin()
    with pytest.raises(NotImplementedError) as excinfo:
        a.gen_code_cabana(None)
    assert "mixin does not implement required function gen_code_cabana" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_input_cabana(None, None)
    assert "mixin does not implement required function call_input_cabana" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_output_cabana(None, None)
    assert "mixin does not implement required function call_output_cabana" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.get_includes_cabana()
    assert "mixin does not implement required function get_includes_cabana" in str(excinfo.value)

