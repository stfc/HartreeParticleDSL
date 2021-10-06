from HartreeParticleDSL.backends.FDPS_backend.FDPS_IO_Mixin import FDPS_IO_Mixin
import pytest

def test_FDPS_IO_Mixin():
    '''Tests the abstract methods required for the mixin class throw
       exceptions if not defined'''
    class mixin(FDPS_IO_Mixin):
        pass

    a = mixin()
    with pytest.raises(NotImplementedError) as excinfo:
        a.gen_code_fdps(None)
    assert "mixin does not implement required function gen_code_fdps" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_input_fdps(None, None)
    assert "mixin does not implement required function call_input_fdps" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_output_fdps(None, None)
    assert "mixin does not implement required function call_output_fdps" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.get_includes_fdps()
    assert "mixin does not implement required function get_includes_fdps" in str(excinfo.value)
