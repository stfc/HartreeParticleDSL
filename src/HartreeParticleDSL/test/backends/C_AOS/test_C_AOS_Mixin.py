from HartreeParticleDSL.backends.C_AOS.C_AOS_IO_Mixin import C_AOS_IO_Mixin
import pytest

def test_C_AOS_IO_Mixin():
    '''Tests the abstract methods required for the mixin class throw
       exceptions if not defined'''
    class mixin(C_AOS_IO_Mixin):
        pass

    a = mixin()
    with pytest.raises(NotImplementedError) as excinfo:
        a.gen_code_c(None)
    assert "mixin does not implement required function gen_code_c" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_input_c(None, None)
    assert "mixin does not implement required function call_input_c" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.call_output_c(None, None)
    assert "mixin does not implement required function call_output_c" in str(excinfo.value)
    with pytest.raises(NotImplementedError) as excinfo:
        a.get_includes_c()
    assert "mixin does not implement required function get_includes_c" in str(excinfo.value)

