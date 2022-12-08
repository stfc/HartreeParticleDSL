import pytest
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.coupled_systems.FDTD_MPI_Kokkos.FDTD_MPI_Kokkos import FDTD_MPI_Kokkos, TOPHAT, PERIODIC, IllegalInterpolatorError
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
from HartreeParticleDSL.Particle_IR.datatypes.datatype import reset_part_and_config, reset_type_mapping_str

def test_FDTD_MPI_Kokkos():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_part_and_config()
    reset_type_mapping_str()
    assert a.dimensionality == 1
    assert a.x_bc_max == PERIODIC
    assert a.x_bc_min == PERIODIC
    assert a.y_bc_max == PERIODIC
    assert a.z_bc_min == PERIODIC
    assert a.y_bc_max == PERIODIC
    assert a.z_bc_min == PERIODIC
    
    correct = '''StructureType<(field: StructureType<(hdt: Scalar<FLOAT, DOUBLE>), (hdtx: Scalar<FLOAT, DOUBLE>), (cnx: Scalar<FLOAT, DOUBLE>), (fac: Scalar<FLOAT, DOUBLE>), (field_order: Scalar<INTEGER, SINGLE>), (fng: Scalar<FLOAT, DOUBLE>), (cfl: Scalar<FLOAT, DOUBLE>), (x_grid_min_local: Scalar<FLOAT, DOUBLE>), (x_grid_max_local: Scalar<FLOAT, DOUBLE>), (x_min_local: Scalar<FLOAT, DOUBLE>), (x_max_local: Scalar<FLOAT, DOUBLE>)>), (nxglobal: Scalar<INTEGER, SINGLE>), (nx: Scalar<INTEGER, SINGLE>), (min_local_cell: Scalar<INTEGER, SINGLE>), (max_local_cell: Scalar<INTEGER, SINGLE>), (ng: Scalar<INTEGER, SINGLE>), (jng: Scalar<INTEGER, SINGLE>), (dx: Scalar<FLOAT, DOUBLE>)>'''
    assert str(backend.structures["field"]) == correct

    assert "\"FDTD_MPI_field.hpp\"" in a._includes
    assert "\"FDTD_MPI_init.hpp\"" in a._includes
    assert "\"FDTD_MPI_boundaries.hpp\"" in a._includes
    assert "\"FDTD_MPI_interpolation.hpp\"" in a._includes
    assert "\"FDTD_MPI_step.hpp\"" in a._includes
    assert "\"FDTD_MPI_init.hpp\"" in a._includes
    assert "\"FDTD_MPI_IO_HDF5.hpp\"" in a._includes

    assert "\"FDTD_MPI_field.hpp\"" in a._includes_header

def test_FDTD_MPI_Kokkos_get_includes():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    includes = a.get_includes()
    reset_part_and_config()
    reset_type_mapping_str()
    assert "\"FDTD_MPI_field.hpp\"" in includes
    assert "\"FDTD_MPI_init.hpp\"" in includes
    assert "\"FDTD_MPI_boundaries.hpp\"" in includes
    assert "\"FDTD_MPI_interpolation.hpp\"" in includes
    assert "\"FDTD_MPI_step.hpp\"" in includes
    assert "\"FDTD_MPI_init.hpp\"" in includes
    assert "\"FDTD_MPI_IO_HDF5.hpp\"" in includes

def test_FDTD_MPI_Kokkos_get_includes_header():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    includes = a.get_includes_header()
    reset_part_and_config()
    reset_type_mapping_str()
    assert "\"FDTD_MPI_field.hpp\"" in includes

def test_FDTD_MPI_Kokkos_set_interpolator():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    a.interpolator = None
    a.set_interpolator(TOPHAT)
    reset_type_mapping_str()
    reset_part_and_config()
    assert a.interpolator == TOPHAT
    with pytest.raises(IllegalInterpolatorError) as excinfo:
        a.set_interpolator(123345)
    assert ("Unsupported interpolator, only TOPHAT is currently supported"
            in str(excinfo.value))

def test_FDTD_MPI_Kokkos_set_interpolator():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    assert a.call_init_grid() == ""

def test_FDTD_MPI_Kokkos_setup_testcase():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.setup_testcase("Myfile.hdf5")
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''field.ng = 4;
field.jng = 4;
load_grid_hdf5(field, "Myfile.hdf5", myrank, nranks, config.config_host(0).space.box_dims);
update_e_field_functor _efield_func(field, field.nx);
update_b_field_functor _bfield_func(field, field.nx);
auto _rp = Kokkos::RangePolicy<>(0, field.nx + 2 * field.ng);
'''
    assert out == correct

def test_FDTD_MPI_Kokkos_call_cleanup_grid():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_cleanup_grid()
    reset_type_mapping_str()
    reset_part_and_config()
    assert out == "kokkos_fdtd_cleanup_1D(field);\n"

def test_FDTD_MPI_Kokkos_call_finish_initialisation_grid():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_finish_initialisation_grid()
    reset_type_mapping_str()
    reset_part_and_config()
    assert out == "bfield_final_bcs(field.bx, field.by, field.bz, field.nx, field.ng);\n"

def test_FDTD_MPI_Kokkos_call_eb_field_first_halfstep():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_eb_fields_first_halfstep()
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''update_eb_fields_half_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,
                         _efield_func, _bfield_func, _rp);
'''
    assert out == correct

def test_FDTD_MPI_Kokkos_call_eb_fields_final_halfstep():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_eb_fields_final_halfstep()
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''update_eb_fields_final_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,
                          _efield_func, _bfield_func, _rp);
'''
    assert correct == out

def test_FDTD_MPI_Kokkos_call_reset_current():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_reset_current()
    reset_type_mapping_str()
    reset_part_and_config()
    assert out == "current_start(field, field.nx, field.ng);\n"

def test_FDTD_MPI_Kokkos_call_finish_current():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.call_finish_current()
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''Kokkos::Experimental::contribute(field.jx, field.scatter_jx);
Kokkos::Experimental::contribute(field.jy, field.scatter_jy);
Kokkos::Experimental::contribute(field.jz, field.scatter_jz);
field.scatter_jx.reset();
field.scatter_jy.reset();
field.scatter_jz.reset();
current_finish(field.jx, field.jy, field.jz,
               field.nx, field.ng);
'''
    assert out == correct
