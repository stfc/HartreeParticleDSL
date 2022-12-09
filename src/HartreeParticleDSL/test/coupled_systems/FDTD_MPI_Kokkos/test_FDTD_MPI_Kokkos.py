import pytest

import ast
import inspect
import textwrap
import os

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import part, config
from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import pir_perpart_visitor
from HartreeParticleDSL.coupled_systems.FDTD_MPI_Kokkos.FDTD_MPI_Kokkos import FDTD_MPI_Kokkos, TOPHAT, PERIODIC, IllegalInterpolatorError
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import Cabana_PIR_Visitor
from HartreeParticleDSL.Particle_IR.datatypes.datatype import reset_part_and_config, reset_type_mapping_str

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, ScalarType, \
        StructureType, ArrayType, PointerType, \
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, reset_part_and_config

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

def test_FDTD_MPI_Kokkos_output_grid():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.output_grid("out.hdf5")
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''
{
char filename[300] = "out.hdf5";
grid_hdf5_output( field, filename, myrank, nranks);
}
'''
    assert out == correct

def test_FDTD_MPI_Kokkos_get_extra_symbols():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    out = a.get_extra_symbols([])
    reset_type_mapping_str()
    reset_part_and_config()
    assert len(out) == 0

    out = a.get_extra_symbols(["call_interpolate_to_particles"])
    assert len(out) == 33
    assert out[32][0] == "cell_x1"
    assert out[0][0] == "idt"

def test_FDTD_MPI_Kokkos_call_interpolate_to_particles():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    # Create a perpart_kernel
    v = pir_perpart_visitor()
    type_mapping_str["part"].add("subpart", INT_TYPE)
    new_type = ArrayType(INT_TYPE, [3])
    type_mapping_str["part"].add("subarray", new_type)
    def x(arg: part, c: config):
        arg.subpart = 1
        arg.subarray[0] = 1
        arg.core_part.position[0] = 2.0

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = v.visit(c)
    backend._current_kernel = pir


    a = FDTD_MPI_Kokkos()
    out = a.call_interpolate_to_particles("part_weight", "part_charge", "part_mass", "part_p[0]", "part_p[1]", "part_p[2]", "dx", "dt")
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''idt = 1.0 / dt;
idx = 1.0 / dx;
dto2 = dt / 2.0;
dtco2 = c * dto2;
dtfac = 0.5 * dt;
idtf = idt;
idxf = idx;
double gxarray[4] = {0.0, 0.0, 0., 0.};
double hxarray[4] = {0.0, 0.0, 0.0, 0.0};
double* gx = &gxarray[1];
double* hx = &hxarray[1];
part_weight = part_weight;
fcx = idtf * part_weight;
fcy = idxf * part_weight;
part_q = part_charge;
part_m = part_mass;
part_mc = c * part_m;
ipart_mc = 1.0 / part_mc;
cmratio = part_q * dtfac * ipart_mc;
ccmratio = c * cmratio;
//Copy out the particle properties
part_x = _core_part_position.access(i, a, 0) - field.field.x_grid_min_local;
part_p_x = part_p[0];
part_p_y = part_p[1];
part_p_z = part_p[2];
part_ux = part_p[0] * ipart_mc;
part_uy = part_p[1] * ipart_mc;
part_uz = part_p[2] * ipart_mc;
//Calculate v(t) from p(t)
gamma_rel = sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0);
root = dtco2 / gamma_rel;

//Move particles to half timestep position (first order)
part_x = part_x + part_ux * root;
cell_x_r = part_x * idx - 0.5;
ex_part = 0.0;
ey_part = 0.0;
ez_part = 0.0;
bx_part = 0.0;
by_part = 0.0;
bz_part = 0.0;
cell_x1 = 0;
interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,
&part_x, part_p_x, part_p_y,
part_p_z, &ex_part,
&ey_part, &ez_part,
&bx_part, &by_part,
&bz_part, cell_x_r, &cell_x1,
gx, hx, field.ex, field.ey, field.ez, field.bx, field.by, field.bz,
idt, idx, dtco2, idtf, idxf, field.nx, fcx, fcy, field.ng);
'''
    assert out == correct

def test_FDTD_MPI_Kokkos_gather_forces_to_grid():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    # Create a perpart_kernel
    v = pir_perpart_visitor()
    type_mapping_str["part"].add("subpart", INT_TYPE)
    new_type = ArrayType(INT_TYPE, [3])
    type_mapping_str["part"].add("subarray", new_type)
    def x(arg: part, c: config):
        arg.subpart = 1
        arg.subarray[0] = 1
        arg.core_part.position[0] = 2.0

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = v.visit(c)
    backend._current_kernel = pir

    out = a.gather_forces_to_grid("delta_x", "part_vy", "part_vz")
    reset_type_mapping_str()
    reset_part_and_config()
    correct = '''
//Gathering forces to grid
part_x = _core_part_position.access(i, a, 0) - field.field.x_grid_min_local;
GatherForcesToGrid_1D(part_weight, part_q, part_x, delta_x,
cell_x1, gx, hx, field.scatter_jx, field.scatter_jy, field.scatter_jz, idt, part_vy, part_vz, idx, dtco2, idtf, idxf,
field.nx, fcx, fcy, field.ng);
'''
    assert correct == out

def test_FDTD_MPI_Kokkos_has_preferred_decomposition():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    assert a.has_preferred_decomposition() == True

def test_FDTD_MPI_Kokkos_get_preferred_decomposition():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    assert a.get_preferred_decomposition("box") == "store_domain_decomposition(field, box);\n"

def test_FDTD_MPI_Kokkos_copy_files():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    a.copy_files()
    files = ['FDTD_MPI_IO_HDF5.cpp',
             'FDTD_MPI_IO_HDF5.hpp',
             'FDTD_MPI_boundaries.cpp',
             'FDTD_MPI_boundaries.hpp',
             'FDTD_MPI_field.hpp',
             'FDTD_MPI_init.cpp',
             'FDTD_MPI_init.hpp',
             'FDTD_MPI_interpolation.hpp',
             'FDTD_MPI_step.cpp',
             'FDTD_MPI_step.hpp']
    success = True
    for f in files:
        success = success and os.path.exists(f)
        if os.path.exists(f):
            os.remove(f)
    assert success

def test_FDTD_MPI_Kokkos_compilation_files():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    out = a.compilation_files()
    assert out[0] == "FDTD_MPI_IO_HDF5.cpp"
    assert out[1] == "FDTD_MPI_boundaries.cpp"
    assert out[2] == "FDTD_MPI_init.cpp"
    assert out[3] == "FDTD_MPI_step.cpp"
    assert len(out) == 4

def test_FDTD_MPI_Kokkos_get_required_packages():
    backend = Cabana_PIR()
    HartreeParticleDSL.set_backend(backend)
    a = FDTD_MPI_Kokkos()
    reset_type_mapping_str()
    reset_part_and_config()
    out = a.get_required_packages()
    assert out[0] == "Kokkos"
    assert len(out) == 1
