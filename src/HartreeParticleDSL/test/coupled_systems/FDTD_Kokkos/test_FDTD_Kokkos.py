import pytest
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.coupled_systems.FDTD_Kokkos.FDTD_Kokkos import FDTD_Kokkos, TOPHAT, PERIODIC, IllegalInterpolatorError
from HartreeParticleDSL.backends.Cabana_backend.Cabana import Cabana
from HartreeParticleDSL.backends.C_AOS.C_AOS import C_AOS

def test_FDTD_Kokkos():
    cabana = Cabana()
    HartreeParticleDSL.set_backend(cabana)
    config = HartreeParticleDSL.Config()
    a = FDTD_Kokkos(0.0, 2.0, 25, config, y_min = -1.0, y_max = 2.5, z_min = 0.0, z_max = 1.0)
    assert a.x_min == 0.0
    assert a.x_max == 2.0
    assert a.y_min == -1.0
    assert a.y_max == 2.5
    assert a.z_min == 0.0
    assert a.z_max == 1.0

    assert a.x_bc_max == 1
    assert a.x_bc_min == 1
    assert a.y_bc_max == 1
    assert a.y_bc_min == 1
    assert a.z_bc_max == 1
    assert a.z_bc_min == 1

    # Test get_includes here
    includes = a.get_includes()
    assert "\"FDTD_field.hpp\"" in includes
    assert "\"FDTD_init.hpp\"" in includes
    assert "\"FDTD_interpolation.hpp\"" in includes
    assert "\"FDTD_boundaries.hpp\"" in includes
    assert "\"FDTD_step.hpp\"" in includes
    assert "\"FDTD_init_cabana.hpp\"" in includes
    assert "\"FDTD_IO_HDF5.hpp\"" in includes
    
    # Test get_includes_header here
    header_includes = a.get_includes_header()
    assert "\"FDTD_field.hpp\"" in header_includes

    with pytest.raises(NotImplementedError) as excinfo:
        b = FDTD_Kokkos(0, 2, 25, config, dimensionality=3)
    assert "Only 1D currently supported" in str(excinfo.value)

    a.set_interpolator(TOPHAT)
    assert a.interpolator == TOPHAT

    with pytest.raises(IllegalInterpolatorError) as excinfo:
        a.set_interpolator("Not valid")
    assert "Unsupported interpolator, only TOPHAT is currently supported." in str(excinfo.value)

    # Test call_init_grid
    res = a.call_init_grid()
    correct = '''field.nx = 25;
field.ng = 4;
kokkos_fdtd_initialize_1D(field, field.nx, field.ng);
update_e_field_functor _efield_func(field, field.nx);
update_b_field_functor _bfield_func(field, field.nx);
auto _rp = Kokkos::RangePolicy<>(0, field.nx + 2 * field.ng);
'''
    assert correct == res

    # Test setup_testcase
    res = a.setup_testcase()
    correct = '''kokkos_fdtd_initialize_example_1D(field, field.nx, field.ng);
fdtd_init_particles<decltype(particle_aosoa_host)>(particle_aosoa_host, config, field);
Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
Kokkos::deep_copy(config.config, config.config_host);
'''
    assert correct == res

    # Test call_cleanup_grid
    res = a.call_cleanup_grid()
    correct = "kokkos_fdtd_cleanup_1D(config);\n"
    assert correct == res

    # Test call_eb_fields_first_halfstep()
    res = a.call_eb_fields_first_halfstep()
    correct = '''update_eb_fields_half_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,
                         _efield_func, _bfield_func, _rp);
'''
    assert correct == res

    # Test call_eb_fields_final_halfstep
    res = a.call_eb_fields_final_halfstep()
    correct = '''update_eb_fields_final_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,
                          _efield_func, _bfield_func, _rp);
'''
    assert correct == res

    # Test call_reset_current
    res = a.call_reset_current()
    correct = '''current_start(field, field.nx, field.ng);
'''
    assert correct == res

    # Test call_finish_current
    res = a.call_finish_current()
    correct = '''current_finish(field.jx, field.jy, field.jz,
               field.nx, field.ng);
'''
    assert correct == res

    # Test output_grid
    res = a.output_grid("abc.def")
    correct = '''{
        char filename[300] = "abc.def";
        grid_hdf5_output( field, filename);
        }
'''
    assert correct == res

    res = a.output_grid("abc.def", variable="x")
    correct = '''{
        char filename[300];
        sprintf(filename, "abc.def%.4d.hdf5", x);
                grid_hdf5_output( field, filename);
        }
'''
    assert correct == res

    #Test call_interpolate_to_particles
    res = a.call_interpolate_to_particles("weight", "charge", "mass", "mom_x", "mom_y", "mom_z", "dx", "dt")
    correct = '''double idt = 1.0 / dt;
double idx = 1.0 / dx;
double dto2 = dt / 2.0;
double dtco2 = c * dto2;
double dtfac = 0.5 * dt;
double idtf = idt;
double idxf = idx;
double gxarray[4] = {0.0, 0.0, 0., 0.};
double hxarray[4] = {0.0, 0.0, 0.0, 0.0};
double* gx = &gxarray[1];
double* hx = &hxarray[1];
double part_weight = _weight.access(i, a);
double fcx = idtf * part_weight;
double fcy = idxf * part_weight;
double part_q = _charge.access(i, a);
double part_m = _mass.access(i, a);
double part_mc = c * part_m;
double ipart_mc = 1.0 / part_mc;
double cmratio = part_q * dtfac * ipart_mc;
double ccmratio = c * cmratio;
//Copy out the particle properties
double part_x = _core_part_position.access(i, a, 0) - field.field.x_grid_min_local;
double part_p_x = _mom_x.access(i, a);
double part_p_y = _mom_y.access(i, a);
double part_p_z = _mom_z.access(i, a);
double part_ux = _mom_x.access(i, a)* ipart_mc;
double part_uy = _mom_y.access(i, a)* ipart_mc;
double part_uz = _mom_z.access(i, a)* ipart_mc;
//Calculate v(t) from p(t)
double gamma_rel = sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0);
double root = dtco2 / gamma_rel;

//Move particles to half timestep position (first order)
part_x = part_x + part_ux * root;
double cell_x_r = part_x * idx - 0.5;
double ex_part = 0.0;
double ey_part = 0.0;
double ez_part = 0.0;
double bx_part = 0.0;
double by_part = 0.0;
double bz_part = 0.0;
int32_t cell_x1 = 0;

interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,
&(part_x), part_p_x, part_p_y,
part_p_z, &(ex_part),
&(ey_part), &(ez_part),
&(bx_part), &(by_part),
&(bz_part), cell_x_r, &(cell_x1),
gx, hx, field.ex, field.ey, field.ez, field.bx, field.by, field.bz,
idt, idx, dtco2, idtf, idxf, field.nx, fcx, fcy, field.ng);
'''
    assert correct == res

    # Test gather_forces_to_grid
    res = a.gather_forces_to_grid("del_x", "part_vy", "part_vz")
    correct = '''
//Gathering forces to grid
part_x = _core_part_position.access(i, a, 0) - field.field.x_grid_min_local;
GatherForcesToGrid_1D(part_weight, part_q, part_x, del_x,
cell_x1, gx, hx, field.jx, field.jy, field.jz, idt, part_vy, part_vz, idx, dtco2, idtf, idxf,
field.nx, fcx, fcy, field.ng);
'''
    assert correct == res

    # Reset backend
    HartreeParticleDSL.set_backend(C_AOS())
