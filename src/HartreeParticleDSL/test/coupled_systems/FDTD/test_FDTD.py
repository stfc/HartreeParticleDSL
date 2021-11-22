import pytest
from HartreeParticleDSL.coupled_systems.FDTD.FDTD import *
from HartreeParticleDSL.HartreeParticleDSL import *
from HartreeParticleDSL.backends.C_AOS.C_AOS import *

def test_FDTD_init():
    conf = Config()
    s = FDTD(0.0, 1.0, 256, conf)
    assert s.x_min == 0.0
    assert s.x_max == 1.0
    assert s.ncells == 256
    assert "field" in conf.config_type
    with pytest.raises(NotImplementedError) as excinfo:
        s = FDTD(0.0, 1.0, 256, conf, dimensionality = 3)
    assert "Only 1D currently supported" in str(excinfo.value)

def test_set_interoplator():
    conf = Config()
    s = FDTD(0.0, 1.0, 256, conf)
    s.set_interpolator(TOPHAT)
    assert s.interpolator == TOPHAT
    with pytest.raises(IllegalInterpolatorError) as excinfo:
        s.set_interpolator("stringy")
    assert ("Unsupported interpolator, only TOPHAT is currently supported.") in str(excinfo.value)
    
def test_call_init_grid():
    conf = Config()
    s = FDTD(0.0, 1.0, 256, conf)
    rval = s.call_init_grid(current_indent=4)
    correct = '''    config.field = (struct *FDTD_field) malloc(sizeof(struct FDTD_field));
    fdtd_initialize(config.field, config.nx, config.ng);\n'''
    assert correct == rval

def test_setup_testcase():
    conf = Config()
    s = FDTD(0.0, 1.0, 256, conf)
    rval = s.setup_testcase()
    correct = '''fdtd_initialize_example_1D(config.field, config.nx, config.ng);\n'''
    assert correct == rval

def test_call_cleanup_grid():
    conf = Config()
    s = FDTD(0.0, 1.0, 256, conf)
    rval = s.call_cleanup_grid()
    correct = '''FDTD_cleanup_1D(config.field);
free(config.field);\n'''
    assert correct == rval

def test_call_interpolate_to_particles():
    conf = Config()
    backend = C_AOS()
    backend.disable_variable_checks()
    set_backend(backend)
    s = FDTD(0.0, 1.0, 256, conf)
    rval = s.call_interpolate_to_particles("part_weight", "part_charge", "part_mass",
            "part_momentum_x", "part_momentum_y", "part_momentum_z", "dx", "dt")
    correct = '''double idt = 1.0 / dt;
double idx = 1.0 / dx;
double dto2 = dt / 2.0;
double dtco2 = c * dto2;
double dtfac = 0.5 * dt;
double idft = idt;
double idxt = idx;
double gxarray[4] = {0.0, 0.0, 0., 0.};
double hxarray[4] = {0.0, 0.0, 0.0, 0.0};
double* gx = &gxarray[1];
double* hx = &hxarray[1];
double part_weight = part1->part_weight;
double fcx = idtf * part_weight;
double fcy = idxf * part_weight;
double part_q = part1->part_charge;
double part_m = part1->part_mass;
double part_mc = c * part_m;
double ipart_mc = 1.0 / part_mc;
double cmratio = part_q * dtfac * ipart_mc;
double ccmratio = c * cmratio;
//Copy out the particle properties
double part_x = part1->core_part.position[0] - config.field->x_grid_min_local;
double part_ux = part1->part_momentum_x;
double part_uy = part1->part_momentum_y;
double part_uz = part1->part_momentum_z;
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
int cell_x1 = 0;

interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,
&part_x, part_p_x, part_p_y,
part_p_z, &ex_part,
&ey_part, &ez_part,
&bx_part, &by_part,
&bz_part, cell_x_r, &cell_x1,
config.field, idt, idx, dtco2, idtf, idxf, config.field->nx, fcx, fcy);\n'''
    assert correct in rval

def test_gather_forces_to_grid():
    conf = Config()
    backend = C_AOS()
    backend.disable_variable_checks()
    set_backend(backend)
    s = FDTD(0.0, 1.0, 256, conf)
    rval = s.gather_forces_to_grid("delta_x", "part_vy", "part_vz")
    correct = '''
//Gathering forces to grid
part_x = part1->core_part.position[0] - config.field->x_grid_min_local;
gather_forces_to_grid_tophat_1D(part_weight, part_q, part_x, delta_x,
cell_x1, config.field, idt, part_vy, part_vz, idx, dtco2, idtf, idxf,
config.field->nx, fcx, fcy);\n'''
    print(correct)
    print(rval)
    assert correct == rval
