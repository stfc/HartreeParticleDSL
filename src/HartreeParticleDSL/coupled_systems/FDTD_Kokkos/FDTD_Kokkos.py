from HartreeParticleDSL.coupled_systems.generic_force_solver.force_solver import force_solver
from HartreeParticleDSL.HartreeParticleDSL import get_backend

PERIODIC = 1
OTHER=2

TOPHAT=1
OTHER=2

class IllegalInterpolatorError(Exception):
    pass

class FDTD_Kokkos(force_solver):
    def __init__(self, x_min, x_max, ncells, config_type,
                 y_min=0.0, y_max=0.0,
                 z_min = 0.0, z_max = 0.0, dimensionality = 1,
                 x_bc_max = PERIODIC, x_bc_min = PERIODIC,
                 y_bc_max = PERIODIC, y_bc_min = PERIODIC,
                 z_bc_max = PERIODIC, z_bc_min = PERIODIC):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max
        self.ncells = ncells
        if dimensionality != 1:
            raise NotImplementedError("Only 1D currently supported")
        self.dimensionality = dimensionality
        self.x_bc_max = x_bc_max
        self.x_bc_min = x_bc_min
        self.y_bc_max = y_bc_max
        self.y_bc_min = y_bc_min
        self.z_bc_max = z_bc_max
        self.z_bc_min = z_bc_min
        self.interpolator = TOPHAT
        get_backend().add_type("FDTD_field", "FDTD_field")
        config_type.add_element("field", "FDTD_field")

    def set_interpolator(self, interpolator):
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def call_init_grid(self, current_indent=0, indent=0):
        assert self.dimensionality == 1
        in_str = " " * current_indent
        code = in_str + "config.field.nx = " + str(self.ncells) + ";\n"
        code = code + in_str + "config.field.ng = 4;\n"
        code = code + in_str + "kokkos_fdtd_initialize_1D(config.field, config.field.nx, config.field.ng);\n"
        code = code + in_str + "update_e_field_functor _efield_func(config.field, config.field.nx);\n"
        code = code + in_str + "update_b_field_functor _bfield_func(config.field, config.field.nx);\n"
        code = code + in_str + "auto _rp = Kokkos::RangePolicy<>(0, config.field.nx + 2 * config.field.ng);\n"
        return code

    def setup_testcase(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "kokkos_fdtd_initialize_example_1D(config.field, config.field.nx, config.field.ng);\n"
        return code

    def call_cleanup_grid(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "kokkos_fdtd_cleanup_1D(config.field);\n"
        return code

    def call_eb_fields_first_halfstep(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_half_1D(config.field, config.field.nx, config.field.ng, config.dt, config.dx,\n"
        code = code + in_str + "                         _efield_func, _bfield_func, _rp);\n"
        return code

    def call_eb_fields_final_halfstep(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_final_1D(config.field, config.field.nx, config.field.ng, config.dt, config.dx,\n"
        code = code + in_str + "                          _efield_func, _bfield_func, _rp);\n"
        return code

    def call_reset_current(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = code + in_str + "current_start(config.field, config.field.nx, config.field.ng);\n"
        return code

    def call_interpolate_to_particles(self, part_weight, part_charge, part_mass,
                                      part_momentum_x, part_momentum_y,
                                      part_momentum_z, dx, dt, current_indent=0,
                                      indent=0):
        '''
        '''
        # FIXME
        assert self.dimensionality == 1
        assert self.interpolator == TOPHAT
        backend = get_backend()
        double_type = backend._type_map["c_double"]
        int32_type = backend._type_map["c_int32_t"]
        in_str = " " * current_indent
        # Do all the setup required beforehand
        code = f"{in_str}" + backend.create_variable("c_double", "idt", f"1.0 / {dt}")
        code = code + f"{in_str}" + backend.create_variable("c_double", "idx", f"1.0 / {dx}")
        code = code + f"{in_str}" + backend.create_variable("c_double", "dto2", f"{dt} / 2.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "dtco2", f"c * dto2")
        code = code + f"{in_str}" + backend.create_variable("c_double", "dtfac", f"0.5 * {dt}")
        code = code + f"{in_str}" + backend.create_variable("c_double", "idft", "idt")
        code = code + f"{in_str}" + backend.create_variable("c_double", "idxt", "idx")
        # Can't yet create ararys with create_variable calls.
        code = code + f"{in_str}{double_type} gxarray[4] = " + "{0.0, 0.0, 0., 0.};\n"
        code = code + f"{in_str}{double_type} hxarray[4] = " + "{0.0, 0.0, 0.0, 0.0};\n"
        code = code + f"{in_str}{double_type}* gx = &gxarray[1];\n"
        code = code + f"{in_str}{double_type}* hx = &hxarray[1];\n"
        # Add gx and hx into scope (since we can't currently use create_variable to do this).
        backend.variable_scope.add_variable("gx", "c_double", True)
        backend.variable_scope.add_variable("hx", "c_double", True)
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_weight", backend.get_particle_access("part1", part_weight))
        code = code + f"{in_str}" + backend.create_variable("c_double", "fcx", "idtf * part_weight")
        code = code + f"{in_str}" + backend.create_variable("c_double", "fcy", "idxf * part_weight")
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_q", backend.get_particle_access("part1", part_charge))
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_m", backend.get_particle_access("part1", part_mass))
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_mc", "c * part_m")
        code = code + f"{in_str}" + backend.create_variable("c_double", "ipart_mc", "1.0 / part_mc")
        code = code + f"{in_str}" + backend.create_variable("c_double", "cmratio", "part_q * dtfac * ipart_mc")
        code = code + f"{in_str}" + backend.create_variable("c_double", "ccmratio", "c * cmratio")
        code = code + in_str + "//Copy out the particle properties\n"
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_x",  backend.get_particle_access("part1", backend.get_particle_position("x")) + " - config.field->x_grid_min_local")
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_ux", backend.get_particle_access("part1", part_momentum_x))
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_uy", backend.get_particle_access("part1", part_momentum_y))
        code = code + f"{in_str}" + backend.create_variable("c_double", "part_uz", backend.get_particle_access("part1", part_momentum_z))
        code = code + in_str + "//Calculate v(t) from p(t)\n"
        code = code + f"{in_str}" + backend.create_variable("c_double", "gamma_rel", "sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0)")
        code = code + f"{in_str}" + backend.create_variable("c_double", "root", "dtco2 / gamma_rel") + "\n"
        code = code + in_str + "//Move particles to half timestep position (first order)\n"
        code = code + f"{in_str}part_x = part_x + part_ux * root;\n"
        code = code + f"{in_str}" + backend.create_variable("c_double", "cell_x_r", "part_x * idx - 0.5")
        code = code + f"{in_str}" + backend.create_variable("c_double", "ex_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "ey_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "ez_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "bx_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "by_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_double", "bz_part", "0.0")
        code = code + f"{in_str}" + backend.create_variable("c_int32_t", "cell_x1", "0") + "\n"

        # Call the interpolation function
        code = code + f"{in_str}interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,\n"
        code = code + in_str*2 + backend.get_pointer("part_x") + ", part_p_x, part_p_y,\n"
        code = code + in_str*2 + "part_p_z, " + backend.get_pointer("ex_part") + ",\n"
        code = code + in_str*2 + backend.get_pointer("ey_part") + ", " + backend.get_pointer("ez_part") + ",\n"
        code = code + in_str*2 + backend.get_pointer("bx_part") + ", " + backend.get_pointer("by_part") + ",\n"
        code = code + in_str*2 + backend.get_pointer("bz_part") + ", cell_x_r, " + backend.get_pointer("cell_x1") + ",\n"
        code = code + in_str*2 + "gx, hx, _field.ex, _field.ey, _field.ez, _field.bx, _field.by, _field.bz,\n"
        code = code + in_str*2 + "idt, idx, dtco2, idtf, idxf, _nx, fcx, fcy, _field.ng);\n"
        return code

    def gather_forces_to_grid(self, delta_x, part_vy, part_vz, current_indent=0, indent=0):
        # FIXME
        assert self.dimensionality == 1
        assert self.interpolator == TOPHAT
        backend = get_backend()
        double_type = backend._type_map["c_double"]
        int32_type = backend._type_map["c_int32_t"]
        #Remove excess " from strings
        in_str = " " * current_indent
        code = f"\n{in_str}//Gathering forces to grid\n"
        code = code + f"{in_str}part_x = " + backend.get_particle_access("part1", backend.get_particle_position("x")) + " - config.field->x_grid_min_local;\n"
        code = code + f"{in_str}GatherForcesToGrid_1D(part_weight, part_q, part_x, {delta_x},\n"
        code = code + f"{in_str}{in_str}" + f"cell_x1, _field.jx, _field.jy, _field.jz, idt, {part_vy}, {part_vz}, idx, dtco2, idtf, idxf,\n"
        code = code + f"{in_str}{in_str}" + "_nx, fcx, fcy, _field.ng);\n"
        return code 
