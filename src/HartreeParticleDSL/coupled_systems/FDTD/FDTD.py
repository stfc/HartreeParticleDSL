from HartreeParticleDSL.force_solvers.generic_force_solver.force_solver import force_solver
from HartreeParticleDSL.HartreeParticleDSL import get_backend

PERIODIC = 1
OTHER=2

TOPHAT=1
OTHER=2

class IllegalInterpolatorError(Exception):
    pass

class FDTD(force_solver):
    def __init__(self, x_min, x_max, ncells, y_min=0.0, y_max=0.0,
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
        if dimensionality != 1:
            raise NotImplementedError("Only 1D currently supported")
        self.dimensionality = dimensionality
        self.x_bc_max = x_bc_max
        self.x_bc_min = x_bc_min
        self.y_bc_max = y_bc_max
        self.y_bc_min = y_bc_min
        self.z_bc_max = z_bc_max
        self.z_bc_min = z_bc_min
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def set_interpolator(self, interpolator):
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def call_init_grid(self):
        pass

    def call_interpolate_to_particles(self, part_weight, part_charge, part_mass, 
                                      part_position, part_momentum_x, part_momentum_y,
                                      part_momentum_z, dx, dt):
        '''
        '''
        assert dimensionality == 1
        assert self.interpolator == TOPHAT
        backend = get_backend()
        # TODO Get indentation
        double_type = backend._type_map["c_double"]
        int32_type = backend._type_map["c_int32_t"]
        # Do all the setup required beforehand
        code = f"{double_type} idt = 1.0 / {dt};\n"
        code = code + f"{double_type} idx = 1.0 / {dx};\n"
        code = code + f"{double_type} dto2 = {dt} / 2.0;\n"
        code = code + f"{double_type} dtco2 = c * dto2;\n"
        code = code + f"{double_type} dtfac = 0.5 * dt;\n"
        code = code + f"{double_type} idtf = idt;\n"
        code = code + f"{double_type} idxf = idx;\n"
        code = code + f"{double_type} gxarray[4] = {0.0, 0.0, 0., 0.};\n"
        code = code + f"{double_type} hxarray[4] = {0.0, 0.0, 0.0, 0.0};\n"
        code = code + f"{double_type}* gx = &gxarray[1];\n"
        code = code + f"{double_type}* hx = &hxarray[1];\n"
        code = code + f"{double_type} part_weight = " + backend.particle_access("part1", part_weight) + ";\n"
        code = code + f"{double_type} fcx = idtf * part_weight;\n"
        code = code + f"{double_type} fcy = idxf * part_weight;\n"
        code = code + f"{double_type} part_q = " + backend.particle_access("part1", part_charge) + ";\n"
        code = code + f"{double_type} part_m = " + backend.particle_access("part1", part_mass) + ";\n"
        code = code + f"{double_type} part_mc = c * part_m;\n"
        code = code + f"{double_type} ipart_mc = 1.0 / part_m;\n"
        code = code + f"{double_type} cmratio = part_q * dtfac * ipart_mc;\n"
        code = code + f"{double_type} ccmratio = c * cmratio;\n"
        code = code + "//Copy out the particle properties\n"
        code = code + f"{double_type} part_x = " + backend.particle_access("part1", backend.get_position("x")) + " - field->x_grid_min_local;\n"
        code = code + f"{double_type} part_ux = " + backend.particle_access("part1", part_momentum_x) + ";\n"
        code = code + f"{double_type} part_uy = " + backend.particle_access("part1", part_momentum_y) + ";\n"
        code = code + f"{double_type} part_uz = " + backend.particle_access("part1", part_momentum_z) + ";\n"
        code = code + "//Calculate v(t) from p(t)\n"
        code = code + f"{double_type} gamma_rel = sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0);\n"
        code = code + f"{double_type} root = dtco2 / gamma_rel;\n\n"
        code = code + "//Move particles to half timestep position (first order)\n"
        code = code + f"part_x = part_x + part_ux * root;\n"
        code = code + f"{double_type} cell_x_r = part_x * idx - 0.5;\n"
        code = code + f"{double_type} ex_part = 0.0;\n"
        code = code + f"{double_type} ey_part = 0.0;\n"
        code = code + f"{double_type} ez_part = 0.0;\n"
        code = code + f"{double_type} bx_part = 0.0;\n"
        code = code + f"{double_type} by_part = 0.0;\n"
        code = code + f"{double_type} bz_part = 0.0;\n"
        code = code + f"{int32_type} cell_x1 = 0;\n\n"

        # Call the interpolation function
        code = code + f"interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,\n"
        code = code + "    " + backend.get_pointer("part_x") + ", part_p_x, part_p_y, \n"
        code = code + "    part_p_z, " + backend.get_pointer("ex_part") + ",\n"
        code = code + "    " + backend.get_pointer("ey_part") + ", " + backend.get_pointer("ez_part") + ",\n"
        code = code + "    " + backend.get_pointer("bx_part") + ", " + backend.get_pointer("by_part") + ",\n"
        code = code + "    " + backend.get_pointer("bz_part") + ", cell_x_r, " + backend.get_pointer("cell_x1") + ",\n"
        code = code + "    field, idt, idx, dtco2, idtf, idxf, field->nx, fcx, fcy);\n"

        return code

    def gather_forces_to_grid(delta_x, part_vy, part_vz):
        assert dimensionality == 1
        assert self.interpolator == TOPHAT
        backend = get_backend()
        # TODO Get indentation
        double_type = backend._type_map["c_double"]
        int32_type = backend._type_map["c_int32_t"]
        code = f"{double_type} part_x = " + backend.particle_access("part1", backend.get_position("x")) + " - field->x_grid_mind_local;\n"
        code = code + f"gather_forces_to_grid_tophat_1D(part_weight, part_q, part_x, {delta_x}, \n"
        code = code + f"    cell_x1, field, idt, {part_vy}, {part_vz}, idx, dtco2, idtf, idxf, \n"
        code = code + f"    field->nx, fcx, fcy);\n"

        return code
