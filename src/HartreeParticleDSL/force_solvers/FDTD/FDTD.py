from HartreeParticleDSL.force_solvers.generic_force_solver.force_solver import force_solver
from HartreeParticleDSL.HartreeParticleDSL import get_backend

PERIODIC = 1
OTHER=2

TOPHAT=1
OTHER=2

class IllegalInterpolatorError(Exception):
    pass

class FDTD(force_solver):
    def __init__(self, x_min, x_max, y_min=0.0, y_max=0.0,
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
        self.dimensionality = dimensionality
        self.x_bc_max = x_bc_max
        self.x_bc_min = x_bc_min
        self.y_bc_max = y_bc_max
        self.y_bc_min = y_bc_min
        self.z_bc_max = z_bc_max
        self.z_bc_min = z_bc_min
        self.interpolator = TOPHAT

    def set_interpolator(self, interpolator):
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def call_init_grid(self):
        pass

    def call_interpolate_to_particles(self):
        # Tied to C still - better this is some standalone functionality/library system but don't have
        # an idea for that yet.
        # Best idea so far is to have a function call, which takes in pointers to the relevant particle values
        # (i.e. position(s), pressures, mass, charge and EM values) and pointers to the grid, plus the grid 
        # spacing and timestep. The function call should update the relevant values using the pointers supplied
        # and silently finish
        assert self.dimensionality == 1
        assert self.interpolator == TOPHAT
        backend = get_backend()
        code = ""
        #Setup some values
        idx = "double idx = 1.0 / dx;\n"
        code = code + idx
        dto2 = "double dto2 = dt / 2.0;\n"
        code = code + dto2
        dtco2 = "double dtco2 = c * dto2;\n"
        code = code + dtco2
        # Loop over each particle
        loop = backend.per_particle_loop_start("part")
        code = code + loop
        # Setup some values
        part_mc = "double part_mc = {} * c;\n".format(backend.particle_access("part", "mass"))
        code = code + part_mc
        ipart_mc = "double ipart_mc = 1.0 / part_mc;\n"
        code = code + ipart_mc
        part_x = "double part_x = {};\n".format(backend.particle_access("part", "core_part.position[0]"))
        code = code + part_x
        part_ux = "double part_ux = {} * ipart_mc;\n".format(backend.particle_access("part", "part_p[0]"))
        code = code + part_ux
        part_uy = "double part_uy = {} * ipart_mc;\n".format(backend.particle_access("part", "part_p[1]"))
        code = code + part_uy
        part_uz = "double part_uz = {} * ipart_mc;\n".format(backend.particle_access("part", "part_p[2]"))
        code = code + part_uz
        # Calculate v(t) from p(t)
        gamma_rel = "double gamma_rel = sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0);\n"
        code = code + gamma_rel
        root = "double root = dtco2 / gamma_rel;\n"
        code = code + root

        # Move particles to half timestep position
        code = code + "part_x = part_x + part_ux * root;\n"

        # Ignore WORK_DONE_INTEGRATED for now...
        # Grid position as a fraction

        
        pass

    def call_transfer_to_grid(self):
        pass
