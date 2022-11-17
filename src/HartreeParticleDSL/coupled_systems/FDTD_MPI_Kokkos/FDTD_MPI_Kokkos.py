from __future__ import annotations

from HartreeParticleDSL.coupled_systems.generic_force_solver.force_solver import force_solver
from HartreeParticleDSL.HartreeParticleDSL import get_backend

from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import Cabana_PIR_Visitor

from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType, DOUBLE_TYPE, INT_TYPE,\
        INT32_TYPE

from HartreeParticleDSL.Particle_IR.nodes.particle_position_reference import ParticlePositionReference

PERIODIC = 1
OTHER=2

TOPHAT=1
OTHER=2

class IllegalInterpolatorError(Exception):
    pass

class FDTD_MPI_Kokkos(force_solver):
    def __init__(self, dimensionality=1,
                 x_bc_max = PERIODIC, x_bc_min = PERIODIC,
                 y_bc_max = PERIODIC, y_bc_min = PERIODIC,
                 z_bc_max = PERIODIC, z_bc_min = PERIODIC):
        self.dimensionality = dimensionality
        self.x_bc_max = x_bc_max
        self.x_bc_min = x_bc_min
        self.y_bc_max = y_bc_max
        self.y_bc_min = y_bc_min
        self.z_bc_max = z_bc_max
        self.z_bc_min = z_bc_min
        self.interpolator = TOPHAT

        # TODO Need to construct the appropriate StructureType here.
        sub_field_struc = StructureType.create([("hdt", DOUBLE_TYPE),
                                                ("hdtx", DOUBLE_TYPE),
                                                ("cnx", DOUBLE_TYPE),
                                                ("fac", DOUBLE_TYPE),
                                                ("field_order", INT_TYPE),
                                                ("fng", DOUBLE_TYPE),
                                                ("cfl", DOUBLE_TYPE),
                                                ("x_grid_min_local", DOUBLE_TYPE),
                                                ("x_grid_max_local", DOUBLE_TYPE),
                                                ("x_min_local", DOUBLE_TYPE),
                                                ("x_max_local", DOUBLE_TYPE)])
       
        # For now, the Kokkos fields are not included in the structure, and are
        # thus inaccessible for the general code, which I think is fine?
        # I might be wrong though.
        fdtd_field = StructureType.create([("field", sub_field_struc),
                                           ("nxglobal", INT_TYPE),
                                           ("nx", INT_TYPE),
                                           ("min_local_cell", INT_TYPE),
                                           ("max_local_cell", INT_TYPE),
                                           ("ng", INT_TYPE),
                                           ("jng", INT_TYPE),
                                           ("dx", DOUBLE_TYPE)])

        get_backend().add_type("FDTD_field", fdtd_field)
        get_backend().add_structure(fdtd_field, "field")

        self._includes = []
        self._includes.append("\"FDTD_MPI_field.hpp\"")
        self._includes.append("\"FDTD_MPI_init.hpp\"")
        self._includes.append("\"FDTD_MPI_boundaries.hpp\"")
        self._includes.append("\"FDTD_MPI_interpolation.hpp\"")
        self._includes.append("\"FDTD_MPI_step.hpp\"")
        self._includes.append("\"FDTD_MPI_init_cabana.hpp\"")
        self._includes.append("\"FDTD_MPI_IO_HDF5.hpp\"")

        self._includes_header = []
        self._includes_header.append("\"FDTD_MPI_field.hpp\"")

    def get_includes(self):
        return self._includes

    def get_includes_header(self):
        return self._includes_header

    def set_interpolator(self, interpolator):
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def call_init_grid(self, current_indent=0, indent=0):
        assert self.dimensionality == 1
        return ""

    def setup_testcase(self, filename, current_indent=0, indent=0):
        in_str = " " * current_indent
        # Can filename be a variable? Need to check
        code = in_str + "load_grid(field, \"" + f"{filename}" + "\", myrank, nranks, box);\n"
        # Letting the particle IO load the particles for now.
        code = code + in_str + "update_e_field_functor _efield_func(field, field.nx);\n"
        code = code + in_str + "update_b_field_functor _bfield_func(field, field.nx);\n"
        code = code + in_str + "auto _rp = Kokkos::RangePolicy<>(0, field.nx + 2 * field.ng);\n"
        return code

    def call_cleanup_grid(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "kokkos_fdtd_cleanup_1D(field);\n"
        return code

    def call_eb_fields_first_halfstep(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_half_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,\n"
        code = code + in_str + "                         _efield_func, _bfield_func, _rp);\n"
        return code

    def call_eb_fields_final_halfstep(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_final_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,\n"
        code = code + in_str + "                          _efield_func, _bfield_func, _rp);\n"
        return code

    def call_reset_current(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "current_start(field, field.nx, field.ng);\n"
        return code

    def call_finish_current(self, current_indent=0, indent=0):
        in_str = " " * current_indent
        code = in_str + "current_finish(field.jx, field.jy, field.jz,\n"
        code = code + in_str + "               field.nx, field.ng);\n"
        return code

    def output_grid(self, filename, variable=None, current_indent=0, indent=0):
        code = "{\n        "
        if variable is not None:
            code = code + "char filename[300];\n"
            code = code + "        sprintf(filename, \"" + f"{filename}%.4d.hdf5" + "\", " + f"{variable});\n        "
        else:
            code = code + "char filename[300]" + f" = \"{filename}\";\n"
        code = code + f"        grid_hdf5_output( field, filename, myrank, nranks);\n"
        code = code + "        }\n"
        return code

    def get_extra_symbols(self, function_list):
        symbols = []
        if "call_interpolate_to_particles" in function_list:
            symbols.append( ("idt", ScalarTypeSymbol("idt", DOUBLE_TYPE)) )
            symbols.append( ("idx", ScalarTypeSymbol("idx", DOUBLE_TYPE)) )
            symbols.append( ("dto2", ScalarTypeSymbol("dto2", DOUBLE_TYPE)) )
            symbols.append( ("dtco2", ScalarTypeSymbol("dtco2", DOUBLE_TYPE)) )
            symbols.append( ("dtfac", ScalarTypeSymbol("dtfac", DOUBLE_TYPE)) )
            symbols.append( ("idtf", ScalarTypeSymbol("idtf", DOUBLE_TYPE)) )
            symbols.append( ("idxf", ScalarTypeSymbol("idxf", DOUBLE_TYPE)) )
            symbols.append( ("part_weight", ScalarTypeSymbol("part_weight", DOUBLE_TYPE)) )
            symbols.append( ("fcx", ScalarTypeSymbol("fcx", DOUBLE_TYPE)) )
            symbols.append( ("fcy", ScalarTypeSymbol("fcy", DOUBLE_TYPE)) )
            symbols.append( ("part_q", ScalarTypeSymbol("part_q", DOUBLE_TYPE)) )
            symbols.append( ("part_m", ScalarTypeSymbol("part_m", DOUBLE_TYPE)) )
            symbols.append( ("part_mc", ScalarTypeSymbol("part_mc", DOUBLE_TYPE)) )
            symbols.append( ("ipart_mc", ScalarTypeSymbol("ipart_mc", DOUBLE_TYPE)) )
            symbols.append( ("cmratio", ScalarTypeSymbol("cmratio", DOUBLE_TYPE)) )
            symbols.append( ("ccmratio", ScalarTypeSymbol("ccmratio", DOUBLE_TYPE)) )
            symbols.append( ("part_x", ScalarTypeSymbol("part_x", DOUBLE_TYPE)) )
            symbols.append( ("part_p_x", ScalarTypeSymbol("part_p_x", DOUBLE_TYPE)) )
            symbols.append( ("part_p_y", ScalarTypeSymbol("part_p_y", DOUBLE_TYPE)) )
            symbols.append( ("part_p_z", ScalarTypeSymbol("part_p_z", DOUBLE_TYPE)) )
            symbols.append( ("part_ux", ScalarTypeSymbol("part_ux", DOUBLE_TYPE)) )
            symbols.append( ("part_uy", ScalarTypeSymbol("part_uy", DOUBLE_TYPE)) )
            symbols.append( ("part_uz", ScalarTypeSymbol("part_uz", DOUBLE_TYPE)) )
            symbols.append( ("gamma_rel", ScalarTypeSymbol("gamma_rel", DOUBLE_TYPE)) )
            symbols.append( ("root", ScalarTypeSymbol("root", DOUBLE_TYPE)) )
            symbols.append( ("cell_x_r", ScalarTypeSymbol("cell_x_r", DOUBLE_TYPE)) )
            symbols.append( ("ex_part", ScalarTypeSymbol("ex_part", DOUBLE_TYPE)) )
            symbols.append( ("ey_part", ScalarTypeSymbol("ey_part", DOUBLE_TYPE)) )
            symbols.append( ("ez_part", ScalarTypeSymbol("ez_part", DOUBLE_TYPE)) )
            symbols.append( ("bx_part", ScalarTypeSymbol("bx_part", DOUBLE_TYPE)) )
            symbols.append( ("by_part", ScalarTypeSymbol("by_part", DOUBLE_TYPE)) )
            symbols.append( ("bz_part", ScalarTypeSymbol("bz_part", DOUBLE_TYPE)) )
            symbols.append( ("cell_x1", ScalarTypeSymbol("cell_x1", INT32_TYPE)) )

        return symbols


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
        code = f"{in_str}" + f"idt = 1.0 / {dt};\n"
        code = code + f"{in_str}" + f"idx = 1.0 / {dx};\n"
        code = code + f"{in_str}" + f"dto2 = {dt} / 2.0;\n"
        code = code + f"{in_str}" + f"dtco2 = c * dto2;\n"
        code = code + f"{in_str}" + f"dtfac = 0.5 * {dt};\n"
        code = code + f"{in_str}" + f"idtf = idt;\n"
        code = code + f"{in_str}" + f"idxf = idx;\n"
        # Can't yet create ararys with create_variable calls.
        code = code + f"{in_str}{double_type} gxarray[4] = " + "{0.0, 0.0, 0., 0.};\n"
        code = code + f"{in_str}{double_type} hxarray[4] = " + "{0.0, 0.0, 0.0, 0.0};\n"
        code = code + f"{in_str}{double_type}* gx = &gxarray[1];\n"
        code = code + f"{in_str}{double_type}* hx = &hxarray[1];\n"
        # Add gx and hx into scope (since we can't currently use create_variable to do this).
#        backend.variable_scope.add_variable("gx", "c_double", True)
#        backend.variable_scope.add_variable("hx", "c_double", True)

        # To get the correct particle access we need to go through the Particle_IR and generate the code. For now
        # its fixed to Cabana_PIR but we should potentially decouple this?

        # Pull the current kernel from the backend so we can grab the symbol table
        kernel = backend.get_current_kernel()
        sym_tab = kernel.symbol_table
        particle_symbols = sym_tab.find_particle_symbols()
        # We assume the first particle we find is particle 1
        part1_sym = list(particle_symbols.items())[0]

        # create the visitor we need
        c_vis = Cabana_PIR_Visitor(backend)
        # Nasty workaround
        c_vis._in_kernel = True

        code = code + f"{in_str}" + f"part_weight = " + part_weight + ";\n"
        code = code + f"{in_str}" + f"fcx = idtf * part_weight;" + "\n"
        code = code + f"{in_str}" + f"fcy = idxf * part_weight;" + "\n"

        code = code + f"{in_str}" + f"part_q = " + part_charge + ";\n"
        code = code + f"{in_str}" + f"part_m = " + part_mass + ";\n"
        code = code + f"{in_str}" + f"part_mc = c * part_m;" + "\n"
        code = code + f"{in_str}" + f"ipart_mc = 1.0 / part_mc;" + "\n"
        code = code + f"{in_str}" + f"cmratio = part_q * dtfac * ipart_mc;" + "\n"
        code = code + f"{in_str}" + f"ccmratio = c * cmratio;" + "\n"
        code = code + in_str + "//Copy out the particle properties\n"
        x_pos_ref = ParticlePositionReference(part1_sym[1], 0)
        code = code + f"{in_str}" + f"part_x = " + c_vis(x_pos_ref) + " - field.field.x_grid_min_local;\n"

        code = code + f"{in_str}" + f"part_p_x = " + part_momentum_x + ";\n"
        code = code + f"{in_str}" + f"part_p_y = " + part_momentum_y + ";\n"
        code = code + f"{in_str}" + f"part_p_z = " + part_momentum_z + ";\n"

        code = code + f"{in_str}" + f"part_ux = " + part_momentum_x + " * ipart_mc;\n"
        code = code + f"{in_str}" + f"part_uy = " + part_momentum_y + " * ipart_mc;\n"
        code = code + f"{in_str}" + f"part_uz = " + part_momentum_z + " * ipart_mc;\n"
        code = code + in_str + "//Calculate v(t) from p(t)\n"
        code = code + f"{in_str}" + f"gamma_rel = sqrtf(part_ux*part_ux + part_uy*part_uy + part_uz*part_uz + 1.0);" + "\n"
        code = code + f"{in_str}" + f"root = dtco2 / gamma_rel;" + "\n\n"
        code = code + in_str + "//Move particles to half timestep position (first order)\n"
        code = code + f"{in_str}part_x = part_x + part_ux * root;\n"

        code = code + f"{in_str}" + f"cell_x_r = part_x * idx - 0.5;" + "\n"
        code = code + f"{in_str}" + f"ex_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"ey_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"ez_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"bx_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"by_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"bz_part = 0.0;" + "\n"
        code = code + f"{in_str}" + f"cell_x1 = 0;" + "\n"

        # TODO
        # Should we add the symbols we "created" here to the symbol table?
        # Need to check if the symbol checking is in Particle_IR to output or
        # in the AST -> Particle_IR stage.

        # Call the interpolation function
        code = code + f"{in_str}interpolate_from_grid_tophat_1D(part_weight, part_q, part_m,\n"
        code = code + in_str*2 + "&part_x" + ", part_p_x, part_p_y,\n"
        code = code + in_str*2 + "part_p_z, " + "&ex_part" + ",\n"
        code = code + in_str*2 + "&ey_part" + ", " + "&ez_part" + ",\n"
        code = code + in_str*2 + "&bx_part" + ", " + "&by_part" + ",\n"
        code = code + in_str*2 + "&bz_part" + ", cell_x_r, " + "&cell_x1" + ",\n"
        code = code + in_str*2 + "gx, hx, field.ex, field.ey, field.ez, field.bx, field.by, field.bz,\n"
        code = code + in_str*2 + "idt, idx, dtco2, idtf, idxf, field.nx, fcx, fcy, field.ng);\n"
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

        # Pull the current kernel from the backend so we can grab the symbol table
        kernel = backend.get_current_kernel()
        sym_tab = kernel.symbol_table
        particle_symbols = sym_tab.find_particle_symbols()
        # We assume the first particle we find is particle 1
        part1_sym = list(particle_symbols.items())[0]
        x_pos_ref = ParticlePositionReference(part1_sym[1], 0)
        # create the visitor we need
        c_vis = Cabana_PIR_Visitor(backend)
        # Nasty workaround
        c_vis._in_kernel = True

        code = code + f"{in_str}part_x = " + c_vis(x_pos_ref) + " - field.field.x_grid_min_local;\n"
        code = code + f"{in_str}GatherForcesToGrid_1D(part_weight, part_q, part_x, {delta_x},\n"
        code = code + f"{in_str}{in_str}" + f"cell_x1, gx, hx, field.jx, field.jy, field.jz, idt, {part_vy}, {part_vz}, idx, dtco2, idtf, idxf,\n"
        code = code + f"{in_str}{in_str}" + "field.nx, fcx, fcy, field.ng);\n"
        return code

    def has_preferred_decomposition(self):
        return True

    def get_preferred_decomposition(self, field_str, current_indent=0, indent=0) -> str:
        '''
        Returns the call to the coupled system which places the preferred
        domain decomposition into the system's boundary object.
        '''
        return " " * current_indent + f"store_domain_decomposition({field_str}, box);\n"
