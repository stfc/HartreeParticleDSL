from __future__ import annotations

from shutil import copy
import inspect
import os

from HartreeParticleDSL.coupled_systems.generic_force_solver.force_solver import force_solver
from HartreeParticleDSL.HartreeParticleDSL import get_backend
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL

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
    '''
    A coupled system supporting Finite Difference Time Domain as part of
    HartreeParticleDSL.

    Currently only supports 1D testcases and periodic boundary conditions.

    All arguments should be left at default for now, but are implemented to
    enable future support

    :param int dimensionality: The dimensionality of the problem, default 1.
    :param x_bc_max: Default periodic.
    :param x_bc_min: Default periodic.
    :param y_bc_max: Default periodic.
    :param y_bc_min: Default periodic.
    :param z_bc_max: Default periodic.
    :param z_bc_min: Default periodic.
    '''
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
        self._includes.append("\"FDTD_MPI_IO_HDF5.hpp\"")

        self._includes_header = []
        self._includes_header.append("\"FDTD_MPI_field.hpp\"")

    def get_includes(self):
        '''
        :returns: The list of includes required for this coupled system.
        :rtype: List of str.
        '''
        return self._includes

    def get_includes_header(self):
        '''
        :returns: The list of includes required for the haeder file for this \
                coupled system.
        :rtype: List of str.
        '''
        return self._includes_header

    def set_interpolator(self, interpolator):
        '''
        Set the interpolator type to be used. This has to be set in the user
        script if they want to use a non-tophat interpolator.

        Currently only tophat is supported anyway.

        :raises IllegalInterpolatorError: If an unsupported interpolator is \
                provided.    
        '''
        if interpolator != TOPHAT:
            raise IllegalInterpolatorError("Unsupported interpolator, only "
                                           "TOPHAT is currently supported.")
        self.interpolator = TOPHAT

    def call_init_grid(self, current_indent=0, indent=0):
        '''
        Call to use from a main script to initialise the grid.

        :returns: the code to call to initialise the grid for simulation.
        :rtype: str
        '''
        return ""

    def setup_testcase(self, filename, current_indent=0, indent=0):
        '''
        Call to use to setup the testcase.

        The filename argument is a required argument and is used to read an
        input file to initialise the grid. The file must be a HDF5 file.

        :param str filename: The HDF5 file to use to setup the testcase.
        '''
        in_str = " " * current_indent
        # For now these are fixed, but will eventually depend on interpolation type
        code = in_str + "field.ng = 4;\n"
        code = code + in_str + "field.jng = 4;\n"
        # Can filename be a variable? Need to check
        if "\"" in filename:
            code = code + in_str + "load_grid_hdf5(field, " + f"{filename}" + ", myrank, nranks, config.config_host(0).space.box_dims);\n"
        else:
            code = code + in_str + "load_grid_hdf5(field, \"" + f"{filename}" + "\", myrank, nranks, config.config_host(0).space.box_dims);\n"
        # Letting the particle IO load the particles for now.
        code = code + in_str + "update_e_field_functor _efield_func(field, field.nx);\n"
        code = code + in_str + "update_b_field_functor _bfield_func(field, field.nx);\n"
        code = code + in_str + "auto _rp = Kokkos::RangePolicy<>(0, field.nx + 2 * field.ng);\n"
        return code

    def call_cleanup_grid(self, current_indent=0, indent=0):
        '''
        Call to use to cleanup the grid solver at the end of the run.
        '''
        in_str = " " * current_indent
        code = in_str + "kokkos_fdtd_cleanup_1D(field);\n"
        return code

    def call_finish_initialisation_grid(self, current_indent=0, indent=0):
        '''
        Final call to initialise the grid, which must be called before
        starting the main section of the simulation.
        '''
        in_str = " " * current_indent
        code = in_str + "bfield_final_bcs(field.bx, field.by, field.bz, field.nx, field.ng);\n"
        return code

    def call_eb_fields_first_halfstep(self, current_indent=0, indent=0):
        '''
        Computes the update to the electric and magnetic field for the first
        halfstep based on the computed current.

        Used in user main function definition.
        '''
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_half_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,\n"
        code = code + in_str + "                         _efield_func, _bfield_func, _rp);\n"
        return code

    def call_eb_fields_final_halfstep(self, current_indent=0, indent=0):
        '''
        Computes the update to the electric and magnetic field for the
        second halfstep based on the computed current.

        Used in user main function definition.
        '''
        in_str = " " * current_indent
        code = in_str + "update_eb_fields_final_1D(field, field.nx, field.ng, config.config_host(0).dt, config.config_host(0).dx,\n"
        code = code + in_str + "                          _efield_func, _bfield_func, _rp);\n"
        return code

    def call_reset_current(self, current_indent=0, indent=0):
        '''
        Resets the current to zero to start a new step.

        Used in user main function definition.
        '''
        in_str = " " * current_indent
        code = in_str + "current_start(field, field.nx, field.ng);\n"
        return code

    def call_finish_current(self, current_indent=0, indent=0):
        '''
        Gathers the current for the current step before any field updates
        occur.

        Used in user main function definition
        '''
        in_str = " " * current_indent
        code = in_str + "Kokkos::Experimental::contribute(field.jx, field.scatter_jx);\n"
        code = code + in_str + "Kokkos::Experimental::contribute(field.jy, field.scatter_jy);\n"
        code = code + in_str + "Kokkos::Experimental::contribute(field.jz, field.scatter_jz);\n"
        code = code + in_str + "field.scatter_jx.reset();\n"
        code = code + in_str + "field.scatter_jy.reset();\n"
        code = code + in_str + "field.scatter_jz.reset();\n"
        code = code + in_str + "current_finish(field.jx, field.jy, field.jz,\n"
        code = code + in_str + "               field.nx, field.ng);\n"
        return code

    def output_grid(self, filename, variable=None, current_indent=0, indent=0):
        '''
        Outputs the current grid system to the file specified in filename as
        a HDF5 file. If the file already exists it will be appended to the
        file.

        If variable is specified the filename will be constructed as
        filename_{variable}.hdf5 where variable is expressed as a 4 length
        integer.

        :param str filename: The filename to output the grid to.
        :param str variable: An optional variable to use as part of the \
                filename construction.
        '''
        code = "\n" + current_indent * " " + "{\n" 
        current_indent = current_indent + indent
        code = code + current_indent * " "
        if variable is not None:
            code = code + "char filename[300];\n"
            if "\"" in filename:
                rep_filename = filename.replace("\"", "")
            else:
                rep_filename = filename
            code = code + current_indent* " " + "sprintf(filename, \"" + f"{rep_filename}%.4d.hdf5" + "\", " + f"{variable});\n"
        else:
            if "\"" in filename:
                code = code + "char filename[300]" + f" = {filename};\n"
            else:
                code = code + "char filename[300]" + f" = \"{filename}\";\n"
        code = code + current_indent * " " +  f"grid_hdf5_output( field, filename, myrank, nranks);\n"
        current_indent = current_indent - indent
        code = code + current_indent * " " + "}\n"
        return code

    def get_extra_symbols(self, function_list):
        '''
        Returns a list of extra symbols required if functions are used inside
        a kernel.

        :param function_list: A list of function names used in a kernel.
        :type function_list: list of str.

        :returns: A list of string, symbol tuples containing the name and \
                symbols required to use a specified function.
        :rtype: list of Tuple( str, Symbol)
        '''
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
        The call to use in kernels to obtain a the magnetic and electric
        field values at the position of a particle.
        After this call, the field values will be stored in ``ex_part``, 
        ``ey_part``, ``ez_part``, ``bx_part``, ``by_part``, and ``bz_part``.

        The inputs to the functions should be the full particle access
        to use each part of the particle, e.g. for part_weight it may be
        ``part1.part_weight``, and dx and dt should be the variables containing
        the cell width and timestep.

        :param part_weight: Access to the particle's weight value.
        :param part_charge: Access to the particle's charge value.
        :param part_mass: Access to the particle's mass value.
        :param part_momentum_x: Access to the particle's x momentum.
        :param part_momentum_y: Access to the particle's y momentum.
        :param part_momentum_z: Access to the particle's z momentum.
        :param dx: Access to the cell width value.
        :param dt: Access to the timestep value.
        '''
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
        '''
        The call to use in kernels to obtain a the magnetic and electric
        field values at the position of a particle.

        :param delta_x: The movement of the particle in the x direction during this step.
        :param part_vy: The part_vy variable.
        :param part_vz: The part_vz variable.
        '''
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
        code = code + f"{in_str}{in_str}" + f"cell_x1, gx, hx, field.scatter_jx, field.scatter_jy, field.scatter_jz, idt, {part_vy}, {part_vz}, idx, dtco2, idtf, idxf,\n"
        code = code + f"{in_str}{in_str}" + "field.nx, fcx, fcy, field.ng);\n"
        return code

    def has_preferred_decomposition(self):
        return True

    def get_preferred_decomposition(self, box_str, current_indent=0, indent=0) -> str:
        '''
        Returns the call to the coupled system which places the preferred
        domain decomposition into the system's boundary object.
        '''
        return " " * current_indent + f"store_domain_decomposition(field, {box_str});\n"

    def copy_files(self):
        '''
        Copies the files required to use this coupled system into the output
        directory.
        '''
        if HartreeParticleDSL.get_mpi():
            subfolder = "MPI_version"
        else:
            subfolder = "serial_version"
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
        BASEPATH = os.path.dirname(inspect.getfile(FDTD_MPI_Kokkos))
        for f in files:
            copy(os.path.join(BASEPATH, subfolder, f), f)

    def compilation_files(self):
        '''
        :returns: A list of strings containing the files required to compile \
                with this coupled system.
        :rtype: list of str.
        '''
        files = ['FDTD_MPI_IO_HDF5.cpp',
                 'FDTD_MPI_boundaries.cpp',
                 'FDTD_MPI_init.cpp',
                 'FDTD_MPI_step.cpp']
        return files

    def get_required_packages(self):
        '''
        :returns: A list of strings containing packages required to compile \
                with CMake.
        :rtype: list of str.
        '''
        return ["Kokkos"]
