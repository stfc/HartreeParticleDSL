import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import config
"""
When we import the config object here, we have the default config object.
In C-style syntax, this looks like:
    struct boundary_type{
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        double z_min;
        double z_max;
        double local_x_min;
        double local_x_max;
        double local_y_min;
        double local_y_max;
        double local_z_min;
        double local_z_max;
    };

    struct space_type{
        struct boundary_type box_dims;
    };

    struct config{
        struct space_type space;
        int64_t nparts;
    }

    In general, the members of the space and boundary structures are visible
    to enable boundary condition support (or other spatially relevant computation),
    with the local_ prefixed members only
    being relevant to MPI computation, and in general these should not be used outside
    of backend code.
"""

# Add a new integer to the config to represent a step count.
# We can access the step count using config.step.
part.add_element("step", "int")

# We can also add array properties
part.add_element("energies", "double[3]")

# Finally we tell the DSL to use this config object for this code.
HartreeParticleDSL.set_config_type(config)
