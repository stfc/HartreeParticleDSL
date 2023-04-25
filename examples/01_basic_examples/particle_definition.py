import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import part
"""
When we import the part object here, we have the default part object.
In C-style syntax, this looks like:
    struct core_part_space{
        double[3] position;
        double[3] velocity;
    };

    struct part{
        struct core_part_space core_part;
    }

We can add new members to the particle that we might want for our computation
"""

# Add a new integer to the particle to represent an "ID"
# We can access the ID from a kernel as part.ID.
part.add_element("ID", "int")

# We can also add array members to represent multidimensional properties
part.add_element("momentum", "double[3]")

# Finally we tell the DSL to use this part object for this code.
HartreeParticleDSL.set_particle_type(part)
