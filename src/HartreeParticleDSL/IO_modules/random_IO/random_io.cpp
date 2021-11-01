#include "random_io.hpp"

void random_io(PS::ParticleSystem<FullParticle>& particle_system, config_type& config){
    //Set up a 1x1x1 cube for now
    PS::DomainInfo dinfo;
    config.space.box_dims.x_min = 0.0;
    config.space.box_dims.x_max = 1.0;
    config.space.box_dims.y_min = 0.0;
    config.space.box_dims.y_max = 1.0;
    config.space.box_dims.z_min = 0.0;

    dinfo.initialize();
    dinfo.setBoundaryCondition ( PS::BOUNDARY_CONDITION_PERIODIC_XYZ );
    dinfo.setPosRootDomain(PS::F64vec(config.space.box_dims.x_min,
                                      config.space.box_dims.y_min,
                                      config.space.box_dims.z_min),
                           PS::F64vec(config.space.box_dims.x_max,
                                      config.space.box_dims.y_max,
                                      config.space.box_dims.z_max));

    for(PS::S32 i = 0 ; i < particle_system.getNumberOfParticleLocal() ; i++){
        particle_system[i].core_part.position[0] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[i].core_part.position[1] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[i].core_part.position[2] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[i].core_part.velocity[0] = 0.0;
        particle_system[i].core_part.velocity[1] = 0.0;
        particle_system[i].core_part.velocity[2] = 0.0;
    }
}
