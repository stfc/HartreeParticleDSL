#include <random_io.hpp>
#include <stdlib>

void random_io(PS::ParticleSystem<FullParticle>& particle_system, config_type& config){
    //Set up a 1x1x1 cube for now
    PS::DomainInfo dinfo;
    config.box.x_min = 0.0;
    config.box.x_max = 1.0;
    config.box.y_min = 0.0;
    config.box.y_max = 1.0;
    config.box.z_min = 0.0;

    dinfo.initialize();
    dinfo.setBoundaryCondition ( PS::BOUNDARY_CONDITION_PERIODIC_XYZ );
    dinfo.setPosRootDomain(PS::F64vec(config.box.x_min,
                                      config.box.y_min,
                                      config.box.z_min),
                           PS::F64vec(config.box.x_max,
                                      config.box.y_max,
                                      config.box.z_max));

    for(PS::S32 i = 0 ; i < pic_system.getNumberOfParticleLocal() ; i++){
        particle_system[p].core_part.position[0] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[p].core_part.position[1] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[p].core_part.position[2] = (double)(rand()) / (double)(RAND_MAX);
        particle_system[p].core_part.velocity[0] = 0.0;
        particle_system[p].core_part.velocity[1] = 0.0;
        particle_system[p].core_part.velocity[2] = 0.0;
    }
}
