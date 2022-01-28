#include "random_io_cabana.hpp"

template void <class aosoa_class>( aosoa_class &particle_system, config_type& config){

    config.config_host(0).boundary.x_min = 0.0;
    config.config_host(0).boundary.x_max = 1.0;
    config.config_host(0).boundary.y_min = 0.0;
    config.config_host(0).boundary.y_max = 1.0;
    config.config_host(0).boundary.z_min = 0.0;
    config.config_host(0).boundary.z_max = 1.0;

    auto core_part_slice = Cabana::slice<core_part_space>(particle_system);
    for(int i = 0; i < particle_system.size(); i++){
        core_part_slice(i).position[0] = (double)(rand()) / (double)(RAND_MAX);
        core_part_slice(i).position[1] = (double)(rand()) / (double)(RAND_MAX);
        core_part_slice(i).position[2] = (double)(rand()) / (double)(RAND_MAX);
        core_part_slice(i).velocity[0] = 0.0;
        core_part_slice(i).velocity[1] = 0.0;
        core_part_slice(i).velocity[2] = 0.0;
    }

    Kokkos::deep_copy(config.config, config.config_host);
}
