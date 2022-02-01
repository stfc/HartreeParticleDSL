#ifndef RANDOM_IO_HPP
#define RANDOM_IO_HPP
#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "random_io_cabana.hpp"

template <class aosoa_class> void random_io( aosoa_class particle_system, config_type& config){

    config.config_host(0).space.box_dims.x_min = 0.0;
    config.config_host(0).space.box_dims.x_max = 1.0;
    config.config_host(0).space.box_dims.y_min = 0.0;
    config.config_host(0).space.box_dims.y_max = 1.0;
    config.config_host(0).space.box_dims.z_min = 0.0;
    config.config_host(0).space.box_dims.z_max = 1.0;

    auto core_part_position_slice = Cabana::slice<core_part_position>(particle_system);
    auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_system);
    for(int i = 0; i < particle_system.size(); i++){
        core_part_position_slice(i,0) = (double)(rand()) / (double)(RAND_MAX);
        core_part_position_slice(i,1) = (double)(rand()) / (double)(RAND_MAX);
        core_part_position_slice(i,2) = (double)(rand()) / (double)(RAND_MAX);
        core_part_velocity_slice(i,0) = 0.0;
        core_part_velocity_slice(i,1) = 0.0;
        core_part_velocity_slice(i,2) = 0.0;
    }

    Kokkos::deep_copy(config.config, config.config_host);
}
#endif
