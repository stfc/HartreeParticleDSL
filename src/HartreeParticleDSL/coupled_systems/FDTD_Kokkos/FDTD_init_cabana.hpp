#ifndef HPDSL_KOKKOS_FDTD_INIT_CABANA_H
#define HPDSL_KOKKOS_FDTD_INIT_CABANA_H

#include "part.h"
#include <random>

std::default_random_engine generator;
double momentum_from_temperature(double mass, double temperature, double drift){                                                                                             double stdev, mu;
  const double kb = 1.3806488e-23;  // J/K
  stdev = sqrt(temperature * mass * kb);
  mu = drift;
  //random_box_muller
  std::normal_distribution<double> d{mu, stdev};
  double x = d(generator);
  return x;
}

template <class aosoa_class> void fdtd_init_particles( aosoa_class particle_aosoa, config_type& config, FDTD_field &field){

    int i0 = 1; //-ng, number ghost cells, which we're setting to 0
    int i1 = 0; // 1 - i0  

    int *npart_per_cell_array = (int*)malloc(sizeof(int) * 2); //i0:nx+i1 in fortran
    for (int i = 0; i < 2; i++){
        npart_per_cell_array[i] = 0;
    }
    
    int nxglobal = field.nx;
    int ncell = nxglobal;
    config.config_host(0).space.box_dims.x_max = 2.0;
    config.config_host(0).space.box_dims.x_min = 0.0;
    double cell_size = ( config.config_host(0).space.box_dims.x_max - config.config_host(0).space.box_dims.x_min ) / (double)(nxglobal);
    double dx = 2.0 / (double)(nxglobal);

    //For each specie
    //Grab species initial condition
    double species_charge[2];
    species_charge[0] = -1.0;
    species_charge[1] = 1.0;
    double species_mass[2];
    species_mass[0] = 1.0;
    species_mass[1] = 3670.48294;
    double species_density[2];
    species_density[0] = 1e19;
    species_density[1] = (1.0 - 2.0 * 1e-3) * 1e19;

    double *density = (double*)malloc(sizeof(double) * nxglobal);
    bool *density_map = (bool*)malloc(sizeof(bool) * nxglobal);

    double species_count[2];
    species_count[0] = 1425000;
    species_count[1] = 1425000;
    double counter = 0;
    auto id_slice = Cabana::slice<id>(particle_aosoa);
    auto mass_slice = Cabana::slice<mass>(particle_aosoa);
    auto charge_slice = Cabana::slice<charge>(particle_aosoa);
    auto weight_slice = Cabana::slice<weight>(particle_aosoa);
    auto part_pos_slice = Cabana::slice<core_part_position>(particle_aosoa);
    auto part_p_slice = Cabana::slice<part_p>(particle_aosoa);
    for(int i = 0; i < 2; i++){
        double npart_per_cell_average = (double)(species_count[i]) / (double)(ncell);
        if(npart_per_cell_average <= 0) continue;
        for(int ix = 0; ix < ncell; ix++){
            density[ix] = species_density[i];
            density_map[ix] = true;
        }
        int ipart = 0;
        int start = counter;
        for(int ix = 0;  ix < ncell; ix++){
            int npart_per_cell = round(npart_per_cell_average); //Assuming density for every cell is the same for now
            for(int ipart = 0; ipart < npart_per_cell; ipart++){
                if (counter > particle_aosoa.size()) break;
                id_slice(counter) = counter;
                charge_slice(counter) = species_charge[i] * 1.602176565e-19;
                mass_slice(counter) = species_mass[i] * 9.10938291e-31;
                double x = (double)(ix+1) * cell_size;
                double r = (double) ((double)rand() / ((double)RAND_MAX));
                part_pos_slice(counter, 0) = x + (r - 0.5) * cell_size;
                counter++;

            }
        }
        //Load particles finished
        int *npart_in_cell = (int*)malloc(sizeof(int) * nxglobal);
        for(int ii = 0; ii < nxglobal; ii++){
            npart_in_cell[ii] = 0;
        }
        //Loop over particles set
        for(int ipart = start; ipart < counter; ipart++){
            double cell_x_r = part_pos_slice(ipart, 0) / dx - 0.5;
            int cell_x = floor(cell_x_r + 0.5);
            double cell_frac_x = (double)(cell_x) - cell_x_r;
            cell_x = cell_x + 1;
            double gx[2];
            gx[0] = 0.5 + cell_frac_x;
            gx[1] = 0.5 - cell_frac_x;
            //      Calculate density at the particle position
            double wdata = 0.0;
            for(int isubx = 0; isubx < 2; isubx++){
                int ii = cell_x + isubx;
                if (ii >= nxglobal) ii -= nxglobal;
                if (!density_map[ii]) ii = cell_x + 1 - isubx;
                wdata = wdata + gx[isubx] * density[ii];
            }
            weight_slice(ipart) = wdata;
            if( gx[1] > gx[0]) cell_x = cell_x + 1;
            if (cell_x >= nxglobal) cell_x -= nxglobal;
            npart_in_cell[cell_x] = npart_in_cell[cell_x] + 1;
        }
        //Loop again to normalise the weights
        double wdata = dx;
        for(int ipart = start; ipart < counter; ipart++){
            int cell_x = floor(part_pos_slice( ipart, 0) / dx + 1.5);
            if(cell_x >= nxglobal) cell_x -= nxglobal;
            double lweight = weight_slice(ipart);
            weight_slice(ipart) = lweight * wdata / npart_in_cell[cell_x];
        }
        free(npart_in_cell);
        std::cout << counter << " particles loaded\n";

        //FIX BCS here
        double **species_temp = (double**) malloc(sizeof(double*) * 3);
        species_temp[0] = (double*) malloc(sizeof(double) * nxglobal);
        species_temp[1] = (double*) malloc(sizeof(double) * nxglobal);
        species_temp[2] = (double*) malloc(sizeof(double) * nxglobal);
        for(int n = 0; n < 3; n++){
            for(int ix = 0; ix < nxglobal; ix++){                                                                                                                                          //We use a specific temperature everywhere for this example
                species_temp[n][ix] = 11594200.000;
            }
        }
        //setup_ic_drift
        //Drift is all 0 for now.
        for(int n = 0; n < 3; n++){
            for(int ipart = start; ipart < counter; ipart++){
                double lmass = mass_slice(ipart);
                //include particle_to_grid.inc
                double cell_x_r = (part_pos_slice(ipart, 0)) / dx - 0.5;
                int cell_x = floor(cell_x_r + 0.5);
                double cell_frac_x = (double)(cell_x) - cell_x_r;
                double gx[2];
                gx[0] = 0.5 + cell_frac_x;
                gx[1] = 0.5 - cell_frac_x;
                double temp_local = 0.0;
                double drift_local = 0.0;
                for(int ix = 0; ix < 2; ix++){
                    int index = cell_x + ix;
                    if (index >= nxglobal) index -= nxglobal; //All values are the same so doesn't matter
                    temp_local += gx[ix] * species_temp[n][index];
                    drift_local += gx[ix] * 0.0; // Drift is 0
                }
                //Calculate momentum, directions are 0/1/2 to x/y/z
                part_p_slice(ipart, n) = momentum_from_temperature(lmass, temp_local, drift_local);
            }
        }
        free(species_temp[0]);
        free(species_temp[1]);
        free(species_temp[2]);
        free(species_temp);
    }
    free(npart_per_cell_array);
    free(density);
    free(density_map);
}

#endif
