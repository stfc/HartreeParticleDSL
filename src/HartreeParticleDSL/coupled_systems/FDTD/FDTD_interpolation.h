#ifndef HPDSL_FDTD_H
#define HPDSL_FDTD_H

#include "FDTD_field.h"

void interpolate_from_grid_tophat_1D(double part_weight, double part_q, double part_m,
                                     double *part_x, double part_p_x, double part_p_y,
                                     double part_p_z, double *ex_part,
                                     double *ey_part, double *ez_part, double *bx_part,
                                     double *by_part, double *bz_part, double cell_x_r,
                                     int32_t *cell_x1, struct FDTD_field *field,
                                     double idt, double idx, double dtco2, double idtf,
                                     double idxf, int32_t nx, double fcx, double fcy);

void gather_forces_to_grid_tophat_1D(double part_weight, double part_q,
                                     double part_x, double delta_x,
                                     int32_t cell_x1, struct FDTD_field *field,
                                     double idt, double part_vy, double part_vz,
                                     double idx, double dtco2, double idtf, double idxf,
                                     int32_t nx, double fcx, double fcy);

#endif
