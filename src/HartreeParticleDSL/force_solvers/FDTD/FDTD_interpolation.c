#include "FDTD_interpolation.h"

void interpolate_from_grid_tophat_1D(double part_weight, double part_q, double part_m,
                                     double *part_x, double part_p_x, double part_p_y,
                                     double part_p_z, double *ex_part,
                                     double *ey_part, double *ez_part, double *bx_part,
                                     double *by_part, double *bz_part, double cell_x_r,
                                     int32_t *cell_x1, struct FDTD_field *field,
                                     double idt, double idx, double dtco2, double idtf,
                                     double idxf, int32_t nx, double fcx, double fcy){

        //Each particle has its own mass and charge
        double part_mc = c * part_m;
        double ipart_mc = 1.0 / part_mc;

        double part_ux = part_p_x * ipart_mc;
        double part_uy = part_p_y * ipart_mc;
        double part_uz = part_p_z * ipart_mc;

        //Round cell position to nearest cell
        *cell_x1 = floor(cell_x_r + 0.5);
        //Calculate fraction of cell between nearest cell boundary and particle
        double cell_frac_x = ((double)(*cell_x1)) - cell_x_r;
        *cell_x1 = (*cell_x1) + 1;

        field->gx[0] = 0.5 + cell_frac_x;
        field->gx[1] = 0.5 - cell_frac_x;

        int32_t cell_x2 = floor(cell_x_r);
        cell_frac_x = ((double)(cell_x2)) - cell_x_r + 0.5;
        cell_x2 = cell_x2 + 1;

        int32_t dcellx = 0;

        field->hx[dcellx]   = 0.5 + cell_frac_x;
        field->hx[dcellx+1] = 0.5 - cell_frac_x;

        *ex_part = field->hx[0] * field->ex[cell_x2] + field->hx[1] * field->ex[cell_x2 + 1];
        *ey_part = field->gx[0] * field->ey[*cell_x1] + field->gx[1] * field->ey[(*cell_x1) + 1];
        *ez_part = field->gx[0] * field->ez[*cell_x1] + field->gx[1] * field->ez[(*cell_x1) + 1];

        *bx_part = field->gx[0] * field->bx[*cell_x1] + field->gx[1] * field->bx[(*cell_x1) + 1];
        *by_part = field->hx[0] * field->by[cell_x2] + field->hx[1] * field->by[cell_x2 + 1];
        *bz_part = field->hx[0] * field->bz[cell_x2] + field->hx[1] * field->bz[cell_x2 + 1];
}

void gather_forces_to_grid_tophat_1D(double part_weight, double part_q,
                                     double part_x, double delta_x,
                                     int32_t cell_x1, struct FDTD_field *field,
                                     double idt, double part_vy, double part_vz,
                                     double idx, double dtco2, double idtf, double idxf,
                                     int32_t nx, double fcx, double fcy){


    part_x = part_x + delta_x;
    double cell_x_r = part_x * idx - 0.5;
    int32_t cell_x3 = floor(cell_x_r + 0.5);
    double cell_frac_x = ((double)(cell_x3)) - cell_x_r;
    cell_x3 = cell_x3 + 1;

    for(int i = -1; i < 3; i++){
        field->hx[i] = 0.0;
    }

    int32_t dcellx = cell_x3 - cell_x1;

    field->hx[dcellx] = 0.5 + cell_frac_x;
    field->hx[dcellx+1] = 0.5 - cell_frac_x;

    for(int i = -1; i < 3; i++){
        field->hx[i] = field->hx[i] - field->gx[i];
    }

    int32_t xmin = 0 + (dcellx -1)/ 2; //sf_min = 0, sf_max = 1
    int32_t xmax = 1 + (dcellx + 1)/ 2;

    double fjx = fcx * part_q;
    double fjy = fcy * part_q * part_vy;
    double fjz = fcy * part_q * part_vz;

    double jxh = 0.0;
    for( int32_t ix = xmin; ix <= xmax; ix++){
        int32_t cx = cell_x1 + ix;
//        if(cx >= nx) cx -= nx;
//        if(cx < 0) cx += nx;
        double wx = field->hx[ix];
        double wy = field->gx[ix] + 0.5 * field->hx[ix];

        //This is the bit that actually solves d(rho)/dt = -div(J)
        jxh = jxh - fjx * wx;
        double jyh = fjy * wy;
        double jzh = fjz * wy;
        field->jx[cx] = field->jx[cx] + jxh;
        field->jy[cx] = field->jy[cx] + jyh;
        field->jz[cx] = field->jz[cx] + jzh;
    }

}
