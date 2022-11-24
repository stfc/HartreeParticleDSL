#ifndef HPDSL_KOKKOS_FDTD_MPI_INTERPOLATION_H
#define HPDSL_KOKKOS_FDTD_MPI_INTERPOLATION_H

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_ScatterView.hpp>

#include "FDTD_MPI_field.hpp"
#include "FDTD_MPI_step.hpp"

KOKKOS_INLINE_FUNCTION
void interpolate_from_grid_tophat_1D(double part_weight, double part_q, double part_m,
                                     double *part_x, double part_px, double part_py, double part_pz,
                                     double *ex_part, double *ey_part, double *ez_part,
                                     double *bx_part, double *by_part, double *bz_part,
                                     double cell_x_r, int *cell_x1, double *gx, double *hx,
                                     field_type ex, field_type ey, field_type ez,
                                     field_type bx, field_type by, field_type bz, double idt,
                                     double idx, double dtco2, double idtf, double idxf,
                                     int nx, double fcx, double fcy, int ng){

    double part_mc = c * part_m;
    double ipart_mc = 1.0 / part_mc;

    double part_ux = part_px * ipart_mc;
    double part_uy = part_py * ipart_mc;
    double part_uz = part_pz * ipart_mc;

    *cell_x1 = floor(cell_x_r + 0.5);
    double cell_frac_x = ((double)(*cell_x1)) - cell_x_r;
    *cell_x1 = (*cell_x1) + 1;

    //First tophat stuff, gx.inc
    gx[0] = 0.5 + cell_frac_x;
    gx[1] = 0.5 - cell_frac_x;

    int cell_x2 = floor(cell_x_r);
    cell_frac_x = ((double)(cell_x2)) - cell_x_r + 0.5;
    cell_x2 = cell_x2+1;

    int dcellx = 0;
    hx[dcellx] = 0.5 + cell_frac_x;
    hx[dcellx+1] = 0.5 - cell_frac_x;

    //Shift all index accesses up by ng compared to the FDPS version
    *ex_part = hx[0] * ex(cell_x2 + ng) + hx[1] * ex(cell_x2 + 1 + ng);
    *ey_part = gx[0] * ey((*cell_x1) + ng) + gx[1] * ey( (*cell_x1) + 1 + ng);
    *ez_part = gx[0] * ez((*cell_x1) + ng) + gx[1] * ez( (*cell_x1) + 1 + ng);

    *bx_part = gx[0] * bx((*cell_x1) + ng) + gx[1] * bx((*cell_x1) + 1 + ng);
    *by_part = hx[0] * by(cell_x2 + ng) + hx[1] * by(cell_x2 + 1 + ng);
    *bz_part = hx[0] * bz(cell_x2 + ng) + hx[1] * bz(cell_x2 + 1 + ng);

}

KOKKOS_INLINE_FUNCTION
void GatherForcesToGrid_1D(double part_weight, double part_q,
                           double part_x, double delta_x,
                           int cell_x1, double *gx, double *hx,
                           scatter_field_type scatter_jx, scatter_field_type scatter_jy, scatter_field_type scatter_jz,
                           double idt, double part_vy, double part_vz,
                           double idx, double dtco2, double idtf, double idxf,
                           int nx, double fcx, double fcy, int jng){

    //Move particle to t + 1.5dt;
    part_x = part_x + delta_x;
    double cell_x_r = part_x * idx - 0.5;
    int cell_x3 = floor(cell_x_r + 0.5);
    double cell_frac_x = ((double)(cell_x3)) - cell_x_r;
    cell_x3 = cell_x3 + 1;

    for(int i = -1; i < 3; i++){
        hx[i] = 0.0;
    }

    int dcellx = cell_x3 - cell_x1;//Centered on 0 like epoch
    hx[dcellx] = 0.5 + cell_frac_x;
    hx[dcellx+1] = 0.5 - cell_frac_x;

    for(int i = -1; i < 3; i++){
        hx[i] = hx[i] - gx[i];
    }

    int xmin = 0 + (dcellx -1)/ 2;
    int xmax = 1 + (dcellx + 1)/ 2;

    double fjx = fcx * part_q;
    double fjy = fcy * part_q * part_vy;
    double fjz = fcy * part_q * part_vz;

    double jxh = 0.0;

    for(int ix = xmin; ix <= xmax; ix++){
        int cx = cell_x1 + ix;

        double wx = hx[ix];
        double wy = gx[ix] + 0.5 * hx[ix];

        //This is the bit that actually solves d(rho)/dt = -div(J)
        jxh = jxh - fjx * wx;

        double jyh = fjy * wy;
        double jzh = fjz * wy;
        //Scatterview is used here
        auto jx = scatter_jx.access();
        auto jy = scatter_jy.access();
        auto jz = scatter_jz.access();
        jx(cx+jng) += jxh;
        jy(cx+jng) += jyh;
        jz(cx+jng) += jzh;
    }
}

#endif
