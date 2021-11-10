#ifndef HPDSL_FDTD_FIELD_H
#define HPDSL_FDTD_FIELD_H

struct FDTD_field{
    double *ex;
    double *ey;
    double *ez;
    double *bx;
    double *by;
    double *bz;
    double *jx;
    double *jy;
    double *jz;
    double *hx;
    double *gx;

    int nx;
    int ng;

    double hdt;
    double hdtx;
    double cnx;
    double fac;
    int field_order;
    double fng;
    double cfl;
    double x_grid_min_local;
    double x_grid_max_local;
};

#endif
