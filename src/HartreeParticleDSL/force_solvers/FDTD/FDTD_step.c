#include "FDTD_step.h"

const PS::F64 c = 299792458.00000000;
const PS::F64 epsilon0 = 8.854187817620389850536563031710750e-12;

void update_e_field_1D(struct FDTD_field *field){

    double cpml_x;
    double c1, c2, c3;
    double cx1, cx2, cx3;

    int nx = field->nx;

    //Only non cpml for now
    if(field->field_order == 2){
        cx1 = field->cnx;
        for(int ix = 0; ix <= nx; ix++){
            field->ex[ix] = field->ex[ix] - field->fac * field->jx[ix];
            int m1 = ix-1;
            field->ey[ix] = field->ey[ix] - cx1 * (field->bz[ix] - field->bz[m1]) - field->fac * field->jy[ix];
            field->ez[ix] = field->ez[ix] + cx1 * (field->by[ix] - field->by[m1]) - field->fac * field->jz[ix];
        }
    }else if(field->field_order == 4){
        c1 = 9.0/8.0;
        c2 = -1.0 / 24.0;
        cx1 = c1 * field->cnx;
        cx2 = c2 * field->cnx;
        for(int ix = 0; ix < nx; ix++){
            field->ex[ix] = field->ex[ix] - field->fac * field->jx[ix];
            int m1 = ix-1;
            int m2 = ix-2;
            int p1 = ix+1;
            field->ey[ix] = field->ey[ix] - cx1 * (field->bz[ix] - field->bz[m1]) - cx2 * (field->bz[p1] - field->bz[m2]) - field->fac * field->jy[ix];
            field->ez[ix] = field->ez[ix] + cx1 * (field->by[ix] - field->by[m1]) + cx2 * (field->by[p1] - field->by[m2]) - field->fac * field->jz[ix];
        }
    }else{ //field_order is 6
        c1 = 75.0/64.0;
        c2 = -25.0 / 384.0;
        c3 = 3.0 / 640.0;
        cx1 = c1 * field->cnx;
        cx2 = c2 * field->cnx;
        cx3 = c3 * field->cnx;

        for(int ix = 0; ix < nx; ix++){
            int m1 = ix-1;
            int m2 = ix-2;
            int m3 = ix-3;
            int p1 = ix+1;
            int p2 = ix+2;
            field->ex[ix] = field->ex[ix] - field->fac * field->jx[ix];
            field->ey[ix] = field->ey[ix] - cx1 * (field->bz[ix] - field->bz[m1])
                   - cx2 * (field->bz[p1] - field->bz[m2]) - cx3 * (field->bz[p2] - field->bz[m3]) - field->fac * field->jy[ix];
            field->ez[ix] = field->ez[ix] + cx1 * (field->by[ix] - field->by[m1])
                   + cx1 * (field->by[p1] - field->by[m2]) + cx3 * (field->by[p2] - field->by[m3]) - field->fac * field->jy[ix];
        }
    }

}

void field_bcs(double *field, int, nx, int ng){

   double *temp = (double*) malloc(sizeof(double)*ng);

   //Update top
   for(int i = 0; i < ng; i++){
       temp[i] = field[i+1];
   }
   for(int i = 0; i < ng; i++){
       field[nx+i+1] = temp[i];
   }

   //Update bottom
   for(int i = 0; i < ng; i++){
       temp[i] = field[nx+i+1-ng];
   }
   for(int i = 0; i < ng; i++){
       field[1+i-ng] = temp[i];
   }

   free(temp);

}

void efield_bcs_1D(struct FDTD_field *field){

    field_bcs(field->ex, field->nx, field->ng);
    field_bcs(field->ey, field->nx, field->ng);
    field_bcs(field->ez, field->nx, field->ng);
}


void bfield_bcs_1D(struct FDTD_field *field, bool mpi_only){

    field_bcs(field->bx, field->nx, field->ng);
    field_bcs(field->by, field->nx, field->ng);
    field_bcs(field->bz, field->nx, field->ng);

}

void bfield_final_bcs_1D(struct FDTD_field *field){
    bfield_bcs_1D(field, true);
    bfield_bcs_1D(field, false);
}

void update_b_field_1D(struct FTDT_field *field){


    double cpml_x;
    double c1, c2, c3;
    double cx1, cx2, cx3;
    int nx = field->nx;
    //Assuming non cpml and using maxwell_solver_yee i guess
    if(field->field_order == 2){
        //Yee solver
        cx1 = field->hdtx;
        for(int ix = 0; ix <= nx; ix++){
            int p1 = ix+1;
            field->by[ix] = field->by[ix] + cx1 * (field->ez[p1] - field->ez[ix]);
            field->bz[ix] = field->bz[ix] - cx1 * (field->ey[p1] - field->ey[ix]);
        }
    }else if(field->field_order == 4){
        c1 = 9.0/8.0;
        c2 = -1.0 / 24.0;

        cx1 = c1 * field->hdtx;
        cx2 = c2 * field->hdtx;

        for(int ix = 0; ix < nx; ix++){
            int m1 = ix-1;
            int p1 = ix+1;
            int p2 = ix+2;
            field->by[ix] = field->by[ix] + cx1 * (field->ez[p1] - field->ez[ix]) + cx2 * (field->ez[p2] - field->ez[m1]);
            field->bz[ix] = field->bz[ix] - cx1 * (field->ey[p1] - field->ey[ix]) - cx2 * (field->ey[p2] - field->ey[m1]);
        }

    }else{ //field_order is 6
        c1 = 75.0/64.0;
        c2 = -25.0 / 384.0;
        c3 = 3.0 / 640.0;

        cx1 = c1 * field->hdtx;
        cx2 = c2 * field->hdtx;
        cx3 = c3 * field->hdtx;

        for(int ix = 0; ix < nx; ix++){
            int m1 = ix-1;
            int m2 = ix-2;
            int p1 = ix+1;
            int p2 = ix+2;
            int p3 = ix+3;
            field->by[ix] = field->by[ix] + cx1 * (field->ez[p1] - field->ez[ix]) 
                + cx2 * (field->ez[p2] - field->ez[m1]) + cx3 * (field->ez[p3] - field->ez[m2]);
            field->bz[ix] = field->bz[ix] - cx1 * (field->ey[p1] - field->ey[ix]) 
                - cx2 * (field->ey[p2] - field->ey[m1]) - cx3 * (field->ey[p3] - field->ey[m2]);
        }
    }

}


void update_eb_fields_half_1D(struct FDTD_field *field, const int nx, const double dt, const double dx){

    field->hdt = dt / 2.0;
    field->hdtx = field->hdt / dx;
    field->cnx = field->hdtx * c*c;
    field->fac = field->hdt / epsilon0;

    update_e_field_1D(field);
    efield_bcs_1D(field);
    update_b_field_1D(field);
    bfield_bcs_1D(field);

}

void update_eb_fields_final_1D(struct FDTD_field *field, const int nx, const double dt, const double dx){

    field->hdt = dt / 2.0;
    field->hdtx = field->hdt / dx;
    field->cnx = field->hdtx * c*c;
    field->fac = field->hdt / epsilon0;

    update_b_field_1D(field);
    bfield_final_bcs_1D(field);
    update_e_field_1D(field);
    efield_bcs_1D(field);
}
