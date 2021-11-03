#include "FDTD_init.h"

void fdtd_initialize_1D(struct *FDTD_field field, int nx, int ng){

    double *ex_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *ey_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *ez_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *bx_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *by_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *bz_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *jx_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *jy_base = (double*) malloc(sizeof(double) * (nx + 2*ng));
    double *jz_base = (double*) malloc(sizeof(double) * (nx + 2*ng));

    double *hx = (double*) malloc(sizeof(double) * 4);
    double *gx = (double*) malloc(sizeof(double) * 4);

    for(int i = 0; i < nx+2*ng; i++){
        ex_base[i] = 0.0;
        ey_base[i] = 0.0;
        ez_base[i] = 0.0;
        bx_base[i] = 0.0;
        by_base[i] = 0.0;
        bz_base[i] = 0.0;
        jx_base[i] = 0.0;
        jy_base[i] = 0.0;
        jz_base[i] = 0.0;
    }

    field->nx = nx;
    field->ng = ng;

    //Gonna do non-0 bounds for now
    field->ex = &ex_base[ng-1];
    field->ey = &ey_base[ng-1];
    field->ez = &ez_base[ng-1];
    field->bx = &bx_base[ng-1];
    field->by = &by_base[ng-1];
    field->bz = &bz_base[ng-1];
    field->jx = &jx_base[ng-1];
    field->jy = &jy_base[ng-1];
    field->jz = &jz_base[ng-1];

    field->hx = &hx[1]; 
    field->gx = &gx[1]; 
}

void FDTD_cleanup_1D(struct *FDTD_field field){

    int ng = field->ng;
    free(&(field->ex[-ng+1]));
    free(&(field->ey[-ng+1]));
    free(&(field->ez[-ng+1]));
    free(&(field->bx[-ng+1]));
    free(&(field->by[-ng+1]));
    free(&(field->bz[-ng+1]));
    free(&(field->jx[-ng+1]));
    free(&(field->jy[-ng+1]));
    free(&(field->jz[-ng+1]));
    
    free(&(hx[-1]));
    free(&(gx[-1]));
}


void fdtd_initialize_example_1D(struct *FDTD_field field, int nx, int ng){

    //Everything is 0 but bz
    for(int i = -ng+1; i < ng+nx; i++){
        field->bz[i] = 2.1;
    }

    field->field_order = 2;
    field->fng = (double)(field->field_order) / 2.0;
    if (field->field_order == 2){
       field->cfl = 1.0;
    }else if(field->field_order == 4){
        field->cfl = 6.0/7.0;
    }else{
        field->cfl = 120.0/149.0;
    }
}
