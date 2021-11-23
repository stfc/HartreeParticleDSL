#ifndef HPDSL_FDTD_STEP_H
#define HPDSL_FDTD_STEP_H

#include "FDTD_field.h";

void update_eb_fields_half_1D(struct FDTD_field *field, const int nx, const double dt, const double dx);

void update_eb_fields_final_1D(struct FDTD_field *field, const int nx, const double dt, const double dx);

#endif
