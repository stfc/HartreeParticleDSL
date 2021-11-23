#ifndef HPDSL_FDTD_INIT_H
#define HPDSL_FDTD_INIT_H

#include "FDTD_field.h"

void fdtd_initialize_1D(struct *FDTD_field field, int nx, int ng);

void FDTD_cleanup_1D(struct *FDTD_field field);

void fdtd_initialize_example_1D(struct *FDTD_field field, int nx, int ng);

#endif
