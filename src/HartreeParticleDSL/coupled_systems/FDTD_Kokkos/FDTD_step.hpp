#ifndef HPDSL_KOKKOS_FDTD_STEP_H
#define HPDSL_KOKKOS_FDTD_STEP_H

#include "FDTD_field.hpp"

const double c = 299792458.00000000;
const double epsilon0 = 8.854187817620389850536563031710750e-12;

struct update_e_field_functor{

    //TODO
    FDTD_field _field;
    int _nx;

    KOKKOS_INLINE_FUNCTION
    update_e_field_functor(FDTD_field field, int nx) : _field(field),
                           _nx(nx) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        double cpml_x;
        double c1, c2, c3;
        double cx1, cx2, cx3;

        //Assuming non cpml
        if(_field.field(0).field_order == 2){
           cx1 = _field.field(0).cnx;
           _field.ex(ix) = _field.ex(ix) - _field.field(0).fac * _field.jx(ix);
           int m1 = ix-1;
           _field.ey(ix) = _field.ey(ix) - cx1 * (_field.bz(ix) - _field.bz(m1)) - _field.field(0).fac * _field.jy(ix);
           _field.ez(ix) = _field.ez(ix) + cx1 * (_field.by(ix) - _field.by(m1)) - _field.field(0).fac * _field.jz(ix);
        }else if(_field.field(0).field_order == 4){
            c1 = 9.0/8.0;
            c2 = -1.0 / 24.0;
            cx1 = c1 * _field.field(0).cnx;
            cx2 = c2 * _field.field(0).cnx;
            _field.ex(ix) = _field.ex(ix) - _field.field(0).fac * _field.jx(ix);
            int m1 = ix-1;
            int m2 = ix-2;
            int p1 = ix+1;
            _field.ey(ix) = _field.ey(ix) - cx1 * (_field.bz(ix) - _field.bz(m1)) - cx2 *
                             (_field.bz(p1) - _field.bz(m2)) - _field.field(0).fac * _field.jy(ix);
            _field.ez(ix) = _field.ez(ix) + cx1 * (_field.by(ix) - _field.by(m1)) + cx2 *
                             (_field.by(p1) - _field.by(m2)) - _field.field(0).fac * _field.jz(ix);
        }else{
            c1 = 75.0/64.0;
            c2 = -25.0 / 384.0;
            c3 = 3.0 / 640.0;
            cx1 = c1 * _field.field(0).cnx;
            cx2 = c2 * _field.field(0).cnx;
            cx3 = c3 * _field.field(0).cnx;

            int m1 = ix-1;
            int m2 = ix-2;
            int m3 = ix-3;
            int p1 = ix+1;
            int p2 = ix+2;
            _field.ex(ix) = _field.ex(ix) - _field.field(0).fac * _field.jx(ix);
            _field.ey(ix) = _field.ey(ix) - cx1 * (_field.bz(ix) - _field.bz(m1))
                   - cx2 * (_field.bz(p1) - _field.bz(m2)) - cx3 * (_field.bz(p2) - _field.bz(m3)) -
                   _field.field(0).fac * _field.jy(ix);
            _field.ez(ix) = _field.ez(ix) + cx1 * (_field.by(ix) - _field.by(m1))
                   + cx1 * (_field.by(p1) - _field.by(m2)) + cx3 * (_field.by(p2) - _field.by(m3)) -
                   _field.field(0).fac * _field.jy(ix);
        }
    }
};

struct update_b_field_functor{
    //TODO
    FDTD_field _field;
    int _nx;

    KOKKOS_INLINE_FUNCTION
    update_b_field_functor(FDTD_field field, int nx) : _field(field), _nx(nx) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix) const{
        double cpml_x;
        double c1, c2, c3;
        double cx1, cx2, cx3;
        //Assuming non cpml and using maxwell_solver_yee
        if(_field.field(0).field_order == 2){
            //Yee solver
            cx1 = _field.field(0).hdtx;
            int p1 = ix+1;
            _field.by(ix) = _field.by(ix) + cx1 * (_field.ez(p1) - _field.ez(ix));
            _field.bz(ix) = _field.bz(ix) - cx1 * (_field.ey(p1) - _field.ey(ix));
        }else if(_field.field(0).field_order == 4){
            c1 = 9.0/8.0;
            c2 = -1.0 / 24.0;

            cx1 = c1 * _field.field(0).hdtx;
            cx2 = c2 * _field.field(0).hdtx;

            int m1 = ix-1;
            int p1 = ix+1;
            int p2 = ix+2;
            _field.by(ix) = _field.by(ix) + cx1 * (_field.ez(p1) - _field.ez(ix)) + cx2 * (_field.ez(p2) - _field.ez(m1));
            _field.bz(ix) = _field.bz(ix) - cx1 * (_field.ey(p1) - _field.ey(ix)) - cx2 * (_field.ey(p2) - _field.ey(m1));

        }else{
            c1 = 75.0/64.0;
            c2 = -25.0 / 384.0;
            c3 = 3.0 / 640.0;

            cx1 = c1 * _field.field(0).hdtx;
            cx2 = c2 * _field.field(0).hdtx;
            cx3 = c3 * _field.field(0).hdtx;

            int m1 = ix-1;
            int m2 = ix-2;
            int p1 = ix+1;
            int p2 = ix+2;
            int p3 = ix+3;
            _field.by(ix) = _field.by(ix) + cx1 * (_field.ez(p1) - _field.ez(ix)) + cx2 * (_field.ez(p2) - _field.ez(m1)) +
                              cx3 * (_field.ez(p3) - _field.ez(m2));
            _field.bz(ix) = _field.bz(ix) - cx1 * (_field.ey(p1) - _field.ey(ix)) - cx2 * (_field.ey(p2) - _field.ey(m1)) -
                              cx3 * (_field.ey(p3) - _field.ey(m2));
        }
    }
};

void update_eb_fields_half_1D(FDTD_field &field, const int nx, const int ng, const double dt, const double dx,
                              update_e_field_functor &update_e_field,
                              update_b_field_functor &update_b_field,
                              Kokkos::RangePolicy<> rangepolicy);

void update_eb_fields_final_1D(FDTD_field &field, const int nx, const int ng, const double dt, const double dx,
                              update_e_field_functor &update_e_field,
                              update_b_field_functor &update_b_field,
                              Kokkos::RangePolicy<> rangepolicy);

void current_start(FDTD_field &field, const int nx, const int jng); 
#endif
