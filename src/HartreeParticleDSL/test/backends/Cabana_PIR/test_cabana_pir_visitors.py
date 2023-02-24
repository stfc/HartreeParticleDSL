from __future__ import annotations

from HartreeParticleDSL.HartreeParticleDSL import set_mpi

from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import *
from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import *
from HartreeParticleDSL.backends.AST_to_Particle_IR.ast_to_pir_visitors import *

from HartreeParticleDSL.inbuilt_kernels.boundaries.periodic_boundaries import periodic_boundaries

from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, ScalarType, \
        StructureType, ArrayType, PointerType, \
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, reset_part_and_config

import ast
import inspect
import textwrap
import pytest

def test_pir_cabana_visit_break():
    def a():
        while(True):
            break
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    a = ast_to_pir_visitor()
    pir = a.visit(c)

    b = Cabana_PIR_Visitor(None)
    out = b(pir)
    correct = '''void a(){
    while(true){
        break;
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_return():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        return
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    return;
}
'''
    assert correct == out

    def b():
        return 1
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''int b(){
    return 1;
}
'''
    assert correct == out

    def e():
        return 1.0
    c = ast.parse(textwrap.dedent(inspect.getsource(e)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''double e(){
    return 1.0;
}
'''
    assert correct == out

    def f():
        return False
    c = ast.parse(textwrap.dedent(inspect.getsource(f)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''bool f(){
    return false;
}
'''
    assert correct == out

    def g():
        return "string"
    c = ast.parse(textwrap.dedent(inspect.getsource(g)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''char* g(){
    return "string";
}
'''
    assert correct == out

    def d(a: c_float):
        return a
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''float d(float a){
    return a;
}
'''
    assert correct == out

    def h():
        return a_call()

    c = ast.parse(textwrap.dedent(inspect.getsource(h)))

    pir = astpir.visit(c)
    with pytest.raises(NotImplementedError):
        out = cpir(pir)

    temp = ScalarType(ScalarType.Intrinsic.INTEGER, 128)
    type_mapping_str["temp"] = temp
    def k(a: temp):
        return a
    c = ast.parse(textwrap.dedent(inspect.getsource(k)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''temp k(temp a){
    return a;
}
'''
    assert correct == out
    del(type_mapping_str["temp"])


    temp = PointerType(INT_TYPE)
    type_mapping_str["temp"] = temp
    def m(b: c_int, a: temp):
        return a
    c = ast.parse(textwrap.dedent(inspect.getsource(m)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''temp m(int b, *temp a){
    return *a;
}
'''
    assert correct == out
    del(type_mapping_str["temp"])

    temp = ArrayType(INT_TYPE, [3,3])
    type_mapping_str["temp"] = temp
    def x(b: c_int):
        create_variable(temp, f)
        return b

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''int x(int b){
    int[3][3] f;
    return b;
}
'''
    assert correct == out
    del(type_mapping_str["temp"])

    temp = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC, ArrayType.Extent.DYNAMIC])
    type_mapping_str["temp"] = temp
    def xx(b: c_int):
        create_variable(temp, f)
        return b

    c = ast.parse(textwrap.dedent(inspect.getsource(xx)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''int xx(int b){
    int** f;
    return b;
}
'''
    assert correct == out
    del(type_mapping_str["temp"])


    temp = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC, 3])
    type_mapping_str["temp"] = temp
    def yy(b: c_int):
        create_variable(temp, f)
        return b

    c = ast.parse(textwrap.dedent(inspect.getsource(yy)))
    pir = astpir.visit(c)
    with pytest.raises(UnsupportedCodeError) as excinfo:
        out = cpir(pir)
    assert ("Got an ArraySymbol with a mixed dynamic and fixed typing, "
            "which is not supported." in str(excinfo.value))
    del(type_mapping_str["temp"])

    temp = ArrayType(INT_TYPE, [3, ArrayType.Extent.DYNAMIC])
    type_mapping_str["temp"] = temp
    def zz(b: c_int):
        create_variable(temp, f)
        return b

    c = ast.parse(textwrap.dedent(inspect.getsource(zz)))
    pir = astpir.visit(c)
    with pytest.raises(UnsupportedCodeError) as excinfo:
        out = cpir(pir)
    assert ("Got an ArraySymbol with a mixed dynamic and fixed typing, "
            "which is not supported." in str(excinfo.value))
    del(type_mapping_str["temp"])

def test_pir_cabana_visit_ifdef():
    def a():
        create_variable(c_int, b)
        if True:
            b = 1
        else:
            b = 2

    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)

    correct = '''void a(){
    int b;
    if(true){
        b = 1;
    }else{
        b = 2;
    }
}
'''
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))

    pir = astpir.visit(c)
    out = cpir(pir)
    assert correct == out

    def b():
        create_variable(c_int, b)
        if b < 1:
            b = 1
        elif b < 2:
            b = 2
        elif b < 3:
            b = 3
        else:
            b = 4   
    c = ast.parse(textwrap.dedent(inspect.getsource(b)))

    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void b(){
    int b;
    if((b < 1)){
        b = 1;
    }else if((b < 2)){
        b = 2;
    }else if((b < 3)){
        b = 3;
    }else{
        b = 4;
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_loop():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        create_variable(c_int, b)
        for b in range(0, 3):
            b = b + 1
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    int b;
    for(b = 0; b < 3; b += 1){
        b = (b + 1);
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_while():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a():
        while True:
            break
    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    while(true){
        break;
    }
}
'''
    assert correct == out

def test_pir_cabana_visit_operations():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    def a(b: c_int, c: c_double, d: c_bool):
        b = b + 1
        c = c - 1.0
        b = b * 2
        c = c / 2.0

        d = b < 2
        d = b <= 3
        d = b > 4
        d = b >= 1
        d = b == 2
        d = d and d
        d = d or d

        d = not d
        b = -1


    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(int b, double c, bool d){
    b = (b + 1);
    c = (c - 1.0);
    b = (b * 2);
    c = (c / 2.0);
    d = (b < 2);
    d = (b <= 3);
    d = (b > 4);
    d = (b >= 1);
    d = (b == 2);
    d = (d && d);
    d = (d || d);
    d = (!d);
    b = -1;
}
'''

    assert correct == out

def test_pir_cabana_visit_call():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    # TODO Check backend call functionality.

    def d():
        create_variable(c_int, x)
        mycall(x)
    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void d(){
    int x;
    mycall(x);
}
'''
    assert correct == out


def test_pir_cabana_visit_members():
    astpir = ast_to_pir_visitor()
    cpir = Cabana_PIR_Visitor(None)
    mystruc1 = StructureType()
    mystruc1.add("d", INT_TYPE)
    type_mapping_str["mystruc1"]=mystruc1
    def a():
        create_variable(mystruc1, d)
        d.d = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(a)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void a(){
    struct mystruc1 d;
    d.d = 1;
}
'''
    assert correct == out
    del(type_mapping_str["mystruc1"])

    mystruc2 = StructureType()
    substruc = StructureType()
    substruc.add("e", INT_TYPE)
    mystruc2.add("d", substruc)
    type_mapping_str["mystruc2"]=mystruc2
    def b():
        create_variable(mystruc2, b)
        b.d.e = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(b)))
    pir = astpir.visit(c)
    out = cpir(pir)

    correct = '''void b(){
    struct mystruc2 b;
    b.d.e = 1;
}
'''
    assert correct == out
    del(type_mapping_str["mystruc2"])

    mystruc2 = StructureType()
    substruc1 = ArrayType(INT_TYPE, [ArrayType.Extent.DYNAMIC])
    mystruc2.add("c", substruc1)
    type_mapping_str["mystruc2"]=mystruc2
    def d():
        create_variable(mystruc2, b)
        b.c[0] = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(d)))
    pir = astpir.visit(c)
    out = cpir(pir)
    correct = '''void d(){
    struct mystruc2 b;
    b.c[0] = 1;
}
'''
    assert correct == out
    del(type_mapping_str["mystruc2"])

def test_pir_cabana_pointer_reference():
    a = PointerReference(PointerSymbol("a", PointerType(INT_TYPE)))
    cpir = Cabana_PIR_Visitor(None)
    out = cpir(a)
    assert "*a" == out

def test_pir_cabana_pointer_symbol():
    a = PointerSymbol("a", PointerType(INT_TYPE))
    cpir = Cabana_PIR_Visitor(None)
    out = cpir.visit_pointersymbol_node(a)
    assert "int*" == out

def test_pir_cabana_arrayreference():
    cpir = Cabana_PIR_Visitor(None)
    ar = ArrayReference(ArraySymbol("x", ArrayType(INT_TYPE, [2,2])), [Literal("2", INT_TYPE), Literal("2", INT_TYPE)])
    out = cpir(ar)
    assert out == "x[2][2]"

def test_pir_cabana_config_reference():

    config = type_mapping_str["config"]
    backend = Cabana_PIR()
    cpir = Cabana_PIR_Visitor(backend)

    def x(arg: part, conf: config):
        conf.space = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    v = pir_perpart_visitor()
    pir = v.visit(c)
    out = cpir(pir)
    correct = '''struct x_functor{
    config_struct_type _conf;

    KOKKOS_INLINE_FUNCTION
     x_functor( config_struct_type conf):
    _conf(conf){}

    void operator()(const int i, const int a) const{
        _conf(0).space = 1;
    }
};
'''
    assert out == correct

    def main():
        config.space = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    v = pir_main_visitor()
    pir = v.visit(c)
    out = cpir(pir)

    correct = "config.config_host(0).space = 1;\n"
    assert correct in out


def test_pir_cabana_particle_references_and_perpart():
    backend = Cabana_PIR()
    cpir = Cabana_PIR_Visitor(backend)
    xs = StructureType()
    xs.add("boo", INT_TYPE)
    backend.add_structure(xs, "xs")
    type_mapping_str["abc"] = xs
    # Create a perpart_kernel
    v = pir_perpart_visitor()
    type_mapping_str["part"].add("subpart", INT_TYPE)
    new_type = ArrayType(INT_TYPE, [3])
    type_mapping_str["part"].add("subarray", new_type)
    def x(arg: part, c: c_int):
        arg.subpart = 1
        arg.subarray[0] = 1
        arg.core_part.position[0] = 2.0

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = v.visit(c)
    out = cpir(pir)
    correct = '''template < class SUBPART, class SUBARRAY, class CORE_PART_POSITION >
struct x_functor{
    config_struct_type _c;
    SUBPART _subpart;
    SUBARRAY _subarray;
    CORE_PART_POSITION _core_part_position;
    abc xs;

    KOKKOS_INLINE_FUNCTION
     x_functor( SUBPART subpart, SUBARRAY subarray, CORE_PART_POSITION core_part_position, config_struct_type c, abc XS):
    _subpart(subpart), _subarray(subarray), _core_part_position(core_part_position), xs(XS), _c(c){}

    void update_structs(abc XS){
        xs = XS;
    }

    void operator()(const int i, const int a) const{
        _subpart.access(i, a) = 1;
        _subarray.access(i, a, 0) = 1;
        _core_part_position.access(i, a, 0) = 2.0;
    }
};
'''
    assert correct == out
    del(type_mapping_str["abc"])

    assert pir == backend._kernels["x"]

    def z(arg: part, c: c_int):
        invoke(x)
    c = ast.parse(textwrap.dedent(inspect.getsource(z)))
    pir = v.visit(c)
    with pytest.raises(UnsupportedCodeError) as excinfo:
        out = cpir(pir)
    assert "Per part kernel cannot contain Invokes." in str(excinfo.value)

    def y(arg: part, c: c_int):
        create_variable(c_int, b)
        create_variable(c_int, f)
        b = c
        f = c + b
        arg.core_part.position[0] = f

    c = ast.parse(textwrap.dedent(inspect.getsource(y)))
    pir = v.visit(c)
    out = cpir(pir)
    correct = '''template < class CORE_PART_POSITION >
struct y_functor{
    config_struct_type _c;
    CORE_PART_POSITION _core_part_position;
    xs xs;

    KOKKOS_INLINE_FUNCTION
     y_functor( CORE_PART_POSITION core_part_position, config_struct_type c, xs XS):
    _core_part_position(core_part_position), xs(XS), _c(c){}

    void update_structs(xs XS){
        xs = XS;
    }

    void operator()(const int i, const int a) const{
        int b;
        int f;
        b = c;
        f = (c + b);
        _core_part_position.access(i, a, 0) = f;
    }
};
'''
    assert correct == out

    cpir._in_kernel = False
    fake = ParticlePositionReference(StructureSymbol("part",StructureType()), 0)
    with pytest.raises(NotImplementedError):
        cpir(fake)


def test_pir_cabana_pairwise():
    backend = Cabana_PIR()
    cpir = Cabana_PIR_Visitor(backend)
    # Create a pairwise kernel
    v = pir_pairwise_visitor()
    def x( arg1: part, arg2: part, c: c_int):
        arg1.subpart = 1

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = v.visit(c)
    with pytest.raises(NotImplementedError) as excinfo:
        out = cpir(pir)
    assert "Pairwise kernels not yet supported in cabana." in str(excinfo.value)

def test_pir_cabana_mainkernel():
    backend = Cabana_PIR()
    cpir = Cabana_PIR_Visitor(backend)
    # Create a main kernel
    v = pir_main_visitor()

    def main():
        create_variable(c_int, a)
        while a < 100:
            a = a + 1

    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    pir = v.visit(c)
    out = cpir(pir)

    correct = '''int main( int argc, char* argv[] ){
    int a;
    while((a < 100)){
        a = (a + 1);
    }
}
'''
    assert correct == out

    # Create a perpart_kernel
    v2 = pir_perpart_visitor()
    reset_part_and_config()
    type_mapping_str["part"].add("subpart", INT_TYPE)
    # subpart still exists from previous tests.
    def x(arg: part, c: c_int):
        arg.subpart = 1
        arg.core_part.position[0] = 2.0

    c = ast.parse(textwrap.dedent(inspect.getsource(x)))
    pir = v2.visit(c)
    out = cpir(pir)

    def main():
        invoke(x)
    c = ast.parse(textwrap.dedent(inspect.getsource(main)))
    pir = v.visit(c)
    xs = StructureType()
    xs.add("boo", INT_TYPE)
    backend.add_structure(xs, "xs")
    out = cpir(pir)
    correct = '''int main( int argc, char* argv[] ){
    Kokkos::deep_copy(config.config, config.config_host);
    x.update_structs(xs);
    Cabana::simd_parallel_for(simd_policy, x, "x");
    Kokkos::fence();
}
'''
    assert correct == out

    set_mpi(True)
    pir = v.visit(c)
    out = cpir(pir)
    correct = '''int main( int argc, char* argv[] ){
    Kokkos::deep_copy(config.config, config.config_host);
    x.update_structs(xs);
    Cabana::simd_parallel_for(simd_policy, x, "x");
    Kokkos::fence();
    {
        _rank_update_functor<decltype(core_part_position_slice), decltype(neighbour_part_rank_slice)> ruf(
            config.config_host(0).space.box_dims, core_part_position_slice, neighbour_part_rank_slice,
            config.config_host(0).space.box_dims.x_ranks,
            config.config_host(0).space.box_dims.y_ranks,
            config.config_host(0).space.box_dims.z_ranks, myrank);
        Cabana::SimdPolicy<VectorLength, ExecutionSpace> temp_policy(0, particle_aosoa.size());
        Cabana::simd_parallel_for( temp_policy, ruf, "rank_update_functor");
        Kokkos::fence();
    }

    Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
    _migrator.exchange_data(particle_aosoa_host, neighbors, myrank, particle_aosoa_host.size());
    particle_aosoa.resize(particle_aosoa_host.size());
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
    core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    neighbour_part_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa);
    neighbour_part_old_position_slice = Cabana::slice<neighbour_part_old_position>(particle_aosoa);
    simd_policy = Cabana::SimdPolicy<VectorLength, ExecutionSpace>(0, particle_aosoa.size());

    x = x_functor<decltype(subpart_slice), decltype(core_part_position_slice)>(subpart_slice, core_part_position_slice, config.config,xs);

}
'''
    set_mpi(False)
    assert correct == out
    # TODO Add structure

    set_mpi(True)
    backend.set_boundary_condition(periodic_boundaries)
    pir = v.visit(c)
    out = cpir(pir)
    set_mpi(False)
    correct = '''int main( int argc, char* argv[] ){
    Kokkos::deep_copy(config.config, config.config_host);
    x.update_structs(xs);
    Cabana::simd_parallel_for(simd_policy, x, "x");
    Kokkos::fence();
    {
        _rank_update_functor<decltype(core_part_position_slice), decltype(neighbour_part_rank_slice)> ruf(
            config.config_host(0).space.box_dims, core_part_position_slice, neighbour_part_rank_slice,
            config.config_host(0).space.box_dims.x_ranks,
            config.config_host(0).space.box_dims.y_ranks,
            config.config_host(0).space.box_dims.z_ranks, myrank);
        Cabana::SimdPolicy<VectorLength, ExecutionSpace> temp_policy(0, particle_aosoa.size());
        Cabana::simd_parallel_for( temp_policy, ruf, "rank_update_functor");
        Kokkos::fence();
    }

    x.update_structs(xs);
    Cabana::simd_parallel_for(simd_policy, periodic_boundaries, "periodic_boundaries");
    Kokkos::fence();
    Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
    _migrator.exchange_data(particle_aosoa_host, neighbors, myrank, particle_aosoa_host.size());
    particle_aosoa.resize(particle_aosoa_host.size());
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
    core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    neighbour_part_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa);
    neighbour_part_old_position_slice = Cabana::slice<neighbour_part_old_position>(particle_aosoa);
    simd_policy = Cabana::SimdPolicy<VectorLength, ExecutionSpace>(0, particle_aosoa.size());

    x = x_functor<decltype(subpart_slice), decltype(core_part_position_slice)>(subpart_slice, core_part_position_slice, config.config,xs);

}
'''
    assert correct == out
