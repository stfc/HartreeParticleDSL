import ast
import inspect
import pytest

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSL import Particle, Config, _HartreeParticleDSL, set_mpi, reset_for_tests
from HartreeParticleDSL.HartreeParticleDSL import part, config

from HartreeParticleDSL.backends.Cabana_PIR_backend.cabana_pir import Cabana_PIR
from HartreeParticleDSL.backends.Cabana_PIR_backend.pir_to_cabana_visitor import Cabana_PIR_Visitor

from HartreeParticleDSL.coupled_systems.base_coupler.base_coupler import *

from HartreeParticleDSL.IO_modules.random_IO.random_IO import *
from HartreeParticleDSL.IO_modules.IO_Exceptions import *

from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType, INT_TYPE, type_mapping_str, \
        ArrayType, PointerType

from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.kernels import PairwiseKernel, PerPartKernel
import HartreeParticleDSL.kernel_types.kernels as kernels

from HartreeParticleDSL.HartreeParticleDSLExceptions import InvalidNameError, \
                                                            UnsupportedTypeError, \
                                                            InternalError, \
                                                            IRGenerationError


from HartreeParticleDSL.inbuilt_kernels.boundaries.periodic_boundaries import periodic_boundaries


class coupler_test(base_coupler):
    def __init__(self):
        pass

    def a_function(self, current_indent=0, indent=0):
        return "test_string()"

    def b_function(self, arg):
        return arg

    def get_includes(self):
        return ["\"test.h\""]

    def get_includes_header(self):
        return []

    def copy_files(self):
        pass

    def setup_testcase(self, filename, current_indent=0):
        return current_indent * " " + "setup_testcase();\n"

    def has_preferred_decomposition(self):
        return True

    def get_preferred_decomposition(self,arg1, current_indent=0):
        return current_indent * " " + "preferred_decomposition()"

def test_cabana_pir_init():
    a = Cabana_PIR()
    assert "<Cabana_Core.hpp>" in a._includes
    assert "<Kokkos_Core.hpp>" in a._includes
    assert "<Kokkos_Random.hpp>" in a._includes
    assert "<iostream>" in a._includes
    assert "<cmath>" in a._includes
    assert "<cstdio>" in a._includes
    assert "\"part.h\"" in a._includes
    assert "<math.h>" in a._includes

def test_cabana_pir_register_kernel():
    a = Cabana_PIR()
    Kern1 = Kern()
    Kern2 = Kern()
    a.register_kernel("x", Kern1)
    assert a._kernels["x"] == Kern1
    with pytest.raises(InvalidNameError) as excinfo:
        a.register_kernel("x", Kern2)
    assert "Kernel with name x already exists." in str(excinfo.value)

    a.register_kernel("y", Kern2)
    assert a._kernels["y"] == Kern2

def test_cabana_pir_add_structure():
    a = Cabana_PIR()
    x = StructureType()
    a.add_structure(x, "x")
    assert a._structures["x"] == x

def test_cabana_pir_add_kernel_slice():
    a = Cabana_PIR()
    a.add_kernel_slices("mykern", "slice1")
    assert "slice1" in  a._kernel_slices["mykern"]

    a.add_kernel_slices("mykern", "slice1")
    assert len(a._kernel_slices.keys()) == 1
    assert len(a._kernel_slices["mykern"]) == 1

    a.add_kernel_slices("mykern", "slice2")
    assert len(a._kernel_slices["mykern"]) == 2

def test_cabana_pir_add_type():
    backend = Cabana_PIR()
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_type("t", 123)
    assert ("Attempting to add a new type but got <class 'int'> "
            "but expected a Particle IR DataType instead." in str(excinfo.value))
    backend.add_type("t", INT_TYPE)
    assert type_mapping_str["t"] == INT_TYPE
    with pytest.raises(InvalidNameError) as excinfo:
        backend.add_type("t", INT_TYPE)
    assert ("Attempting to add new type t but a type with that name already "
            "exists." in str(excinfo.value))
    del type_mapping_str["t"]

def test_cabana_pir_create_global_variable():
    backend = Cabana_PIR()
    backend.create_global_variable(INT_TYPE, "hello", "0")
    assert backend._global_values["hello"] == "0"
    assert HartreeParticleDSL.global_symbol_table().lookup("hello").datatype == INT_TYPE
    with pytest.raises(IRGenerationError) as excinfo:
        backend.create_global_variable(INT_TYPE, "hello", "0")
    assert "Tried to add a new symbol hello but it was already present in the symbol table." in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        backend.create_global_variable("stringtype", "hello", "0")
    assert ("Attempting to create global variable but c_type argument is not "
            "a supported datatype. Got <class 'str'>" in str(excinfo.value))
    pt = PointerType(INT_TYPE)
    backend.create_global_variable(pt, "hello2", "&a")
    assert backend._global_values["hello2"] == "&a"
    assert HartreeParticleDSL.global_symbol_table().lookup("hello2").datatype == pt
    at = ArrayType(INT_TYPE, [2,2])
    backend.create_global_variable(at, "hello3")
    assert backend._global_values["hello3"] is None
    assert HartreeParticleDSL.global_symbol_table().lookup("hello3").datatype == at
    st = StructureType()
    backend.create_global_variable(st, "hello4")
    assert backend._global_values["hello4"] is None
    assert HartreeParticleDSL.global_symbol_table().lookup("hello4").datatype == st

def test_cabana_pir_set_io_modules():
    mod = Random_Particles()
    backend = Cabana_PIR()
    backend.set_io_modules(mod, mod)
    assert backend._input_module is mod
    assert backend._output_module is mod
    class temp_module(IO_Module):
        def __init__(self):
            pass
    mod2 = temp_module()
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod2, mod)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)
    with pytest.raises(InvalidIOModuleError) as excinfo:
        backend.set_io_modules(mod, mod2)
    assert "Cabana_PIR backend does not support temp_module" in str(excinfo.value)

def test_cabana_pir_add_include():
    backend = Cabana_PIR()
    backend.add_include("<string.h>")
    assert "<string.h>" in backend._includes
    backend.add_include("<string.h>")
    assert backend._includes.count("<string.h>") is 1

def test_cabana_pir_generate_includes():
    backend = Cabana_PIR()
    mod = Random_Particles()
    backend.set_io_modules(mod, None)
    strin = backend.generate_includes()
    assert "<Cabana_Core.hpp>" in strin
    assert "random_io_cabana.hpp" in strin
    assert "<Cabana_Core.hpp>" in strin
    assert "\"part.h\"" in strin
    assert "<Kokkos_Random.hpp>" in strin
    assert "<iostream>" in strin
    assert "<cmath>" in strin
    assert "<cstdio>" in strin
    assert "<math.h>" in strin
    backend = Cabana_PIR()
    backend.set_io_modules(None, mod)
    strin = backend.generate_includes()
    assert "<Cabana_Core.hpp>" in strin
    assert "random_io_cabana.hpp" in strin
    assert "<Cabana_Core.hpp>" in strin
    assert "\"part.h\"" in strin
    assert "<Kokkos_Random.hpp>" in strin
    assert "<iostream>" in strin
    assert "<cmath>" in strin
    assert "<cstdio>" in strin
    assert "<math.h>" in strin

class coupler_test2(base_coupler):
    def __init__(self):
        pass

    def a_function(self):
        return "test_string"

    def b_function(self, arg):
        return arg

    def get_includes(self):
        return ["a"]

    def get_includes_header(self):
        return ["a"]
    
    def has_preferred_decomposition(self):
        return False

    def setup_testcase(self, filename, current_indent=0):
        return current_indent*" " + "setup2()"

    def get_extra_symbols(self, flist):
        return "x"

class dummy_module(IO_Module, Cabana_PIR_IO_Mixin):

    def __init__(self):
        pass

    def gen_code_cabana(self, part):
        return "hello"

    def get_header_includes_cabana_pir(self):
        return []


def test_cabana_pir_gen_headers(capsys):
    backend = Cabana_PIR()
    part = Particle()
    config = Config()
    coupler = coupler_test2()
    backend.add_coupler(coupler)
    mod = dummy_module()
    backend.set_io_modules(mod, mod)
    _HartreeParticleDSL.the_instance = None
    HartreeParticleDSL.set_mpi(True)
    backend.create_global_variable(INT_TYPE, "hello", "0")
    backend.create_global_variable(INT_TYPE, "hello2")
    backend.gen_headers(config, part)
    f_str = ""
    with open('part.h', 'r') as f:
        f_str = f.readlines()
    f_str = "".join(f_str)
    correct = '''#ifndef PART_H
#define PART_H
#include <Kokkos_Core.hpp>
#include <Cabana_Core.hpp>
#include a
using MemorySpace = Kokkos::HostSpace;
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using DeviceType = Kokkos::Device<Kokkos::DefaultExecutionSpace, MemorySpace>;
using HostType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
const int VectorLength = 16;

extern int hello;
extern int hello2;
struct boundary{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
    double local_x_min, local_x_max;
    double local_y_min, local_y_max;
    double local_z_min, local_z_max;
    int x_ranks, y_ranks, z_ranks;
};

struct space_type{
    boundary box_dims;
    int nparts;
};

struct neighbour_config_type{
};

struct config_view_type{
        struct space_type space;
        struct neighbour_config_type neighbour_config;
};

using config_struct_type = Kokkos::View<struct config_view_type*, MemorySpace>;
using config_struct_host = config_struct_type::HostMirror;
struct config_type{
    config_struct_type config;
    config_struct_host config_host;
};
enum FieldNames{core_part_velocity = 0,
                 core_part_position,
                 neighbour_part_cutoff,
                 neighbour_part_deletion_flag,
                 neighbour_part_rank,
                 neighbour_part_old_position
               };
using DataTypes = Cabana::MemberTypes<double[3],
    double[3],
    double,
    int,
    int,
    double[3]>;


template<class aosoa> class Migrator{

    private:
        Kokkos::View<double***, MemorySpace> pos_space;
        Kokkos::View<double***, MemorySpace> vel_space;
        Kokkos::View<double**, MemorySpace> cutoff_space;
        int _buffer_size;

    public:
        Migrator(int buffer_size, int nr_neighbours){
            _buffer_size = buffer_size;
            pos_space = Kokkos::View<double***, MemorySpace>("temp_pos", nr_neighbours, buffer_size, 3);
            vel_space = Kokkos::View<double***, MemorySpace>("temp_velocity", nr_neighbours, buffer_size, 3);
            cutoff_space = Kokkos::View<double**, MemorySpace>("temp_cutoff", nr_neighbours, buffer_size);
        }

    void exchange_data( aosoa &particle_aosoa, std::vector<int> neighbors, int myrank, int npart){

        auto rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, "rank");
        auto last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, "last_pos");
        auto pos_s = Cabana::slice<core_part_position>(particle_aosoa, "position");
        auto vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, "velocity");
        auto cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, "cutoff");
        int *send_count = (int*) malloc(sizeof(int) * neighbors.size());
        int count_neighbours = 0;
        int end = particle_aosoa.size() - 1;
        for(int i = 0; i < neighbors.size(); i++){
                send_count[i] = 0;
        }

        for(int i = particle_aosoa.size()-1; i>=0; i--){
            if(rank_slice(i) != myrank && rank_slice(i) >= 0){
                int therank = rank_slice(i);
                for(int k = 0; k < neighbors.size(); k++){
                    if(therank == neighbors[k]){
                        therank = k;
                        break;
                    }
                }
            int pos = send_count[therank];
            pos_space(therank, pos, 0) = pos_s(i, 0);
            pos_space(therank, pos, 1) = pos_s(i, 1);
            pos_space(therank, pos, 2) = pos_s(i, 2);
            vel_space(therank, pos, 0) = vel_s(i, 0);
            vel_space(therank, pos, 1) = vel_s(i, 1);
            vel_space(therank, pos, 2) = vel_s(i, 2);
            cutoff_space(therank, pos) = cutoff_s(i);
            send_count[therank]++;

            while(rank_slice(end) != myrank && end > 0){
                end--;
            }
            if(end > i){
                rank_slice(i) = rank_slice(end);
                pos_s(i, 0) = pos_s(end, 0);
                pos_s(i, 1) = pos_s(end, 1);
                pos_s(i, 2) = pos_s(end, 2);
                vel_s(i, 0) = vel_s(end, 0);
                vel_s(i, 1) = vel_s(end, 1);
                vel_s(i, 2) = vel_s(end, 2);
                cutoff_s(i) = cutoff_s(end);
                rank_slice(end) = -1;
            }else{
                rank_slice(i) = -1;
                end++;
            }
            continue;
        }

    }
    // Data collected, need to send information to neighbours to know what to expect
    int *recv_count = (int*) malloc(sizeof(int) * neighbors.size());
    MPI_Request *requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2);
    int req_num = 0;
    for(int i = 0; i < neighbors.size(); i++){
        recv_count[i] = 0;
        if(neighbors[i] == myrank){
            continue;
        }
        MPI_Irecv(&recv_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);
        MPI_Isend(&send_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);
    }
    MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    int total_size = 0;
    for(int i = 0; i < neighbors.size(); i++){
             total_size += recv_count[i];
    }
    Kokkos::View<double***, MemorySpace> r_pos_space("temp_pos", neighbors.size(), total_size, 3);
    Kokkos::View<double***, MemorySpace> r_vel_space("temp_vel", neighbors.size(), total_size, 3);
    Kokkos::View<double**, MemorySpace> r_cutoff_space("temp_cutoff", neighbors.size(), total_size);

    free(requests);
    requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2 * 3);
    req_num = 0;
    int tag = 0;
    for(int i = 0; i < neighbors.size(); i++){
        if(neighbors[i] != myrank){
            tag = 0;
            MPI_Irecv(&r_pos_space.data()[r_pos_space.extent(1)*r_pos_space.extent(2)*i], recv_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_vel_space.data()[r_vel_space.extent(1)*r_vel_space.extent(2)*i], recv_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_cutoff_space.data()[r_cutoff_space.extent(1)*i], recv_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            tag = 0;
            MPI_Isend(&pos_space.data()[pos_space.extent(1)*pos_space.extent(2)*i], send_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&vel_space.data()[vel_space.extent(1)*vel_space.extent(2)*i], send_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&cutoff_space.data()[cutoff_space.extent(1)*i], send_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
        }
    }
    MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
    free(requests);
    int recvd = 0;
    int sent = 0;
    for(int i = 0; i < neighbors.size(); i++){
        recvd += recv_count[i];
        sent += send_count[i];
    }
    int size_change = recvd - sent;
    int current_size =  particle_aosoa.size();
    if(size_change != 0){
        particle_aosoa.resize(current_size+size_change);
    }
    auto new_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, "new_rank");
    for(int i = particle_aosoa.size() - 1; i > end; i--){
        new_rank_slice(i) = -1;
    }
    if(size_change > 0){
        if(sent = 0){
            end = current_size;
        }
        while(end < current_size && end < particle_aosoa.size() && new_rank_slice(end) != -1 ) end++;
        for(int i = 0; i < particle_aosoa.size(); i++){
            new_rank_slice(i) = -1;
        }
    }
    auto new_last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, "new_last_pos");
    auto new_pos_s = Cabana::slice<core_part_position>(particle_aosoa, "new_position");
    auto new_vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, "new_velocity");
    auto new_cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, "new_cutoff");
    int x = 0;
    for(int j = 0; j < neighbors.size(); j++){
        for(int i = 0; i < recv_count[j]; i++){
            new_pos_s(end+x, 0) = r_pos_space(j,i,0);
            new_pos_s(end+x, 1) = r_pos_space(j,i,1);
            new_pos_s(end+x, 2) = r_pos_space(j,i,2);
            new_vel_s(end+x, 0) = r_vel_space(j,i,0);
            new_vel_s(end+x, 1) = r_vel_space(j,i,1);
            new_vel_s(end+x, 2) = r_vel_space(j,i,2);
            new_cutoff_s(end+x) = r_cutoff_space(j,i);
            x++;
        }
    }
    free(recv_count);
    free(send_count);
}};

KOKKOS_INLINE_FUNCTION
int get_oneD_rank(int x_r, int y_r, int z_r, int x_ranks, int y_ranks, int z_ranks){
    int oneD_rank = z_r*(x_ranks*y_ranks) + y_r*(x_ranks) + x_r;
    return oneD_rank;
}
KOKKOS_INLINE_FUNCTION
void get_threeD_rank(int rank, int *x, int *y, int *z, int x_ranks, int y_ranks, int z_ranks){
    int z_r = rank / (x_ranks*y_ranks);
    int y_r = (rank - z_r*x_ranks*y_ranks) / x_ranks;
    int x_r = rank - z_r*x_ranks*y_ranks - y_r*x_ranks;
    *x = x_r;
    *y = y_r;
    *z = z_r;
};

template<class PartPosSlice, class RankSlice>
struct _rank_update_functor{
    boundary _box;
    PartPosSlice _part_pos;
    RankSlice _rank;
    int _myrank;
    int _xranks;
    int _yranks;
    int _zranks;
    int _local_x_rank;
    int _local_y_rank;
    int _local_z_rank;
    int _nranks;

    KOKKOS_INLINE_FUNCTION
    _rank_update_functor(boundary box, PartPosSlice pos, RankSlice rank, int xranks, int yranks, int zranks, int myrank):
        _box(box), _part_pos(pos), _rank(rank), _xranks(xranks), _yranks(yranks), _zranks(zranks), _myrank(myrank){
        get_threeD_rank(myrank, &_local_x_rank, &_local_y_rank, &_local_z_rank, xranks, yranks, zranks);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix, const int ij) const{
        int xr, yr, zr;
        get_threeD_rank(_myrank, &xr, &yr, &zr, _xranks, _yranks, _zranks);
        xr = _local_x_rank;
        yr = _local_y_rank;
        zr = _local_z_rank;
        if(_part_pos.access(ix, ij, 0) >= _box.local_x_max){
            xr = xr + 1;
            if( xr >= _xranks ) xr = 0;
        }
        if(_part_pos.access(ix, ij, 0) < _box.local_x_min){
            xr = xr - 1;
            if( xr < 0 ) xr = _xranks-1;
        }
        if(_part_pos.access(ix, ij, 1) >= _box.local_y_max){
            yr = yr + 1;
            if( yr >= _yranks ) yr = 0;
        }
        if(_part_pos.access(ix, ij, 1) < _box.local_y_min){
            yr = yr - 1;
            if( yr < 0 ) yr = _yranks-1;
        }
        if(_part_pos.access(ix, ij, 2) >= _box.local_z_max){
            zr = zr + 1;
            if( zr >= _zranks ) zr = 0;
        }
        if(_part_pos.access(ix, ij, 2) < _box.local_z_min){
            zr = zr - 1;
            if( zr < 0 ) zr = _zranks-1;
        }
        _rank.access(ix, ij) = get_oneD_rank(xr, yr, zr, _xranks, _yranks, _zranks);
    }};

#endif'''
    HartreeParticleDSL.set_mpi(False)
    assert f_str == correct

    HartreeParticleDSL.set_cuda(True)
    with pytest.raises(NotImplementedError) as excinfo:
        backend.gen_headers(config, part)
    assert "Non-constant global variables are not supported in Cabana_PIR with CUDA enabled." in str(excinfo.value)

    backend = Cabana_PIR()
    backend.gen_headers(config, part)
    f_str = ""
    with open('part.h', 'r') as f:
        f_str = f.readlines()
    f_str = "".join(f_str)
    assert "using MemorySpace = Kokkos::CudaSpace;" in f_str
    HartreeParticleDSL.set_cuda(False)

def kern2(part: part, config: config):
    create_variable(c_double, a, 2.0)
    part.core_part.position[0] = part.core_part.position[0] + 2.0
    part.x = 2

def kern3(part1: part, part2: part, config:config):
    part1.core_part.position[0] = part2.core_part.position[0]

def test_cabana_pir_gen_kernel():
    _HartreeParticleDSL.the_instance = None
    backend = Cabana_PIR()
    config = Config()
    part = Particle()
    part.add_element("x", "int")
    backend.gen_headers(config, part)
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    assert backend._kernels_pir["kern2"] is not None

    kernel2 = kernels.pairwise_interaction(kern3) 
    backend.gen_kernel(kernel2)
    assert backend._kernels_pir["kern3"] is not None

    with pytest.raises(NotImplementedError) as excinfo:
        backend.gen_kernel("oops")
    assert ("Cabana PIR backend doesn't yet support kernel of type <class 'str'>."
            in str(excinfo.value))

def main():
    initialise(particle_count=0, filename="abd")
    create_variable(c_int, a, 0)
    a = a + 1
    invoke(kern2)
    cleanup()

class random_particles_test(Random_Particles):
    def gen_code_cabana_pir(self, a):
        return "void func(){}"

def test_print_main():
    '''Test the print_main function of Cabana'''
    _HartreeParticleDSL.the_instance = None
    backend = Cabana_PIR()
    backend.add_coupler(coupler_test())
    backend.set_boundary_condition(periodic_boundaries)
    HartreeParticleDSL.set_backend(backend)
    mod = random_particles_test()
    mod2 = random_particles_test()
    backend.set_io_modules(mod, mod2)
    config = Config()
    part = Particle()
    part.add_element("x", "int")
    backend.gen_headers(config, part)
    kernel = kernels.perpart_interaction(kern2)
    backend.gen_kernel(kernel)
    m = ast.parse(inspect.getsource(main))
    HartreeParticleDSL.set_mpi(False)
    backend.create_global_variable(INT_TYPE, "hello1", "0")
    backend.create_global_variable(INT_TYPE, "hello2")
    backend.print_main(m)
    with open('code.cpp') as f:
        out = f.read()
        correct ='''int main( int argc, char* argv[] ){
    int a;
    Kokkos::ScopeGuard scope_guard(argc, argv);
{
    config_type config;
    config.config = config_struct_type("config", 1);
    config.config_host = Kokkos::create_mirror_view(config.config);
    int myrank = 0;
    int nranks = 1;
    setup_testcase();
    Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 0);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 0);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);

    Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());
    auto core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    auto neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    auto neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);
    auto x_slice = Cabana::slice<x>(particle_aosoa);
    kern2_functor<decltype(core_part_position_slice), decltype(x_slice)> kern2(core_part_position_slice, x_slice, config.config);
    periodic_boundaries_functor<decltype(core_part_position_slice)> periodic_boundaries(core_part_position_slice, config.config);

    a = 0;
    a = (a + 1);
    Kokkos::deep_copy(config.config, config.config_host);
    Cabana::simd_parallel_for(simd_policy, kern2, "kern2");
    Kokkos::fence();
    Cabana::simd_parallel_for(simd_policy, periodic_boundaries, "periodic_boundaries");
    Kokkos::fence();
    }

}'''
        assert correct in out
        correct = '''template < class CORE_PART_POSITION, class X >
struct kern2_functor{
    config_struct_type _config;
    CORE_PART_POSITION _core_part_position;
    X _x;

    KOKKOS_INLINE_FUNCTION
     kern2_functor( CORE_PART_POSITION core_part_position, X x, config_struct_type config):
    _core_part_position(core_part_position), _x(x), _config(config){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, const int a) const{
        double a;
        a = 2.0;
        _core_part_position.access(i, a, 0) = (_core_part_position.access(i, a, 0) + 2.0);
        _x.access(i, a) = 2;
    }
};'''
        assert correct in out

        correct = '''int hello1 = 0;
int hello2;
void func(){}
void func(){}'''
        assert correct in out
    HartreeParticleDSL.set_cuda(True)
    with pytest.raises(NotImplementedError) as excinfo:
        backend.print_main(m)
    assert "Non-constant global variables are not supported in Cabana_PIR with CUDA enabled." in str(excinfo.value)
    HartreeParticleDSL.set_cuda(False)



def test_cabana_pir_set_cutoff():
    a = Cabana_PIR()
    with pytest.raises(NotImplementedError) as excinfo:
        a.set_cutoff(0)
    assert ("Cabana PIR backend doesn't yet support pairwise interactions" in
            str(excinfo.value))

def test_cabana_pir_initialisation_code():
    mod = Random_Particles()
    backend = Cabana_PIR()
    backend.set_io_modules(mod, mod)
    val = backend.initialisation_code(1000, None)
    correct = '''Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 1000);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 1000);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
'''
    assert val == correct

def test_cabana_pir_gen_particle():
    backend = Cabana_PIR()
    part = Particle()
    part.add_element("thing", "double[4]")
    part.add_element("thing2", "double")
    rval = backend.gen_particle(part)
    correct = '''enum FieldNames{thing = 0,
                 thing2,
                 core_part_velocity,
                 core_part_position,
                 neighbour_part_cutoff,
                 neighbour_part_deletion_flag
               };
using DataTypes = Cabana::MemberTypes<double[4],
    double,
    double[3],
    double[3],
    double,
    int>;
'''
    assert correct == rval

    HartreeParticleDSL.set_mpi(True)
    rval = backend.gen_particle(part)
    HartreeParticleDSL.set_mpi(False)
    assert ",\n                 neighbour_part_rank" in rval
    assert ",\n                 neighbour_part_old_position" in rval
    assert ",\n    int,\n    double[3]" in rval


def test_cabana_pir_gen_config():
    backend = Cabana_PIR()
    conf = Config()
    conf.add_element("thing", "double[4]")
    rval = backend.gen_config(conf)
    correct = '''struct boundary{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
};

struct space_type{
    boundary box_dims;
    int nparts;
};

struct neighbour_config_type{
};

struct config_view_type{
        struct space_type space;
        struct neighbour_config_type neighbour_config;
        double thing[4];
};

using config_struct_type = Kokkos::View<struct config_view_type*, MemorySpace>;
using config_struct_host = config_struct_type::HostMirror;
struct config_type{
    config_struct_type config;
    config_struct_host config_host;
};
'''
    assert correct == rval
    HartreeParticleDSL.set_mpi(True)
    rval = backend.gen_config(conf)
    HartreeParticleDSL.set_mpi(False)
    correct = '''    double local_x_min, local_x_max;
    double local_y_min, local_y_max;
    double local_z_min, local_z_max;
    int x_ranks, y_ranks, z_ranks;
'''
    assert correct in rval

def test_cabana_pir_cleanup():
    backend = Cabana_PIR()
    st = backend.cleanup(current_indent=0)
    assert st == "\n}\n"

def test_cabana_pir_println():
    backend = Cabana_PIR()
    statement = backend.println("\"string 123\"", "total", current_indent=0)
    assert statement == "std::cout << \"string 123\" << total << \"\\n\";\n"

def test_cabana_pir_initialise():
    mod = Random_Particles()
    backend = Cabana_PIR()
    part = Particle()
    config = Config()
    backend.set_io_modules(mod, mod)
    backend.gen_particle(part)
    backend.gen_config(config)
    backend.add_structure("c_int", "mystruct")
    temp = StructureType()
    backend.add_structure(temp, "mystruct2")
    type_mapping_str["test"] = temp
    backend._kernel_slices["slice1"] = ["s1", "s2"]
    out = backend.initialise(123, "myfile", 4)
    correct = '''    Kokkos::ScopeGuard scope_guard(argc, argv);
{
    config_type config;
    config.config = config_struct_type("config", 1);
    config.config_host = Kokkos::create_mirror_view(config.config);
    mystruct mystruct;
    test mystruct2;
    int myrank = 0;
    int nranks = 1;
    Cabana::AoSoA<DataTypes, DeviceType, VectorLength> particle_aosoa( "particle_list", 123);
    Cabana::AoSoA<DataTypes, HostType, VectorLength> particle_aosoa_host( "particle_list_host", 123);
    random_io<decltype(particle_aosoa_host)>(particle_aosoa_host, config);
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);

    Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd_policy( 0, particle_aosoa.size());
    auto core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    auto core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    auto neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    auto neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);
    slice1_functor<decltype(s1_slice), decltype(s2_slice)> slice1(s1_slice, s2_slice, config.config, mystruct, mystruct2);
'''
    assert correct == out

    HartreeParticleDSL.set_mpi(True)
    with pytest.raises(NotImplementedError) as excinfo:
        rval = backend.initialise(123, "myfile", 4)
    assert ("Can't yet do a decomposition when coupled systems didn't do it "
            "for us." in str(excinfo.value))
    coupler = coupler_test()
    backend.add_coupler(coupler)
    coupler2 = coupler_test2()
    backend.add_coupler(coupler2)
    out = backend.initialise(123, "myfile", 4)
    HartreeParticleDSL.set_mpi(False)
    assert "int _provided;\n" in out
    assert "MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &_provided  );\n" in out
    assert "preferred_decomposition()" in out
    assert "setup2()" in out

    backend.add_coupler(coupler_test())
    HartreeParticleDSL.set_mpi(True)
    with pytest.raises(NotImplementedError) as excinfo:
         rval = backend.initialise(123, "myfile", 4) 
    HartreeParticleDSL.set_mpi(False)
    assert ("Can't handle multiple coupled systems with a preferred "
            "decomposition" in str(excinfo.value))

def test_cabana_pir_call_language_function():
    backend = Cabana_PIR()
    with pytest.raises(AttributeError):
        rval1 = backend.call_language_function("a_c_call", "*part", "20")
    rval2 = backend.call_language_function("cleanup", current_indent=2, indent=1)
    assert rval2 == "\n}\n"
    with pytest.raises(AttributeError):
        rval3 = backend.call_language_function("a_c_call", "*part", "20", current_indent=4, indent=1)
    with pytest.raises(NotImplementedError) as excinfo:
        val4 = backend.call_language_function("set_cutoff", "0.5", "C_AOS.CONSTANT")
    assert "Cabana PIR backend doesn't yet support pairwise interactions" in str(excinfo.value)

def test_cabana_pir_add_coupler():
    backend = Cabana_PIR()
    coupler = coupler_test()
    backend.add_coupler(coupler)
    assert coupler in backend._coupled_systems
    with pytest.raises(UnsupportedTypeError) as excinfo:
        backend.add_coupler(32)
    assert ("Can only couple to base_coupler classes or "
            "subclasses. Found int") in str(excinfo.value)
    call = Call("a_function")
    test_visitor = Cabana_PIR_Visitor(backend)
    assert test_visitor._visit(call) == "test_string()"
    call = Call("not_a_function")
    assert test_visitor._visit(call) == "not_a_function()"


def test_cabana_pir_write_output():
    backend = Cabana_PIR()
    class temp_module(IO_Module, Cabana_PIR_IO_Mixin):
        def __init__(self):
            pass
        def call_output_cabana_pir(self, num_parts, filename, variable, current_indent=4):
            return "Success"
    a = temp_module()
    mod = Random_Particles()
    backend.set_io_modules(mod, a)
    assert backend.write_output("a") == "Success\n"


def test_cabana_pir_boundary_conditions():
    backend = Cabana_PIR()
    with pytest.raises(TypeError) as excinfo:
        backend.set_boundary_condition(3)
    assert ("Cannot set boundary condition to a non perpart kernel" 
             in str(excinfo.value))

    assert backend._boundary_condition is None
    assert backend._boundary_condition_tree is None

    backend.set_boundary_condition(periodic_boundaries)

    assert backend._boundary_condition is periodic_boundaries
    assert backend._boundary_condition_tree == periodic_boundaries.get_kernel_tree()

    assert isinstance(backend.boundary_condition, PerPartKernel)

def test_mpi_headers():
    backend = Cabana_PIR()
    part = Particle()
    config = Config()
    part.add_element("testvar", "int[4][2]")
    part.add_element("testvar2", "double")
    part.add_element("testvar3", "float")
    part.add_element("testvar4", "int64_t")
    config.add_element("testcar", "double[4][2]")
    config.add_element("testcar2", "int")
    out = backend.mpi_headers(config, part)
    correct = '''

template<class aosoa> class Migrator{

    private:
        Kokkos::View<int****, MemorySpace> testvar_space;
        Kokkos::View<double**, MemorySpace> testvar2_space;
        Kokkos::View<float**, MemorySpace> testvar3_space;
        Kokkos::View<int64_t**, MemorySpace> testvar4_space;
        Kokkos::View<double***, MemorySpace> pos_space;
        Kokkos::View<double***, MemorySpace> vel_space;
        Kokkos::View<double**, MemorySpace> cutoff_space;
        int _buffer_size;

    public:
        Migrator(int buffer_size, int nr_neighbours){
            _buffer_size = buffer_size;
            testvar_space = Kokkos::View<int****, MemorySpace>("testvar_id", nr_neighbours, buffer_size, 4, 2);
            testvar2_space = Kokkos::View<double**, MemorySpace>("testvar2_id", nr_neighbours, buffer_size);
            testvar3_space = Kokkos::View<float**, MemorySpace>("testvar3_id", nr_neighbours, buffer_size);
            testvar4_space = Kokkos::View<int64_t**, MemorySpace>("testvar4_id", nr_neighbours, buffer_size);
            pos_space = Kokkos::View<double***, MemorySpace>("temp_pos", nr_neighbours, buffer_size, 3);
            vel_space = Kokkos::View<double***, MemorySpace>("temp_velocity", nr_neighbours, buffer_size, 3);
            cutoff_space = Kokkos::View<double**, MemorySpace>("temp_cutoff", nr_neighbours, buffer_size);
        }

    void exchange_data( aosoa &particle_aosoa, std::vector<int> neighbors, int myrank, int npart){

        auto rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, "rank");
        auto last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, "last_pos");
        auto pos_s = Cabana::slice<core_part_position>(particle_aosoa, "position");
        auto vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, "velocity");
        auto cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, "cutoff");
        auto testvar_s = Cabana::slice<testvar>(particle_aosoa, "testvar");
        auto testvar2_s = Cabana::slice<testvar2>(particle_aosoa, "testvar2");
        auto testvar3_s = Cabana::slice<testvar3>(particle_aosoa, "testvar3");
        auto testvar4_s = Cabana::slice<testvar4>(particle_aosoa, "testvar4");
        int *send_count = (int*) malloc(sizeof(int) * neighbors.size());
        int count_neighbours = 0;
        int end = particle_aosoa.size() - 1;
        for(int i = 0; i < neighbors.size(); i++){
                send_count[i] = 0;
        }

        for(int i = particle_aosoa.size()-1; i>=0; i--){
            if(rank_slice(i) != myrank && rank_slice(i) >= 0){
                int therank = rank_slice(i);
                for(int k = 0; k < neighbors.size(); k++){
                    if(therank == neighbors[k]){
                        therank = k;
                        break;
                    }
                }
            int pos = send_count[therank];
            pos_space(therank, pos, 0) = pos_s(i, 0);
            pos_space(therank, pos, 1) = pos_s(i, 1);
            pos_space(therank, pos, 2) = pos_s(i, 2);
            vel_space(therank, pos, 0) = vel_s(i, 0);
            vel_space(therank, pos, 1) = vel_s(i, 1);
            vel_space(therank, pos, 2) = vel_s(i, 2);
            cutoff_space(therank, pos) = cutoff_s(i);
            for(int ii=0; ii < 4; ii++){
                for(int iii=0; iii < 2; iii++){
                    testvar_space(therank, pos, ii, iii) = testvar_s(i, ii, iii);
                }
            }
            testvar2_space(therank, pos) = testvar2_s(i);
            testvar3_space(therank, pos) = testvar3_s(i);
            testvar4_space(therank, pos) = testvar4_s(i);
            send_count[therank]++;

            while(rank_slice(end) != myrank && end > 0){
                end--;
            }
            if(end > i){
                rank_slice(i) = rank_slice(end);
                pos_s(i, 0) = pos_s(end, 0);
                pos_s(i, 1) = pos_s(end, 1);
                pos_s(i, 2) = pos_s(end, 2);
                vel_s(i, 0) = vel_s(end, 0);
                vel_s(i, 1) = vel_s(end, 1);
                vel_s(i, 2) = vel_s(end, 2);
                cutoff_s(i) = cutoff_s(end);
                for(int ii=0; ii < 4; ii++){
                    for(int iii=0; iii < 2; iii++){
                        testvar_s(i, ii, iii) = testvar_s(end, ii, iii);
                    }
                }
                testvar2_s(i) = testvar2_s(end);
                testvar3_s(i) = testvar3_s(end);
                testvar4_s(i) = testvar4_s(end);
                rank_slice(end) = -1;
            }else{
                rank_slice(i) = -1;
                end++;
            }
            continue;
        }

    }
    // Data collected, need to send information to neighbours to know what to expect
    int *recv_count = (int*) malloc(sizeof(int) * neighbors.size());
    MPI_Request *requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2);
    int req_num = 0;
    for(int i = 0; i < neighbors.size(); i++){
        recv_count[i] = 0;
        if(neighbors[i] == myrank){
            continue;
        }
        MPI_Irecv(&recv_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);
        MPI_Isend(&send_count[i], 1, MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &requests[req_num++]);
    }
    MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    int total_size = 0;
    for(int i = 0; i < neighbors.size(); i++){
             total_size += recv_count[i];
    }
    Kokkos::View<double***, MemorySpace> r_pos_space("temp_pos", neighbors.size(), total_size, 3);
    Kokkos::View<double***, MemorySpace> r_vel_space("temp_vel", neighbors.size(), total_size, 3);
    Kokkos::View<double**, MemorySpace> r_cutoff_space("temp_cutoff", neighbors.size(), total_size);
    Kokkos::View<int****, MemorySpace> r_testvar_space("temp_testvar", neighbors.size(), total_size, 4, 2);
    Kokkos::View<double**, MemorySpace> r_testvar2_space("temp_testvar2", neighbors.size(), total_size);
    Kokkos::View<float**, MemorySpace> r_testvar3_space("temp_testvar3", neighbors.size(), total_size);
    Kokkos::View<int64_t**, MemorySpace> r_testvar4_space("temp_testvar4", neighbors.size(), total_size);

    free(requests);
    requests = (MPI_Request*) malloc(sizeof(MPI_Request) * neighbors.size() * 2 * 7);
    req_num = 0;
    int tag = 0;
    for(int i = 0; i < neighbors.size(); i++){
        if(neighbors[i] != myrank){
            tag = 0;
            MPI_Irecv(&r_pos_space.data()[r_pos_space.extent(1)*r_pos_space.extent(2)*i], recv_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_vel_space.data()[r_vel_space.extent(1)*r_vel_space.extent(2)*i], recv_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_cutoff_space.data()[r_cutoff_space.extent(1)*i], recv_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_testvar_space.data()[r_testvar_space.extent(1) * r_testvar_space.extent(2) * r_testvar_space.extent(3) * i], recv_count[i] * 4 * 2, MPI_INT, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_testvar2_space.data()[r_testvar2_space.extent(1)*i], recv_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_testvar3_space.data()[r_testvar3_space.extent(1)*i], recv_count[i], MPI_FLOAT, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Irecv(&r_testvar4_space.data()[r_testvar4_space.extent(1)*i], recv_count[i], MPI_INT64_T, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            tag = 0;
            MPI_Isend(&pos_space.data()[pos_space.extent(1)*pos_space.extent(2)*i], send_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&vel_space.data()[vel_space.extent(1)*vel_space.extent(2)*i], send_count[i]*3, MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&cutoff_space.data()[cutoff_space.extent(1)*i], send_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&testvar_space.data()[testvar_space.extent(1) * r_testvar_space.extent(2) * r_testvar_space.extent(3) * i], send_count[i] * 4 * 2, MPI_INT, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&testvar2_space.data()[testvar2_space.extent(1)*i], send_count[i], MPI_DOUBLE, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&testvar3_space.data()[testvar3_space.extent(1)*i], send_count[i], MPI_FLOAT, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
            MPI_Isend(&testvar4_space.data()[testvar4_space.extent(1)*i], send_count[i], MPI_INT64_T, neighbors[i], tag++, MPI_COMM_WORLD, &requests[req_num++]);
        }
    }
    MPI_Waitall(req_num, requests, MPI_STATUSES_IGNORE);
    free(requests);
    int recvd = 0;
    int sent = 0;
    for(int i = 0; i < neighbors.size(); i++){
        recvd += recv_count[i];
        sent += send_count[i];
    }
    int size_change = recvd - sent;
    int current_size =  particle_aosoa.size();
    if(size_change != 0){
        particle_aosoa.resize(current_size+size_change);
    }
    auto new_rank_slice = Cabana::slice<neighbour_part_rank>(particle_aosoa, "new_rank");
    for(int i = particle_aosoa.size() - 1; i > end; i--){
        new_rank_slice(i) = -1;
    }
    if(size_change > 0){
        if(sent = 0){
            end = current_size;
        }
        while(end < current_size && end < particle_aosoa.size() && new_rank_slice(end) != -1 ) end++;
        for(int i = 0; i < particle_aosoa.size(); i++){
            new_rank_slice(i) = -1;
        }
    }
    auto new_last_pos_s = Cabana::slice<neighbour_part_old_position>(particle_aosoa, "new_last_pos");
    auto new_pos_s = Cabana::slice<core_part_position>(particle_aosoa, "new_position");
    auto new_vel_s = Cabana::slice<core_part_velocity>(particle_aosoa, "new_velocity");
    auto new_cutoff_s = Cabana::slice<neighbour_part_cutoff>(particle_aosoa, "new_cutoff");
    auto new_testvar_s = Cabana::slice<testvar>(particle_aosoa, "new_testvar");
    auto new_testvar2_s = Cabana::slice<testvar2>(particle_aosoa, "new_testvar2");
    auto new_testvar3_s = Cabana::slice<testvar3>(particle_aosoa, "new_testvar3");
    auto new_testvar4_s = Cabana::slice<testvar4>(particle_aosoa, "new_testvar4");
    int x = 0;
    for(int j = 0; j < neighbors.size(); j++){
        for(int i = 0; i < recv_count[j]; i++){
            new_pos_s(end+x, 0) = r_pos_space(j,i,0);
            new_pos_s(end+x, 1) = r_pos_space(j,i,1);
            new_pos_s(end+x, 2) = r_pos_space(j,i,2);
            new_vel_s(end+x, 0) = r_vel_space(j,i,0);
            new_vel_s(end+x, 1) = r_vel_space(j,i,1);
            new_vel_s(end+x, 2) = r_vel_space(j,i,2);
            new_cutoff_s(end+x) = r_cutoff_space(j,i);
            for(int ii=0; ii < 4;ii++){
                for(int iii=0; iii < 2;iii++){
                    new_testvar_s(end+x, ii, iii) = r_testvar_space(j,i, ii, iii);
                }
            }
            new_testvar2_s(end+x) = r_testvar2_space(j,i);
            new_testvar3_s(end+x) = r_testvar3_space(j,i);
            new_testvar4_s(end+x) = r_testvar4_space(j,i);
            x++;
        }
    }
    free(recv_count);
    free(send_count);
}};

KOKKOS_INLINE_FUNCTION
int get_oneD_rank(int x_r, int y_r, int z_r, int x_ranks, int y_ranks, int z_ranks){
    int oneD_rank = z_r*(x_ranks*y_ranks) + y_r*(x_ranks) + x_r;
    return oneD_rank;
}
KOKKOS_INLINE_FUNCTION
void get_threeD_rank(int rank, int *x, int *y, int *z, int x_ranks, int y_ranks, int z_ranks){
    int z_r = rank / (x_ranks*y_ranks);
    int y_r = (rank - z_r*x_ranks*y_ranks) / x_ranks;
    int x_r = rank - z_r*x_ranks*y_ranks - y_r*x_ranks;
    *x = x_r;
    *y = y_r;
    *z = z_r;
};

template<class PartPosSlice, class RankSlice>
struct _rank_update_functor{
    boundary _box;
    PartPosSlice _part_pos;
    RankSlice _rank;
    int _myrank;
    int _xranks;
    int _yranks;
    int _zranks;
    int _local_x_rank;
    int _local_y_rank;
    int _local_z_rank;
    int _nranks;

    KOKKOS_INLINE_FUNCTION
    _rank_update_functor(boundary box, PartPosSlice pos, RankSlice rank, int xranks, int yranks, int zranks, int myrank):
        _box(box), _part_pos(pos), _rank(rank), _xranks(xranks), _yranks(yranks), _zranks(zranks), _myrank(myrank){
        get_threeD_rank(myrank, &_local_x_rank, &_local_y_rank, &_local_z_rank, xranks, yranks, zranks);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ix, const int ij) const{
        int xr, yr, zr;
        get_threeD_rank(_myrank, &xr, &yr, &zr, _xranks, _yranks, _zranks);
        xr = _local_x_rank;
        yr = _local_y_rank;
        zr = _local_z_rank;
        if(_part_pos.access(ix, ij, 0) >= _box.local_x_max){
            xr = xr + 1;
            if( xr >= _xranks ) xr = 0;
        }
        if(_part_pos.access(ix, ij, 0) < _box.local_x_min){
            xr = xr - 1;
            if( xr < 0 ) xr = _xranks-1;
        }
        if(_part_pos.access(ix, ij, 1) >= _box.local_y_max){
            yr = yr + 1;
            if( yr >= _yranks ) yr = 0;
        }
        if(_part_pos.access(ix, ij, 1) < _box.local_y_min){
            yr = yr - 1;
            if( yr < 0 ) yr = _yranks-1;
        }
        if(_part_pos.access(ix, ij, 2) >= _box.local_z_max){
            zr = zr + 1;
            if( zr >= _zranks ) zr = 0;
        }
        if(_part_pos.access(ix, ij, 2) < _box.local_z_min){
            zr = zr - 1;
            if( zr < 0 ) zr = _zranks-1;
        }
        _rank.access(ix, ij) = get_oneD_rank(xr, yr, zr, _xranks, _yranks, _zranks);
    }};

'''
    assert out == correct
    part.add_element("testvar5", "int32_t")
    with pytest.raises(NotImplementedError) as excinfo:
        out = backend.mpi_headers(config, part)
    assert ("Don't know currently how to support element testvar5 with "
            "datatype int32_t for MPI" in str(excinfo.value))

def test_cabana_pir_get_current_kernel():
    backend = Cabana_PIR()
    assert backend.get_current_kernel() is None

def test_cabana_pir_get_extra_symbols():
    a = Cabana_PIR()
    x = StructureType()
    a.add_structure(x, "x")
    a.add_coupler(coupler_test2())
    result = a.get_extra_symbols([])
    assert result[0][0] == "x"
    assert result[0][1].datatype == x
    assert result[1] == "x"

def test_cabana_pir_get_mpi_comm_after_bcs():
    a = Cabana_PIR()
    
    part = Particle()
    part.add_element("test", "int")
    a._particle = part
    out = a.gen_mpi_comm_after_bcs()
    correct = '''    Cabana::deep_copy(particle_aosoa_host, particle_aosoa);
    _migrator.exchange_data(particle_aosoa_host, neighbors, myrank, particle_aosoa_host.size());
    particle_aosoa.resize(particle_aosoa_host.size());
    Cabana::deep_copy(particle_aosoa, particle_aosoa_host);
    core_part_position_slice = Cabana::slice<core_part_position>(particle_aosoa);
    core_part_velocity_slice = Cabana::slice<core_part_velocity>(particle_aosoa);
    neighbour_part_cutoff_slice = Cabana::slice<neighbour_part_cutoff>(particle_aosoa);
    neighbour_part_deletion_flag_slice = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa);
    test_slice = Cabana::slice<test>(particle_aosoa);

    simd_policy = Cabana::SimdPolicy<VectorLength, ExecutionSpace>(0, particle_aosoa.size());

'''
    assert out == correct

def test_cabana_pir_gen_makefile():
    a = Cabana_PIR()
    with pytest.raises(NotImplementedError):
        a._gen_makefile()
