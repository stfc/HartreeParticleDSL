from __future__ import annotations

from HartreeParticleDSL.HartreeParticleDSL import get_mpi

from HartreeParticleDSL.HartreeParticleDSLExceptions import IllegalLoopError, UnsupportedCodeError

from HartreeParticleDSL.Particle_IR.nodes.array_mixin import ArrayMixin
from HartreeParticleDSL.Particle_IR.nodes.assignment import Assignment
from HartreeParticleDSL.Particle_IR.nodes.auto_reference import AutoReference
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.call import Call
from HartreeParticleDSL.Particle_IR.nodes.funcdef import FuncDef
from HartreeParticleDSL.Particle_IR.nodes.invoke import Invoke
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.kernels import PerPartKernel, \
    PairwiseKernel, SourceBoundaryKernel, SinkBoundaryKernel
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal
from psyclone.psyir.nodes import Node
from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.nodes.statement import Return, EmptyStatement
from HartreeParticleDSL.Particle_IR.nodes.operation import BinaryOperation, UnaryOperation
from HartreeParticleDSL.Particle_IR.nodes.particle_position_reference import \
        ParticlePositionReference
from HartreeParticleDSL.Particle_IR.nodes.particle_reference import ParticleReference
from HartreeParticleDSL.Particle_IR.nodes.pointer_reference import PointerReference
from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol
from HartreeParticleDSL.backends.base_backend.pir_visitor import PIR_Visitor

from HartreeParticleDSL.Particle_IR.datatypes.datatype import type_mapping_str, ScalarType, \
        INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE, INT64_TYPE, INT32_TYPE, BOOL_TYPE, STRING_TYPE, \
        BASE_PARTICLE_TYPE, StructureType, ArrayType

class Cabana_PIR_Visitor(PIR_Visitor):
    '''
    Generic Cabana PIR visitor.

    :param str indent_string: The indentation to use in this Visitor. Default \
                              "    ".
    :param int initial_indent_depth: Specifies the initial indentation, default \
                                     is 0.

    '''

    def __init__(self, parent, **kwargs):
        super().__init__(**kwargs)
        self._in_kernel = False
        self._slices = []
        self._parent = parent

    @classmethod
    def get_cpp_datatype(cls, datatype: DataType):
        type_map = {INT_TYPE: "int",
                FLOAT_TYPE: "float",
                DOUBLE_TYPE: "double",
                INT64_TYPE: "int64_t",
                INT32_TYPE: "int32_t",
                BOOL_TYPE: "bool",
                STRING_TYPE: "char*",
                BASE_PARTICLE_TYPE: "part"}
        t = type_map.get(datatype)
        if isinstance(datatype, ArrayType):
            childtype = datatype._datatype
            t = type_map.get(childtype)
        if t is not None:
            return t
        for key in type_mapping_str.keys():
            if isinstance(datatype, StructureType):
                if type_mapping_str[key] == datatype:
                    return "struct " + key
            if type_mapping_str[key] == datatype:
                return key

    def addSlice(self, slice_name: str):
        if slice_name not in self._slices:
            self._slices.append(slice_name)

    def visit_autosymbol_node(self, symbol: AutoSymbol) -> str:
        return f"auto {symbol.name} = {symbol.initial_value}"

    def visit_symbol_node(self, symbol: Symbol) -> str:
        # Find the name in the type_mapping_str
        type_str = ""
        for key in type_mapping_str.keys():
            if type_mapping_str[key] == symbol.datatype:
                type_str = key
        return Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)

    def visit_pointersymbol_node(self, symbol: PointerSymbol) -> str:
        # Find the name in the type_mapping_str
        type_str = ""
        type_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype._datatype)
        return f"{type_str}*"

    def visit_arraysymbol_node(self, symbol: ArraySymbol) -> str:
        # Find the name in the type_mapping_str
        type_str = Cabana_PIR_Visitor.get_cpp_datatype(symbol.datatype)
        all_dynamic = True
        for index, extent in enumerate(symbol.datatype.shape):
            if isinstance(extent, int) and (index == 0 or all_dynamic == False):
                all_dynamic = False
            elif isinstance(extent, int) or not all_dynamic:
                raise UnsupportedCodeError("Got an ArraySymbol with a mixed "
                                           "dynamic and fixed typing, which is "
                                           "not supported.")
        if all_dynamic:
            type_str = type_str + ("*"*len(symbol.datatype.shape))
        else:
            for extent in symbol.datatype.shape:
                type_str = type_str + "[" + f"{extent}" + "]"
        return type_str

    # Implement visitors
    def visit_pairwisekernel_node(self, node: PairwiseKernel) -> str:
        self._in_kernel = True
        self._slices = []
        raise NotImplementedError("Pairwise kernels not yet supported in cabana.")
        # self._parent.register_kernel(node.name, node)

    def visit_perpartkernel_node(self, node: PerPartKernel) -> str:
        self._in_kernel = True
        self._slices = []
        # Check rules. Doesn't contain invokes.
        if len(node.walk(Invoke)) != 0:
            raise UnsupportedCodeError("Per part kernel cannot contain Invokes.")

        # Need to find any random number calls
        calls = node.walk(Call)
        first_random_call = None
        last_random_call = None
        for call in calls:
            if call.func_name == "random_number":
                if first_random_call == None:
                    first_random_call = call
                last_random_call = call
        # If we have a random call then we wrap in the wrapper calls
        # The auto symbol handles the initialisation of the random state so
        # we just add something after the last.
        if first_random_call is not None:
            generator_symbol = node.symbol_table.lookup("_generator")
            post_assign = Call.create("_random_pool.free_state", [AutoReference(generator_symbol)])
            assign_ancestor = last_random_call.ancestor(Assignment)
            if assign_ancestor is not None:
                container = assign_ancestor.parent
                position = assign_ancestor.position
                container.addchild(post_assign, position+1)
            else:
                container = last_random_call.parent
                position = last_random_call.position
                container.addchild(post_assign, position+1)

        # Arg 1 is always a particle
        #TODO Save the particle 1 symbol name.
        self._parent.register_kernel(node.name, node)
        argument_names = []
        for arg in node.arguments:
            # Need the correct typing on these arguments
            argument_names.append(arg.symbol.name)

        self.indent()
        extra_arrays = self._parent.get_writable_arrays()
        # Extra indentation for symbol table
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in (argument_names + list(extra_arrays.keys())):
                if not isinstance(symbol, AutoSymbol):
                    sym = self._visit(symbol)
                    symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
                else:
                    symbol_list = symbol_list + self._nindent + self._visit(symbol) + ";\n"
        self.dedent()

        #Current indentation level isn't quite correct
        kernel_body = self._visit(node.body)
        self.dedent()

        # Start codegen.

        # Create the templated functor
        rval = self._nindent
        if len(self._slices) > 0:
            rval = rval + "template < "
            all_slices = []
            for slices in self._slices:
                all_slices.append("class " + slices.upper())
            classes = ", ".join(all_slices)
            rval = rval + classes + " >\n"
        rval = rval + f"{self._nindent}struct {node.name}_functor" + "{\n"
        self.indent()
        rval = rval + f"{self._nindent}config_struct_type _{node.arguments[1].symbol.name};\n"
        if self._parent.get_require_random():
            rval = rval + f"{self._nindent}Kokkos::Random_XorShift64_Pool<> _random_pool;\n"
        all_slices = []
        for slices in self._slices:
            rval = rval + f"{self._nindent}{slices.upper()} _{slices};\n"

        for name in extra_arrays.keys():
            rval = rval + f"{self._nindent}Kokkos::View<{extra_arrays[name][0]}, MemorySpace> {name};\n"

        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            rval = rval + f"{self._nindent}{typename} {structure};\n"

        # Constructor
        rval = rval + "\n"
        rval = rval + f"{self._nindent}KOKKOS_INLINE_FUNCTION\n"
        rval = rval + f"{self._nindent} {node.name}_functor( "
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"{slices.upper()} {slices}")
        all_slices.append("config_struct_type " + node.arguments[1].symbol.name)
        if self._parent.get_require_random():
            all_slices.append("Kokkos::Random_XorShift64_Pool<> random_pool")
        for name in extra_arrays.keys():
            all_slices.append(f"Kokkos::View<{extra_arrays[name][0]}, MemorySpace> _{name}")
        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            all_slices.append(f"{typename} {structure.upper()}")
        classes = ", ".join(all_slices)
        rval = rval + classes + "):\n"
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"_{slices}({slices})")
        for structure in self._parent.structures:
            all_slices.append(f"{structure}({structure.upper()})")
        if self._parent.get_require_random():
            all_slices.append("_random_pool(random_pool)")
        for name in extra_arrays.keys():
            all_slices.append(f"{name}(_{name})")
        classes = ", ".join(all_slices)
        rval = rval + f"{self._nindent}{classes}"
        if len(all_slices) > 0:
            rval = rval + ", "
        rval = rval + f"_{node.arguments[1].symbol.name}({node.arguments[1].symbol.name})" + "{}\n"

        # Need to have an update_structs call if there are structures
        if len(self._parent.structures) > 0:
            rval = rval + f"\n{self._nindent}void update_structs("
            all_structs = []
            for struct in self._parent.structures:
                if self._parent._structures[structure] in type_mapping_str.values():
                    typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
                else:
                    typename = structure
                all_structs.append(typename + " " + struct.upper())
            rval = rval + ", ".join(all_structs) + "){\n"
            self.indent()
            for struct in self._parent.structures:
                rval = rval + self._nindent + struct + " = " + struct.upper() + ";\n"
            self.dedent()
            rval = rval + self._nindent + "}\n"

        rval = rval + "\n"
        rval = rval + self._nindent + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + self._nindent + "void operator()(const int i, const int a) const{\n"
        rval = rval + symbol_list
        rval = rval + kernel_body
        rval = rval + self._nindent + "}\n"

        self.dedent()
        rval = rval + self._nindent + "};\n"
        for slices in self._slices:
            self._parent.add_kernel_slices(node.name, slices)
        return rval

    def visit_sourceboundarykernel_node(self, node: SourceBoundaryKernel) -> str:
        self._in_kernel = True

        self._slices = []
        # Check rules. Doesn't contain invokes.
        if len(node.walk(Invoke)) != 0:
            raise UnsupportedCodeError("Source boundary kernel cannot contain Invokes.")

        # Need to find any random number calls
        calls = node.walk(Call)
        first_random_call = None
        last_random_call = None
        for call in calls:
            if call.func_name == "random_number":
                if first_random_call == None:
                    first_random_call = call
                last_random_call = call
        # If we have a random call then we wrap in the wrapper calls
        # The auto symbol handles the initialisation of the random state so
        # we just add something after the last.
        if first_random_call is not None:
            generator_symbol = node.symbol_table.lookup("_generator")
            post_assign = Call.create("_random_pool.free_state", [AutoReference(generator_symbol)])
            assign_ancestor = last_random_call.ancestor(Assignment)
            if assign_ancestor is not None:
                container = assign_ancestor.parent
                position = assign_ancestor.position
                container.addchild(post_assign, position+1)
            else:
                container = last_random_call.parent
                position = last_random_call.position
                container.addchild(post_assign, position+1)

        # Arg 1 is always a particle
        #TODO Save the particle 1 symbol name.
        self._parent.register_kernel(node.name, node)
        argument_names = []
        for arg in node.arguments:
            # Need the correct typing on these arguments
            argument_names.append(arg.symbol.name)

        extra_arrays = self._parent.get_writable_arrays()
        self.indent()
        # Extra indentation for symbol table
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in (argument_names + list(extra_arrays.keys())):
                if not isinstance(symbol, AutoSymbol):
                    sym = self._visit(symbol)
                    symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
                else:
                    symbol_list = symbol_list + self._nindent + self._visit(symbol) + ";\n"
        self.dedent()

        #Current indentation level isn't quite correct
        kernel_body = self._visit(node.body)
        self.dedent()

        # Start codegen.

        # Create the templated functor
        rval = self._nindent
        if len(self._slices) > 0:
            rval = rval + "template < "
            all_slices = []
            for slices in self._slices:
                all_slices.append("class " + slices.upper())
            classes = ", ".join(all_slices)
            rval = rval + classes + " >\n"
        rval = rval + f"{self._nindent}struct {node.name}_functor" + "{\n"
        self.indent()
        rval = rval + f"{self._nindent}config_struct_type _{node.arguments[1].symbol.name};\n"
        if self._parent.get_require_random():
            rval = rval + f"{self._nindent}Kokkos::Random_XorShift64_Pool<> _random_pool;\n"
        for name in extra_arrays.keys():
            rval = rval + f"{self._nindent}Kokkos::View<{extra_arrays[name][0]}, MemorySpace> {name};\n"
        all_slices = []
        for slices in self._slices:
            rval = rval + f"{self._nindent}{slices.upper()} _{slices};\n"
        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            rval = rval + f"{self._nindent}{typename} {structure};\n"

        # Constructor
        rval = rval + "\n"
        rval = rval + f"{self._nindent}KOKKOS_INLINE_FUNCTION\n"
        rval = rval + f"{self._nindent} {node.name}_functor( "
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"{slices.upper()} {slices}")
        all_slices.append("config_struct_type " + node.arguments[1].symbol.name)
        if self._parent.get_require_random():
            all_slices.append("Kokkos::Random_XorShift64_Pool<> random_pool")
        for name in extra_arrays.keys():
            all_slices.append(f"Kokkos::View<{extra_arrays[name][0]}, MemorySpace> _{name}")
        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            all_slices.append(f"{typename} {structure.upper()}")
        classes = ", ".join(all_slices)
        rval = rval + classes + "):\n"
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"_{slices}({slices})")
        for structure in self._parent.structures:
            all_slices.append(f"{structure}({structure.upper()})")
        if self._parent.get_require_random():
            all_slices.append("_random_pool(random_pool)")
        for name in extra_arrays.keys():
            all_slices.append(f"{name}(_{name})")
        classes = ", ".join(all_slices)
        rval = rval + f"{self._nindent}{classes}"
        if len(all_slices) > 0:
            rval = rval + ", "
        rval = rval + f"_{node.arguments[1].symbol.name}({node.arguments[1].symbol.name})" + "{}\n"

        # Need to have an update_structs call if there are structures
        if len(self._parent.structures) > 0:
            rval = rval + f"\n{self._nindent}void update_structs("
            all_structs = []
            for struct in self._parent.structures:
                if self._parent._structures[structure] in type_mapping_str.values():
                    typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
                else:
                    typename = structure
                all_structs.append(typename + " " + struct.upper())
            rval = rval + ", ".join(all_structs) + "){\n"
            self.indent()
            for struct in self._parent.structures:
                rval = rval + self._nindent + struct + " = " + struct.upper() + ";\n"
            self.dedent()
            rval = rval + self._nindent + "}\n"

        # Generate inflow count function
        rval = rval + f"\n{self._nindent}int get_inflow_count()" + "{\n"
        self.indent()
        rval = rval + f"{self._nindent}return {node.source_count};\n"
        self.dedent()
        rval = rval + f"{self._nindent}" + "}\n"


        rval = rval + "\n"
        rval = rval + self._nindent + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + self._nindent + "void operator()(const int i, const int a) const{\n"
        rval = rval + symbol_list
        rval = rval + kernel_body
        rval = rval + self._nindent + "}\n"

        self.dedent()
        rval = rval + self._nindent + "};\n"
        for slices in self._slices:
            self._parent.add_kernel_slices(node.name, slices)

        self._in_kernel = False
        return rval

    def visit_sinkboundarykernel_node(self, node: SinkBoundaryKernel) -> str:
        self._in_kernel = True

        self._slices = []
        # Check rules. Doesn't contain invokes.
        if len(node.walk(Invoke)) != 0:
            raise UnsupportedCodeError("Sink boundary kernel cannot contain Invokes.")

        # Need to find any random number calls
        calls = node.walk(Call)
        first_random_call = None
        last_random_call = None
        for call in calls:
            if call.func_name == "random_number":
                if first_random_call == None:
                    first_random_call = call
                last_random_call = call
        # If we have a random call then we wrap in the wrapper calls
        # The auto symbol handles the initialisation of the random state so
        # we just add something after the last.
        if first_random_call is not None:
            generator_symbol = node.symbol_table.lookup("_generator")
            post_assign = Call.create("_random_pool.free_state", [AutoReference(generator_symbol)])
            assign_ancestor = last_random_call.ancestor(Assignment)
            if assign_ancestor is not None:
                container = assign_ancestor.parent
                position = assign_ancestor.position
                container.addchild(post_assign, position+1)
            else:
                container = last_random_call.parent
                position = last_random_call.position
                container.addchild(post_assign, position+1)

        # Arg 1 is always a particle
        #TODO Save the particle 1 symbol name.
        self._parent.register_kernel(node.name, node)
        argument_names = []
        for arg in node.arguments:
            # Need the correct typing on these arguments
            argument_names.append(arg.symbol.name)

        extra_arrays = self._parent.get_writable_arrays()
        self.indent()
        # Extra indentation for symbol table
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in (argument_names + list(extra_arrays.keys())):
                if not isinstance(symbol, AutoSymbol):
                    sym = self._visit(symbol)
                    symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
                else:
                    symbol_list = symbol_list + self._nindent + self._visit(symbol) + ";\n"
        self.dedent()

        #Current indentation level isn't quite correct
        kernel_body = self._visit(node.body)
        self.dedent()

        # Start codegen.

        # Create the templated functor
        rval = self._nindent
        if len(self._slices) > 0:
            rval = rval + "template < "
            all_slices = []
            for slices in self._slices:
                all_slices.append("class " + slices.upper())
            classes = ", ".join(all_slices)
            rval = rval + classes + " >\n"
        rval = rval + f"{self._nindent}struct {node.name}_functor" + "{\n"
        self.indent()
        rval = rval + f"{self._nindent}config_struct_type _{node.arguments[1].symbol.name};\n"
        if self._parent.get_require_random():
            rval = rval + f"{self._nindent}Kokkos::Random_XorShift64_Pool<> _random_pool;\n"
        for name in extra_arrays.keys():
            rval = rval + f"{self._nindent}Kokkos::View<{extra_arrays[name][0]}, MemorySpace> {name};\n"
        all_slices = []
        for slices in self._slices:
            rval = rval + f"{self._nindent}{slices.upper()} _{slices};\n"
        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            rval = rval + f"{self._nindent}{typename} {structure};\n"

        # Constructor
        rval = rval + "\n"
        rval = rval + f"{self._nindent}KOKKOS_INLINE_FUNCTION\n"
        rval = rval + f"{self._nindent} {node.name}_functor( "
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"{slices.upper()} {slices}")
        all_slices.append("config_struct_type " + node.arguments[1].symbol.name)
        if self._parent.get_require_random():
            all_slices.append("Kokkos::Random_XorShift64_Pool<> random_pool")
        for name in extra_arrays.keys():
            all_slices.append(f"Kokkos::View<{extra_arrays[name][0]}, MemorySpace> _{name}")
        for structure in self._parent.structures:
            if self._parent._structures[structure] in type_mapping_str.values():
                typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
            else:
                typename = structure
            all_slices.append(f"{typename} {structure.upper()}")
        classes = ", ".join(all_slices)
        rval = rval + classes + "):\n"
        all_slices = []
        for slices in self._slices:
            all_slices.append(f"_{slices}({slices})")
        for structure in self._parent.structures:
            all_slices.append(f"{structure}({structure.upper()})")
        if self._parent.get_require_random():
            all_slices.append("_random_pool(random_pool)")
        for name in extra_arrays.keys():
            all_slices.append(f"{name}(_{name})")
        classes = ", ".join(all_slices)
        rval = rval + f"{self._nindent}{classes}"
        if len(all_slices) > 0:
            rval = rval + ", "
        rval = rval + f"_{node.arguments[1].symbol.name}({node.arguments[1].symbol.name})" + "{}\n"

        # Need to have an update_structs call if there are structures
        if len(self._parent.structures) > 0:
            rval = rval + f"\n{self._nindent}void update_structs("
            all_structs = []
            for struct in self._parent.structures:
                if self._parent._structures[structure] in type_mapping_str.values():
                    typename = [k for k, v in type_mapping_str.items() if v == self._parent._structures[structure]][0]
                else:
                    typename = structure
                all_structs.append(typename + " " + struct.upper())
            rval = rval + ", ".join(all_structs) + "){\n"
            self.indent()
            for struct in self._parent.structures:
                rval = rval + self._nindent + struct + " = " + struct.upper() + ";\n"
            self.dedent()
            rval = rval + self._nindent + "}\n"

        rval = rval + "\n"
        rval = rval + self._nindent + "KOKKOS_INLINE_FUNCTION\n"
        rval = rval + self._nindent + "void operator()(const int i, const int a) const{\n"
        rval = rval + symbol_list
        rval = rval + kernel_body
        rval = rval + self._nindent + "}\n"

        self.dedent()
        rval = rval + self._nindent + "};\n"
        for slices in self._slices:
            self._parent.add_kernel_slices(node.name, slices)

        self._in_kernel = False
        return rval

    def visit_mainkernel_node(self, node: MainKernel) -> str:
        self._in_kernel = False

        #Current indentation level isn't quite correct
        self.indent()
        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if symbol.name == "config":
                continue
            sym = self._visit(symbol)
            symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"

        self.dedent()
        body = self._visit(node.body)
        rval = self._nindent + f"int {node.name}( int argc, char* argv[] )" + "{\n"
        rval = rval + symbol_list
        # Append the code body
        rval = rval + body
        rval = rval + self._nindent + "}\n"

        return rval

    def visit_body_node(self, node: Body) -> str:
        self.indent()
        string = ""
        for child in node.children:
            if not isinstance(child, EmptyStatement):
                string = string + self._nindent + self._visit(child) + "\n"
        self.dedent()
        return string

    def visit_ifelseblock_node(self, node: IfElseBlock) -> str:
        cond = self._visit(node.condition)
        ifbody = self._visit(node.ifbody)
#        if node.is_else_if():
#            raise NotImplementedError()
        string = f"if({cond})" + "{\n"
        string = string + ifbody
        string = string + self._nindent + "}"
        if len(node.elsebody.children) > 0:
            string = string + "else"
            if node.is_else_if():
                self.dedent()
                child = self._visit(node.elsebody)
                child = child.lstrip()
                child = child.rstrip()
                self.indent()
                string = string + " " + child
            else:
                string = string + "{\n"
                string = string + self._visit(node.elsebody)
                string = string + self._nindent + "}"
        return string

    def visit_invoke_node(self, node: Invoke) -> str:
        rval = ""
        for invoke in node.children:
            # TODO Check if the call updates particle position.
            name = invoke.value
            pir_kernel = self._parent._kernels[name]
            assigns = pir_kernel.walk(Assignment)
            updates_part_pos = False
            for assign in assigns:
                if isinstance(assign.lhs, ParticlePositionReference):
                    updates_part_pos = True
                # TODO If it contains a "Call" node to an unknown call maybe
                # this should just be set to True as well?

            rval = rval + f"Kokkos::deep_copy(config.config, config.config_host);\n"
            if len(self._parent.structures) > 0:
                rval = rval + self._nindent
                rval = rval + f"{invoke.value}.update_structs("
                struct_list = []
                for struct in self._parent.structures:
                    struct_list.append(struct)
                rval = rval + ", ".join(struct_list) + ");\n"
            # TODO Need to do stuff depending on kernel type.
            if isinstance(pir_kernel, PerPartKernel):
                rval = rval + f"{self._nindent}Cabana::simd_parallel_for(simd_policy, {invoke.value}, "
                rval = rval + "\"" + invoke.value + "\");\n"
                # TODO Check if we need to block (i.e. only if kernel dependencies or final kernel)
                rval = rval + self._nindent + "Kokkos::fence();"
            elif isinstance(pir_kernel, SourceBoundaryKernel):
                if get_mpi():
                    rval = rval + f"{self._nindent}if(myrank == 0)" + "{\n"
                else:
                    rval = rval + f"{self._nindent}" + "{\n"
                self.indent()
                rval = rval + f"{self._nindent}int new_parts = {invoke.value}.get_inflow_count();\n"
                rval = rval + f"{self._nindent}int old_size = particle_aosoa.size();\n"
                rval = rval + f"{self._nindent}int new_size = old_size + new_parts;\n"
                # We need to resize both particle structures
                rval = rval + f"{self._nindent}particle_aosoa.resize(new_size);\n"
                rval = rval + f"{self._nindent}particle_aosoa_host.resize(new_size);\n"
                # Need to redo slices and functors
                current_indent = len(self._nindent)
                indent = len(self._indent)
                rval = rval + self._parent.gen_slices_and_functors(current_indent=current_indent, indent=indent)
                # Loop over the new stuff and do the source bc
                rval = rval + f"{self._nindent}Cabana::SimdPolicy<VectorLength, ExecutionSpace> simd(old_size+1, new_size);\n"
                rval = rval + f"{self._nindent}Cabana::simd_parallel_for(simd, {invoke.value}, "
                rval = rval + "\"" + invoke.value + "\");\n"
                self.dedent()
                rval = rval + f"{self._nindent}" + "}\n"
            elif isinstance(pir_kernel, SinkBoundaryKernel):
                # Sink boundaries update the neighbour_part_deletion_flag
                rval = rval + f"{self._nindent}" + "{\n"
                self.indent()
                rval = rval + f"{self._nindent}Cabana::simd_parallel_for(simd_policy, {invoke.value}, "
                rval = rval + "\"" + invoke.value + "\");\n"
                # Positive values are deleted
                rval = rval + f"{self._nindent}auto sort_data = Cabana::binByKey(neighbour_part_deletion_flag_slice, 2);\n"
                rval = rval + f"{self._nindent}Cabana::permute(sort_data, particle_aosoa);\n"
                # TODO Resize the array.
                rval = rval + f"{self._nindent}Cabana::deep_copy(particle_aosoa_host, particle_aosoa);\n"
                rval = rval + f"{self._nindent}auto local_deletion_flags = Cabana::slice<neighbour_part_deletion_flag>(particle_aosoa_host);\n"
                rval = rval + f"{self._nindent}int end = particle_aosoa.size();\n"
                rval = rval + f"{self._nindent}for(int j = particle_aosoa.size()-1; j > 0; j--)" + "{\n"
                self.indent()
                rval = rval + f"{self._nindent}if(local_deletion_flags.size() <= 0)" + "{\n"
                self.indent()
                rval = rval + f"{self._nindent}end = j+1;\n"
                rval = rval + f"{self._nindent}break;\n"
                self.dedent()
                rval = rval + f"{self._nindent}" + "}\n"
                self.dedent()
                rval = rval + f"{self._nindent}" + "}\n"
                rval = rval + f"{self._nindent}particle_aosoa.resize(end);\n"
                rval = rval + f"{self._nindent}particle_aosoa_host.resize(end);\n"
                # Redo slices and functors
                current_indent = len(self._nindent)
                indent = len(self._indent)
                rval = rval + self._parent.gen_slices_and_functors(current_indent=current_indent, indent=indent)
                self.dedent()
                rval = rval + f"{self._nindent}" + "}\n"
            else:
                raise NotImplementedError()
            # If we have MPI and move particles do the pre-boundary condition stuff.
            if updates_part_pos and get_mpi():
                rval = rval + "\n" + self._parent.gen_mpi_comm_before_bcs(len(self._nindent), len(self._indent))


            if updates_part_pos and (self._parent.boundary_condition is not None):
                bound = self._parent.boundary_condition
                if len(self._parent.structures) > 0:
                    rval = rval + "\n" + self._nindent
                    rval = rval + f"{invoke.value}.update_structs("
                    struct_list = []
                    for struct in self._parent.structures:
                        struct_list.append(struct)
                    rval = rval + ", ".join(struct_list) + ");"
                rval = rval + f"\n{self._nindent}Cabana::simd_parallel_for(simd_policy, {bound.name}, "
                rval = rval + "\"" + bound.name + "\");\n"
                rval = rval + self._nindent + "Kokkos::fence();"
            # If we have MPI and move particles do the post-boundary condition stuff
            if updates_part_pos and get_mpi():
                rval = rval + "\n" + self._parent.gen_mpi_comm_after_bcs(len(self._nindent), len(self._indent))
        return rval

    def visit_arraymember_node(self, node: ArrayMember) -> str:
        indices = []
        for index in node.indices:
            indices.append(f"[{self._visit(index)}]")
        return node.name + "".join(indices)

    def visit_arrayreference_node(self, node: ArrayReference) -> str:
        indices = []
        extra_arrays = list(self._parent.get_writable_arrays().keys())
        if node.symbol.name in extra_arrays:
            for index in node.indices:
                indices.append(f"{self._visit(index)}")
            return node.symbol.name + "(" + ", ".join(indices) + ")"
        for index in node.indices:
            indices.append(f"[{self._visit(index)}]")
        return node.symbol.name + "".join(indices)

    def visit_assignment_node(self, node: Assignment) -> str:
        lhs = self._visit(node.lhs)
        rhs = self._visit(node.rhs)
        return f"{lhs} = {rhs};"

    def visit_call_node(self, node: Call) -> str:
        # Need to handle special function calls
        # Use the call in the parent
        func_name = node.func_name
        args = []
        for arg in node.children:
            args.append(self._visit(arg))
        try:
            current_indent = len(self._nindent)
            indent = len(self._indent)
            call_str = self._parent.call_language_function(func_name, *args,
                    current_indent=current_indent, indent=indent).lstrip()
            end = ""
            if isinstance(node.parent, Body) and func_name != "cleanup":
                end = ";\n"
            return call_str + end
        except AttributeError as err:
            pass
        arg_string = ", ".join(args)
        end = ""
        if isinstance(node.parent, Body):
            end = ";"
        return f"{func_name}({arg_string}){end}"

    def visit_funcdef_node(self, node: FuncDef) -> str:
        # Don't expect this to appear much in real code, however
        # its used extensively in testing the rest of this code
        # so the implementation here exists.

        return_type = "void"
        returns = node.walk(Return)
        if len(returns) != 0:
            if len(returns[0].children) != 0:
                if isinstance(returns[0].children[0], Literal):
                    datatype = returns[0].children[0].datatype
                    if datatype.intrinsic == ScalarType.Intrinsic.INTEGER:
                        return_type = "int"
                    elif datatype.intrinsic == ScalarType.Intrinsic.FLOAT:
                        return_type = "double"
                    elif datatype.intrinsic == ScalarType.Intrinsic.BOOLEAN:
                        return_type = "bool"
                    elif datatype.intrinsic == ScalarType.Intrinsic.CHARACTER:
                        return_type = "char*"
                elif isinstance(returns[0].children[0], Reference):
                        return_type = Cabana_PIR_Visitor.get_cpp_datatype(returns[0].children[0].symbol.datatype)
                else:
                    raise NotImplementedError()

        body = self._visit(node.body)

        argument_names = []
        arg_base_names = []
        for arg in node.arguments:
            name = arg.symbol.name
            arg_base_names.append(name)
            typ = Cabana_PIR_Visitor.get_cpp_datatype(arg.symbol.datatype)
            if isinstance(arg, PointerReference):
                argument_names.append(f"*{typ} {name}")
            else:
                argument_names.append(f"{typ} {name}")
        extra_array_names = []
        for name in self._parent.get_writable_arrays().keys():
            extra_array_names.append(name)

        symbols = node.symbol_table.get_symbols()
        sym_strings = []
        symbol_list = ""
        self.indent()
        for pairs in symbols.keys():
            symbol = symbols[pairs]
            if pairs not in (arg_base_names + extra_array_names):
                sym = self._visit(symbol)
                symbol_list = symbol_list + self._nindent + sym + " " + pairs + ";\n"
        self.dedent()

        rval = f"{self._nindent}{return_type} {node.name}("
        rval = rval + ", ".join(argument_names) + "){\n"
        rval = rval + symbol_list
        rval = rval + body
        rval = rval + self._nindent + "}\n"
        return rval


    def visit_literal_node(self, node: Literal) -> str:
        if node.datatype.intrinsic == ScalarType.Intrinsic.CHARACTER:
            return "\"" + node.value + "\""
        if node.datatype.intrinsic == ScalarType.Intrinsic.BOOLEAN:
            if node.value == "True":
                return "true"
            else:
                return "false"
        return node.value

    def visit_loop_node(self, node: Loop) -> str:
        start = self._visit(node.start_expr)
        stop = self._visit(node.stop_expr)
        step = self._visit(node.step_expr)
        body = self._visit(node.body)
        string = f"for({node.variable.name} = {start}; {node.variable.name} < {stop}; {node.variable.name} += {step})" + "{\n"
        string = string + body
        string = string + self._nindent + "}"
        return string

    def visit_member_node(self, node: Member) -> str:
        return node.name

    def visit_binaryoperation_node(self, node: BinaryOperation) -> str:
        lhs = self._visit(node.children[0])
        rhs = self._visit(node.children[1])
        operator = ""
        # This should be a map
        if node.operator == BinaryOperation.BinaryOp.ADDITION:
            operator = "+"
        elif node.operator == BinaryOperation.BinaryOp.SUBTRACTION:
            operator = "-"
        elif node.operator == BinaryOperation.BinaryOp.MULTIPLY:
            operator = "*"
        elif node.operator == BinaryOperation.BinaryOp.DIVISION:
            operator = "/"
        elif node.operator == BinaryOperation.BinaryOp.LESS_THAN:
            operator = "<"
        elif node.operator == BinaryOperation.BinaryOp.LESS_THAN_EQUAL:
            operator = "<="
        elif node.operator == BinaryOperation.BinaryOp.GREATER_THAN:
            operator = ">"
        elif node.operator == BinaryOperation.BinaryOp.GREATER_THAN_EQUAL:
            operator = ">="
        elif node.operator == BinaryOperation.BinaryOp.EQUALITY:
            operator = "=="
        elif node.operator == BinaryOperation.BinaryOp.LOG_AND:
            operator = "&&"
        elif node.operator == BinaryOperation.BinaryOp.LOG_OR:
            operator = "||"

        return f"({lhs} {operator} {rhs})"

    def visit_unaryoperation_node(self, node: UnaryOperation) -> str:
        if node.operator == UnaryOperation.UnaryOp.UNARYSUB:
            return f"-{self._visit(node.children[0])}"
        elif node.operator == UnaryOperation.UnaryOp.LOG_NOT:
            return f"(!{self._visit(node.children[0])})"

    def visit_particlepositionreference_node(self, node: ParticlePositionReference) -> str:
        if not self._in_kernel:
            raise NotImplementedError()
        self.addSlice("core_part_position")
        # TODO We should check which particle this accesses when pariwise is
        # supported.
        return f"_core_part_position.access(i, a, {node.dimension})"

    def visit_particlecorereference_node(self, node: ParticleCoreReference) -> str:
        # TODO We should check which particle this accesses when pariwise is
        # supported.
        slice_name = "core_part_" + node.member.name
        self.addSlice(slice_name)

        rstr = f"_{slice_name}.access(i, a"
        if isinstance(node.member, ArrayMixin):
            indices = []
            for index in node.member.indices:
                indices.append(self._visit(index))
            rstr = rstr + ", " + ", ".join(indices) + ")"
        else:
            rstr = rstr + ")"
        return rstr

    def visit_particlereference_node(self, node: ParticleReference) -> str:
        # TODO We should check which particle this accesses when pariwise is
        # supported.
        slice_name = node.member.name
        self.addSlice(slice_name)

        rstr = f"_{slice_name}.access(i, a"
        if isinstance(node.member, ArrayMixin):
            indices = []
            for index in node.member.indices:
                indices.append(self._visit(index))
            rstr = rstr + ", " + ", ".join(indices) + ")"
        else:
            rstr = rstr + ")"
        return rstr

    def visit_configreference_node(self, node: ConfigReference) -> str:
        # TODO
        if self._in_kernel:
            return f"_{node.symbol.name}(0).{self._visit(node.member)}"
        else:
            return "config.config_host(0)." + self._visit(node.member)

    def visit_pointerreference_node(self, node: PointerReference) -> str:
        return f"*{node.symbol.name}"

    def visit_autoreference_node(self, node: AutoReference) -> str:
        return f"{node.symbol.name}"

    def visit_scalarreference_node(self, node: ScalarReference) -> str:
        return f"{node.symbol.name}"

    def visit_break_node(self, node: Break) -> str:
        return "break;"

    def visit_return_node(self, node: Return) -> str:
        if len(node.children) == 1:
            return f"return {self._visit(node.children[0])};"
        return "return;"

    def visit_structuremember_node(self, node: StructureMember) -> str:
        return f"{node.name}.{self._visit(node.member)}"

    def visit_structurereference_node(self, node: StructureReference) -> str:
        return f"{node.symbol.name}.{self._visit(node.member)}"

    def visit_while_node(self, node: While) -> str:
        cond = self._visit(node.condition)
        body = self._visit(node.body)
        string = f"while({cond})" + "{\n"
        string = string + body
        string = string + self._nindent + "}"
        return string
