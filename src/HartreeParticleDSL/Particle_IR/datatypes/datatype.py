'''
This module contains all of the data type classes used to describe types in
Particle IR.
'''
from __future__ import annotations

import abc
from abc import ABCMeta
from collections import OrderedDict
from enum import Enum
from typing import Union, Dict, List, Tuple
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

import psyclone.psyir.symbols.datatypes as psyDT
from psyclone.psyir.symbols.datatypes import DataType, ScalarType, StructureType
from psyclone.psyir.symbols import Symbol

class NoType(DataType):
    ''' Empty datatype (e.g. void type).'''
    # pylint: disable=too-few-public-methods

    def __str__(self):
        return "NoType"

class PointerType(DataType):
    '''
    Descibes a pointer to a datatype, instead of a basic datatype.
    Implementing this in a backend is not required, but may be useful sometimes.

    :param datatype: The datatype to point to.
    :type datatype: :py:class:`psyclone.Particle_IR.datatypes.datatype.DataType`

    :raises TypeError: if datatype is not a DataType.
    '''
    # pylint: disable=too-few-public-methods

    def __init__(self, datatype: DataType) -> None:
        if not isinstance(datatype, DataType):
            raise TypeError("Attempted to make a PointerType with datatype "
                            f"of type {type(datatype)} instead of a DataType.")
        self._datatype = datatype

    def __str__(self) -> str:
        '''
        :returns: a description of this type.
        :rtype: str
        '''
        return f"PointerType<{self._datatype}>"

    def __eq__(self, other):
        return self is other

#    def __hash__(self):
#        return hash(id(self))

class ArrayType(DataType):
    '''
    Describes an array type.

    :param datatype: The datatype this array type is of.
    :type datatype: :py:class:`psyclone.Particle_IR.datatypes.datatype.DataType`
    :param list shape: shape of the symbol in row-major order (rightmost index is \
                       contiguous in memory - IE c-style). If is \
                       ArrayType.Extent.DYNAMIC it assumes the extent is unknown \
                       and thus is a run-time parameter. All indices are assumed \
                       to start from 0 (C-Style) and run up to index-1.

    :raises TypeError: if any of the arguments are of the wrong type.
    '''
    class Extent(Enum):
        '''
        Enumeration of array shape extents that are unspecified at compile time.
        '''
        DYNAMIC = psyDT.ArrayType.Extent.DEFERRED

    def __eq__(self, other):
        eq = type(self) == type(other)
        if eq:
            eq = eq and self._datatype == other._datatype
            eq = eq and len(self.shape) == len(other.shape)
            if eq:
                for i, el in enumerate(self.shape):
                    eq = eq and el == other.shape[i]
        return eq

#    def __hash__(self):
#        return hash((self._datatype, tuple(self.shape)))

    def __init__(self, datatype: DataType, shape: List[Union[int, Extent]]) -> None:
        if not isinstance(datatype, DataType):
            raise TypeError("Attempted to make an ArrayType with datatype "
                            f"of type {type(datatype)} instead of a DataType.")
        self._datatype = datatype

        self._shape = []
        if not isinstance(shape, list) or len(shape) < 1:
            raise TypeError("shape argument needs to be a nonempty list of int/Extents but found "
                            f"{type(shape)}.")
        for child in shape:
            if not isinstance(child, (int, ArrayType.Extent)):
                raise TypeError("Each member of the shape argument needs to be an int "
                                f"or Extent but found {type(child)}.")
            self._shape.append(child)

    @property
    def shape(self) -> List[Union[int, Extent]]:
        '''
        :returns: the shape of this ArrayType.
        :rtype: List of int or Extent.
        '''
        return self._shape

    def __str__(self) -> str:
        '''
        :returns: a string representation of this array datatype.
        :rtype: str
        '''

        dims = []
        for dimension in self.shape:
            if isinstance(dimension, ArrayType.Extent):
                dims.append("*")
            else:
                dims.append(f"{dimension}")
        dim_string = ", ".join(dims)
        return f"ArrayType<{self._datatype}: [{dim_string}]>"

# Create common datatype
FLOAT_TYPE = ScalarType(ScalarType.Intrinsic.REAL, ScalarType.Precision.SINGLE)
DOUBLE_TYPE = ScalarType(ScalarType.Intrinsic.REAL, ScalarType.Precision.DOUBLE)

INT_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, ScalarType.Precision.SINGLE)
LONG_LONG_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, ScalarType.Precision.DOUBLE)
INT32_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 32)
INT64_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 64)
INT8_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 8)

BOOL_TYPE = ScalarType(ScalarType.Intrinsic.BOOLEAN, ScalarType.Precision.UNDEFINED)

STRING_TYPE = ScalarType(ScalarType.Intrinsic.CHARACTER, ScalarType.Precision.UNDEFINED)

# Particles always have a position and velocity, and its always considered to be a 3D double array
# internally for now.
PARTICLE_POSITION_TYPE = ArrayType(DOUBLE_TYPE, [3])

BASE_PARTICLE_TYPE = StructureType()
_CORE_PART_TYPE = StructureType()
_CORE_PART_TYPE.add("position", PARTICLE_POSITION_TYPE, Symbol.Visibility.PUBLIC)
_CORE_PART_TYPE.add("velocity", ArrayType( DOUBLE_TYPE, [3]), Symbol.Visibility.PUBLIC)
BASE_PARTICLE_TYPE.add("core_part", _CORE_PART_TYPE, Symbol.Visibility.PUBLIC)

BASE_CONFIG_TYPE = StructureType()
_SPACE_TYPE = StructureType()
_BOUNDARY_TYPE = StructureType()
_BOUNDARY_TYPE.add("x_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("x_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("y_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("y_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("z_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("z_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_x_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_x_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_y_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_y_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_z_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_BOUNDARY_TYPE.add("local_z_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
_SPACE_TYPE.add("box_dims", _BOUNDARY_TYPE, Symbol.Visibility.PUBLIC)
BASE_CONFIG_TYPE.add("space", _SPACE_TYPE, Symbol.Visibility.PUBLIC)
BASE_CONFIG_TYPE.add("nparts", INT64_TYPE, Symbol.Visibility.PUBLIC)

type_mapping_str = {"c_int": INT_TYPE,
                    "int": INT_TYPE,
                    "c_double": DOUBLE_TYPE,
                    "double": DOUBLE_TYPE,
                    "c_float": FLOAT_TYPE,
                    "float": FLOAT_TYPE,
                    "c_int64_t": INT64_TYPE,
                    "int64_t": INT64_TYPE,
                    "c_int32_t": INT32_TYPE,
                    "int32_t": INT32_TYPE,
                    "c_int8_t": INT8_TYPE,
                    "int8_t": INT8_TYPE,
                    "c_bool": BOOL_TYPE,
                    "bool": BOOL_TYPE,
                    "char*": STRING_TYPE,
                    "string": STRING_TYPE,
                    "part": BASE_PARTICLE_TYPE,
                    "config": BASE_CONFIG_TYPE}

def reset_part_and_config():
    '''
    This function resets the BASE_PARTICLE_TYPE and BASE_CONFIG_TYPE to their
    defaults.
    '''
    # pylint: disable=global-statement, invalid-name
    global BASE_PARTICLE_TYPE, BASE_CONFIG_TYPE
    BASE_PARTICLE_TYPE = StructureType()
    _CORE_PART_TYPE = StructureType()
    _CORE_PART_TYPE.add("position", PARTICLE_POSITION_TYPE, Symbol.Visibility.PUBLIC)
    _CORE_PART_TYPE.add("velocity", ArrayType( DOUBLE_TYPE, [3]), Symbol.Visibility.PUBLIC)
    BASE_PARTICLE_TYPE.add("core_part", _CORE_PART_TYPE, Symbol.Visibility.PUBLIC)
    
    BASE_CONFIG_TYPE = StructureType()
    _SPACE_TYPE = StructureType()
    _BOUNDARY_TYPE = StructureType()
    _BOUNDARY_TYPE.add("x_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("x_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("y_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("y_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("z_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("z_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_x_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_x_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_y_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_y_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_z_min", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _BOUNDARY_TYPE.add("local_z_max", DOUBLE_TYPE, Symbol.Visibility.PUBLIC)
    _SPACE_TYPE.add("box_dims", _BOUNDARY_TYPE, Symbol.Visibility.PUBLIC)
    BASE_CONFIG_TYPE.add("space", _SPACE_TYPE, Symbol.Visibility.PUBLIC)
    BASE_CONFIG_TYPE.add("nparts", INT64_TYPE, Symbol.Visibility.PUBLIC)
    type_mapping_str["part"] = BASE_PARTICLE_TYPE
    type_mapping_str["config"] = BASE_CONFIG_TYPE


def reset_type_mapping_str():
    '''
    This function resets the type_mapping_str to its initial state.
    '''
    valid_keys = ["c_int","int", "c_double", "double",
            "c_float", "float", "c_int64_t", "int64_t", "c_int32_t",
            "int32_t"," c_int8_t", "int8_t", "c_bool", "bool",
            "char*", "string", "part", "config"]
    dict_keys = type_mapping_str.keys()
    to_remove = []
    for key in dict_keys:
        if key not in valid_keys:
            to_remove.append(key)
    for key in to_remove:
        del type_mapping_str[key]
