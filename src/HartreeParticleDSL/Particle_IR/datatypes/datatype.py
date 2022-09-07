from __future__ import annotations

import abc
from abc import ABCMeta
from collections import OrderedDict
from enum import Enum
import inspect
from typing import Union, Dict, List
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

class DataType(metaclass=ABCMeta):
    '''Abstract base class from which all types are derived.'''

    @abc.abstractmethod
    def __str__(self):
        '''
        :returns: a description of this type.
        :rtype: str
        '''


class NoType(DataType):
    ''' Empty datatype (e.g. void type).'''

    def __str__(self):
        return "NoType"

class ScalarType(DataType):
    '''
    Describes a scalar datatype and its precision.
    
    :param intrinsic: the intrinsic of this scalar type.
    :type intrinsic: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.ScalarType.Intrinsic`
    :param precision: the precision of this scalar type.
    :type precision: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.ScalarType.Precision` or int

    :raises TypeError: if any of the arguments are of the wrong type.
    '''

    class Intrinsic(Enum):
        '''
        Enumeration of the different intrinsic scalar datatypes supported by ParticleIR.
        '''
        INTEGER = 1
        FLOAT = 2
        BOOLEAN = 3

    class Precision(Enum):
        '''
        Enumeration of the different types of precision that may be specified.
        '''
        SINGLE = 1
        DOUBLE = 2
        UNDEFINED = 3

    def __init__(self, intrinsic: Intrinsic, precision: Union[Precision, int]) -> None:
        if not isinstance(intrinsic, ScalarType.Intrinsic):
            raise TypeError(
                f"Expected 'intrinsic' to be of type ScalarType.Intrinsic "
                f"but found {type(intrinsic)}.")
        if not isinstance(precision, (int, ScalarType.Precision)):
            raise TypeError(
                f"Expected 'precision' to be of type ScalarType.Precision or int "
                f"but found {type(precision)}.")


        self._intrinsic = intrinsic
        self._precision = precision

    @property
    def intrinsic(self) -> Intrinsic:
        '''
        :returns: the intrinsic used by this scalar type.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.ScalarType.Intrinsic`
        '''
        return self._intrinsic

    @property
    def precision(self) -> Union[Precision, int]:
        ''''
        :returns: the precision used by this scalar type.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.ScalarType.Precision` or int
        '''
        return self._precision

    def __str__(self) -> str:
        '''
        :returns: a description of this scalar datatype
        :rtype: str
        '''
        precision = ""
        if isinstance(self.precision, ScalarType.Precision):
            precision = self.precision.name
        else:
            precision = f"{self.precision}"
        return f"Scalar<{self.intrinsic.name}, {precision}>"

class StructureType(DataType):
    '''
    Describes a 'structure' type that is composed of a dict of other datatypes.
    '''

    def __init__(self) -> None:
        self._components = OrderedDict()
    
    def __str__(self) -> str:
        '''
        :returns: a string representation of this structure type
        :rtype: str
        '''
        s = "StructureType<"
        comp_strs = []
        for comp in self._components.keys():
            comp_strs.append(f"({comp}: {self._components[comp]})")
        comp_str = ", ".join(comp_strs)
        s = s + comp_str + ">"
        return s

    @property
    def components(self) -> Dict[str, DataType]:
        '''
        :returns: the components that make up this structured type.
        :rtype: :py:class:`collections.OrderedDict`
        '''
        return self._components

    def add(self, name: str, typ: DataType) -> None:
        '''
        Add a child type to this StructureType

        :param name: The name of the member of this structure type.
        :type name: str
        :param typ: The DataType of the member to be added.
        :type typ: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.DataType`

        :raises TypeError: if any of the arguments are the wrong type.
        '''
        if not isinstance(name, str):
            raise TypeError(f"Expected name argument to be a string but got "
                            f"{type(name)}.")
        if not isinstance(typ, DataType):
            raise TypeError(f"Expected typ argument to be a DataType but got "
                            f"{type(typ)}.")
        if name in self._components.keys():
            raise IRGenerationError("names in a StructureType must be unique "
                                    f"but provided duplicate {name}.")
        self._components[name] = typ

    @staticmethod
    def create(components: List[Tuple[str, DataType]]) -> StructureType:
        '''
        Creates a StructureType from the supplied components.

        :param components: the name and type of each component of the new structure.
        :type components: list of 2-tuple(str, DataType)

        :returns: the new StructureType
        :rtype: :py:class:`psyclone.Particle_IR.datatypes.datatype.StructureType`

        :raise TypeError: if the members of components are not length 2.
        '''
        stype = StructureType()
        for component in components:
            if len(component) != 2:
                raise TypeError("Each component must be specified using a 2-tuple "
                                f"of (name, type) but found: {component}")
            stype.add(component[0], component[1])
        return stype

    def lookup(self,name: str) -> Union[DataType, None]:
        '''
        :returns: the DataType describing the named member of this StructureType or None.
        :rtype: :py:class:`psyclone.Particle_IR.datatypes.datatype.DataType` or None
        '''
        return self._components.get(name)

class PointerType(DataType):
    '''
    Descibes a pointer to a datatype, instead of a basic datatype.
    Implementing this in a backend is not required, but may be useful sometimes.

    :param datatype: The datatype to point to.
    :type datatype: :py:class:`psyclone.Particle_IR.datatypes.datatype.DataType`

    :raises TypeError: if datatype is not a DataType.
    '''

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
        DYNAMIC = 1

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
FLOAT_TYPE = ScalarType(ScalarType.Intrinsic.FLOAT, ScalarType.Precision.SINGLE)
DOUBLE_TYPE = ScalarType(ScalarType.Intrinsic.FLOAT, ScalarType.Precision.DOUBLE)

INT_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, ScalarType.Precision.SINGLE)
LONG_LONG_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, ScalarType.Precision.DOUBLE)
INT32_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 32)
INT64_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 64)
INT8_TYPE = ScalarType(ScalarType.Intrinsic.INTEGER, 8)

BOOL_TYPE = ScalarType(ScalarType.Intrinsic.BOOLEAN, ScalarType.Precision.UNDEFINED)

# Particles always have a position, and its always considered to be a 3D double array
# internally for now.
PARTICLE_POSITION_TYPE = ArrayType(DOUBLE_TYPE, [3])
