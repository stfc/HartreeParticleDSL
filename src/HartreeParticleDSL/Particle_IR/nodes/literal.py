'''
This module contains the Literal class.
'''

from __future__ import annotations

import re

from psyclone.psyir.nodes import DataNode
import psyclone.psyir.nodes.literal as psyLit
from psyclone.psyir.symbols.datatypes import ScalarType

class Literal(DataNode):
#class Literal(psyLit.Literal):
    '''
    Node representing a Literal in the Particle IR tree.

    :param str value: The string representation of the value of this Literal.
    :param datatype: The datatype of this Literal.
    :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ScalarType`

    :raises TypeError: If the datatype is not a valid ScalarType.
    :raises TypeError: If the supplied value is not a string.
    :raises ValueError: If the datatype is an INTEGER but the value contains \
                        non-number characters.
    :raises ValueError: If the datatype is a REAL but the value is not a real \
                        representation.
    :raises ValueError: If the datatype is a BOOLEAN but the value is not True \
                        or False.
    '''
    _real_value = r'^[+-]?[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?$'
    _int_value = r'([+-]?[1-9][0-9]*|0)'
    _text_name = "Literal"
    _colour = "Yellow"

    def __init__(self, value: str, datatype: ScalarType) -> None:
        super().__init__()

        if not isinstance(datatype, ScalarType):
            raise TypeError("Literal datatype must be a ScalarType but "
                            f"{type(datatype)} supplied.")

        if not isinstance(value, str):
            raise TypeError("Literal value must be a string but "
                            f"{type(value)} supplied.")

        if (datatype.intrinsic == ScalarType.Intrinsic.INTEGER and
            not re.fullmatch(Literal._int_value, value)):
            raise ValueError("Constructing integer Literal but got a value of "
                             f"'{value}' instead of an integer value.")

        if (datatype.intrinsic == ScalarType.Intrinsic.REAL and
            not re.fullmatch(Literal._real_value, value)):
            raise ValueError("Constructing float Literal but got a value of "
                             f"'{value}' instead of a float value.")


        if (datatype.intrinsic == ScalarType.Intrinsic.BOOLEAN and
            value not in ("True", "False")):
            raise ValueError("Constructing boolean Literal but got a value of "
                             f"'{value}' instead of True or False.")
        if datatype.intrinsic == ScalarType.Intrinsic.REAL:
            self._value = value.replace("E", "e")
        else:
            self._value = value
        self._datatype = datatype

    @property
    def datatype(self) -> ScalarType:
        '''
        :returns: The datatype of this Literal.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ScalarType`
        '''
        return self._datatype

    @property
    def value(self) -> str:
        '''
        :returns: The string representation of the value of this Literal.
        :rtype: str
        '''
        return self._value

    def node_str(self, colour=True) -> str:
        '''
        :returns: The string representation of this node in the tree.
        :rtype: str
        '''
        return self.coloured_name(colour) + f"['{self.value}', {self.datatype}]"
