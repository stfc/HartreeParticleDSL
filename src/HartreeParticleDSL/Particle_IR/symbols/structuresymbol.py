'''
This module contains the StructureSymbol, which represents StructureTypes
in the Particle_IR tree.
'''
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import StructureType

class StructureSymbol(Symbol):
    '''
    Class for those symbols that are to a Structure Type.

    :param str name: name of the symbol.
    :param datatype: datatype of the symbol.
    :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.StructureType`
    :param kwargs: additional arguments provided by the Symbol base class.
    :type kwargs: unwrapped dict.
    '''

    def __init__(self, name: str, datatype : StructureType, **kwargs) -> None:
        self._datatype = None
        super().__init__(name, **kwargs)
        self.datatype = datatype

    @property
    def datatype(self) -> StructureType:
        '''
        :returns: datatype of the StructureSymbol.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.StructureType`
        '''
        return self._datatype

    @datatype.setter
    def datatype(self, value: StructureType) -> None:
        '''
        Setter for the datattype of a StructureSymbol.

        :param value: new value for the datatype.
        :type value: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.StructureType`

        :raises TypeError: if the value is not of the correct type.
        '''

        if not isinstance(value, StructureType):
            raise TypeError(
                    f"The datatype of a {type(self)} must be specified "
                    f"using a StructureType but got {type(value)}.")
        self._datatype = value
