'''
This module contains the PointerSymbol class to represent Pointers in the
Particle_IR tree.
'''
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import PointerType

class PointerSymbol(Symbol):
    '''
    Class for those symbols that are to a Pointer Type.

    :param str name: name of the symbol.
    :param datatype: datatype of the symbol.
    :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.PointerType`
    :param kwargs: additional arguments provided by the Symbol base class.
    :type kwargs: unwrapped dict.
    '''

    def __init__(self, name: str, datatype : PointerType, **kwargs):
        self._datatype = None
        super().__init__(name, **kwargs)
        self.datatype = datatype

    @property
    def datatype(self) -> PointerType:
        '''
        :returns: datatype of the PointerSymbol.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.PointerType`
        '''
        return self._datatype

    @datatype.setter
    def datatype(self, value: PointerType) -> None:
        '''
        Setter for the datattype of a PointerSymbol.

        :param value: new value for the datatype.
        :type value: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.PointerType`

        :raises TypeError: if the value is not of the correct type.
        '''

        if not isinstance(value, PointerType):
            raise TypeError(
                    f"The datatype of a {type(self)} must be specified "
                    f"using a PointerType but got {type(value)}.")
        self._datatype = value
