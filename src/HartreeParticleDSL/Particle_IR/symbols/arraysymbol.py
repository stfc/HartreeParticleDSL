
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import ArrayType

class ArraySymbol(Symbol):
    '''
    Class for those symbols that are to an Array Type.

    :param str name: name of the symbol.
    :param datatype: datatype of the symbol.
    :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ArrayType`
    :param kwargs: additional arguments provided by the Symbol base class.
    :type kwargs: unwrapped dict.
    '''

    def __init__(self, name: str, datatype : ArrayType, **kwargs):
        self._datatype = None
        super().__init__(name, **kwargs)
        self.datatype = datatype

    @property
    def datatype(self) -> ArrayType:
        '''
        :returns: datatype of the StructureSymbol.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ArrayType`
        '''
        return self._datatype

    @datatype.setter
    def datatype(self, value: ArrayType) -> None:
        '''
        Setter for the datattype of a StructureSymbol.

        :param value: new value for the datatype.
        :type value: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ArrayType`

        :raises TypeError: if the value is not of the correct type.
        '''

        # TODO
        if not isinstance(value, ArrayType):
            raise TypeError(
                    f" The datatype of a {type(self)} must be specified "
                    f" using a ArrayType but got {type(value)}")
        self._datatype = value
