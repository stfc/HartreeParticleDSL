'''
This module contains the ScalarTypeSymbol which represents scalar types
in the Particle_IR tree.
'''
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType

class ScalarTypeSymbol(Symbol):
    '''
    Class for symbols to a ScalarType.

    :param str name: name of the symbol.
    :param datatype: datatype of the symbol.
    :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ScalarType`
    :param kwargs: additional arguments provided by the Symbol base class.
    :type kwargs: unwrapped dict.
    '''

    def __init__(self, name: str, datatype : ScalarType, **kwargs):
        self._datatype = None
        super().__init__(name, **kwargs)
        self.datatype = datatype

    @property
    def datatype(self) -> ScalarType:
        '''
        :returns: datatype of the TypedSymbol.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ScalarType`
        '''
        return self._datatype

    @datatype.setter
    def datatype(self, value: ScalarType) -> None:
        '''
        Setter for the datattype of a Typedsymbol.

        :param value: new value for the datatype.
        :type value: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.ScalarType`

        :raises TypeError: if the value is not of the correct type.
        '''
        if not isinstance(value, ScalarType):
            raise TypeError(
                    f"The datatype of a {type(self)} must be specified "
                    f"using a ScalarType but got {type(value)}.")
        self._datatype = value
