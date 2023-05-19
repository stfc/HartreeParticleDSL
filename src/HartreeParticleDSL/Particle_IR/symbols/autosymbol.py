'''
This module contains the AutoSymbol which represents types
that would be expressed by the C++ auto keyword in the Particle IR
tree.
'''
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from psyclone.psyir.symbols.datatypes import ScalarType

class AutoSymbol(Symbol):
    '''
    Class for symbols to an auto type.

    :param str name: name of the symbol.
    :param kwargs: additional arguments provided by the Symbol base class.
    :type kwargs: unwrapped dict.
    '''

    def __init__(self, name: str, initial_value: str, **kwargs):
        self._datatype = None
        self._initial_value = None
        super().__init__(name, **kwargs)
        self.initial_value = initial_value

    @property
    def initial_value(self) -> str:
        '''
        :returns: initial_value of the AutoSymbol.
        :rtype: str
        '''
        return self._initial_value

    @initial_value.setter
    def initial_value(self, value: str) -> None:
        '''
        Setter for the initial value of an AutoSymbol

        :param value: new inital value for the AutoSymbol
        :type value: str

        :raises TypeError: if the value is not of the correct type.
        '''
        if not isinstance(value, str):
            raise TypeError(
                    f"The initial_value of a {type(self)} must be specified "
                    f"using a str but got {type(value)}.")
        self._initial_value = value
