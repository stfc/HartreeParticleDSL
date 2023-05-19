'''
This module contains the AutoReference class.
'''

from psyclone.psyir.nodes import Reference
from HartreeParticleDSL.Particle_IR.symbols.autosymbol import AutoSymbol

class AutoReference(Reference):
    '''
    Contains a AutoReference to a variable in the ParticleIR tree.

    :param symbol: The scalar type symbol that this reference points to.
    :type symbol: \
            :py:class:`HartreeParticleDSL.Particle_IR.symbols.autosymbol.AutoSymbol`

    '''
    _text_name = "AutoReference"

    def __init__(self, symbol: AutoSymbol) -> None:
        super().__init__(symbol=symbol)
        self.symbol = symbol

    @property
    def symbol(self) -> AutoSymbol:
        return super().symbol

    @symbol.setter
    def symbol(self, symbol: AutoSymbol) -> None:
        '''
        Sets the symbol this AutoReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: \
                :py:class:`HartreeParticleDSL.Particle_IR.symbols.autosymbol.AutoSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, AutoSymbol):
            raise TypeError("Attempted to make an AutoReference to a non-AutoSymbol. "
                            f"Got {type(symbol)} as input.")
        # Weird workaround to use the parent setter. See
        # https://stackoverflow.com/questions/1021464/how-to-call-a-property-of-the-base-class-if-this-property-is-being-overwritten-i
        super(AutoReference, self.__class__).symbol.fset(self, symbol)
