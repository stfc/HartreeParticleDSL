'''
This module contains the PointerReference class.
'''

from __future__ import annotations

from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol

class PointerReference(Reference):
    '''
    Contains a PointerReference to a variable in the ParticleIR tree.

    :param symbol: The scalar type symbol that this reference points to.
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.pointersymbol.PointerSymbol`

    '''
    # pylint: disable=undefined-variable

    def __init__(self, symbol: PointerSymbol) -> None:
        super().__init__(symbol=symbol)
        self.symbol = symbol

    @property
    def symbol(self) -> Symbol:
        return super().symbol

    @symbol.setter
    def symbol(self, symbol: PointerSymbol) -> None:
        '''
        Sets the symbol this PointerReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.pointersymbol.PointerSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, PointerSymbol):
            raise TypeError("Attempted to make a PointerReference to a non-PointerSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol
