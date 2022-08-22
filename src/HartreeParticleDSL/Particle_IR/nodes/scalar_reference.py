from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol

class ScalarReference(Reference):
    '''
    Contains a ScalarReference to a variable in the ParticleIR tree.

    :param symbol: The scalar type symbol that this reference points to.
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol.ScalarTypeSymbol`

    '''

    def __init__(self, symbol: ScalarTypeSymbol) -> None:
        super().__init__()
        self.symbol = symbol
    
    @symbol.setter
    def symbol(self, symbol: ScalarTypeSymbol) -> None:
        '''
        Sets the symbol this ScalarReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol.ScalarTypeSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, ScalarTypeSymbol):
            raise TypeError("Attempted to make a ScalarReference to a non-ScalarTypeSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol
