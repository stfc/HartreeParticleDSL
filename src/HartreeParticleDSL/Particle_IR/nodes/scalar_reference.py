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
 
    @property
    def symbol(self) -> ScalarTypeSymbol:
        return super().symbol

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
        # Weird workaround to use the parent setter. See
        # https://stackoverflow.com/questions/1021464/how-to-call-a-property-of-the-base-class-if-this-property-is-being-overwritten-i
        super(ScalarReference, self.__class__).symbol.fset(self, symbol)
