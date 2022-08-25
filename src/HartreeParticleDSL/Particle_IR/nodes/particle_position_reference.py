from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

class ParticlePositionReference(Reference):
    '''
    Contains a ParticlePositionReference in the IR tree.

    :param symbol: The structure type symbol that this reference points to.
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`
    :param int dimension: The dimension of the particle position accessed in this node.
    '''
    def __init__(self, symbol: StructureSymbol, dimension: int) -> None:
        super.__init__()
        self.symbol = symbol
        self.dimension = dimension

    @symbol.setter
    def symbol(self, symbol: StructureSymbol) -> None:
        '''
        Sets the symbol this ParticlePositionReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, StructureSymbol):
            raise TypeError("Attempted to make a ParticlePositionReference to a non-StructureSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol

    @property
    def dimension(self) -> int:
        '''
        :returns: the dimension accessed in this ParticlePositionReference
        :rtype: int
        '''
        return self._dimension

    @dimension.setter
    def dimension(self, dimension: int) -> None:
        '''
        Sets the dimension this ParticlePositionReference accesses.

        :param int dimension: the dimension accessed in this particle position reference.

        :raises TypeError: if the member is not the correct type.
        '''
        if not isinstance(dimension, int):
            raise TypeError("Attempted to make a ParticlePositionReference with a non-int dimension. "
                            f" Got {type(dimension)} as input.")
        self._dimension = dimension

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"ParticlePositionReference[{self.symbol.name}: {self.dimension}]"
        return nodestr
