from typing import List

from HartreeParticleDSL.Particle_IR.nodes.array_mixin import ArrayMixin
from HartreeParticleDSL.Particle_IR.nodes.node import DataNode
from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol

class ArrayReference(ArrayMixin, Reference):
    '''
    Contains an ArrayReference to an ArraySymbol

    :param symbol: The array symbol referenced by this reference
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.arraysymbol.ArraySymbol`
    :param indices: The indices used to access this array
    :type indices: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
    '''

    def __init__(self, symbol: ArraySymbol, indices: List[DataNode]) -> None:
        super().__init__()
        self.symbol = symbol
        
        for index in indices:
            self.addchild(index)

    @symbol.setter
    def symbol(self, symbol: ArraySymbol) -> None:
        '''
        Sets the symbol this ArrayReference points to

        :param symbol: The symbol this Reference references.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.arraysymbol.ArraySymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, ArraySymbol):
            raise TypeError("Attempted to make an ArrayReference to a non-ArraySymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol

    def node_str(self) -> str:
        '''
        :returns: a text description of this assignment
        :rtype: str
        '''
        nodestr = f"ArrayReference[{self.symbol.name}: ("
        indices = []
        for i in self.children:
            indices.append(i.node_str())
        indices_str = ", ".join(indices)
        nodestr = nodestr + f"{indices_str})]"
        return nodestr
