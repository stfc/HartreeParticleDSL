'''
This module contains the abstract Reference class.
'''

from abc import ABCMeta
from HartreeParticleDSL.Particle_IR.nodes.node import DataNode
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol

import psyclone.psyir.nodes.reference as psyRef

class Reference(DataNode, psyRef.Reference, metaclass=ABCMeta):
    '''
    Contains a Reference to a variable in the ParticleIR tree.
    '''
    def __init__(self, symbol=None):
        super().__init__(symbol=symbol)

        self._symbol = None

    @property
    def symbol(self) -> Symbol:
        '''
        :returns: the symbol this Reference refers to.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
        '''
        return self._symbol

    @symbol.setter
    def symbol(self, symbol: Symbol) -> None:
        '''
        Sets the symbol this Reference refers to. This is specialised for
        each subclass of Reference.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
        '''
        self._symbol = symbol

    @property
    def name(self) -> str:
        '''
        :returns: the name of the reference symbol
        :rtype: str
        '''
        return self._symbol.name

    @property
    def is_array(self) -> bool:
        '''
        :returns: if this reference is an array.
        :rtype: bool
        '''
        return False

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"{type(self).__name__}[{self.name}]"
        return nodestr
