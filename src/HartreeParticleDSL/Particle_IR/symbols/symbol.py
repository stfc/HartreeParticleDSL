'''
This module contains the basic Symbol abstract class, which is used as the base
class for other Symbols in Particle_IR.
'''

from __future__ import annotations

from abc import ABCMeta
from enum import Enum
from psyclone.psyir.nodes import Node

import psyclone.psyir.symbols.symbol as psySym

class Symbol(psySym.Symbol, metaclass=ABCMeta):
    '''
    Generic Symbol item for the Symbol Table and PIR References.

    :param str name: name of the symbol.
    :param visibility: visibility of the symbol.
    :type visibility: :py:class`HartreeParticleDSL.Particle_IR.symbols.Symbol.Visibility`

    :raises TypeError: if the name is not a str.
    '''

    class Visibility(Enum):
        '''
        Enumeration of the different visibilty attributes supported in the PIR.
        If no visibility information is supplied for a Symbol then it is given
        the DEFAULT_VISIBILITY.

        LOCAL: the symbol is local to the current kernel or function.
        GLOBAL: The variable is visibile globally.
        '''
        PUBLIC=psySym.Symbol.Visibility.PUBLIC
        PRIVATE=psySym.Symbol.Visibility.PRIVATE
        LOCAL=3
        GLOBAL=4

    def __init__(self, name: str, visibility: Visibility=Visibility.LOCAL) -> None:
        if not isinstance(name, str):
            raise TypeError(f"{type(self)} 'name' attribute should be of type str"
                            f" but {type(name)} found.")
        super().__init__(name=name, visibility=visibility)

    @property
    def visibility(self) -> Visibility:
        '''
        :returns: the visibility of this Symbol.
        :rtype: :py:class`HartreeParticleDSL.Particle_IR.symbols.Symbol.Visibility`
        '''
        return self._visibility

    @visibility.setter
    def visibility(self, value: Visibility) -> None:
        '''
        Setter for the visibility attribute.

        :raises TypeError: if the supplied value is not an instance of Symbol.Visibility.
        '''
        if not isinstance(value, Symbol.Visibility):
            raise TypeError(f"{type(self)} visibility attribute should be of type "
                            f"Particle_IR.symbols.Symbol.Visibility but got {type(value)}.")

        self._visibility = value

    @property
    def name(self) -> str:
        '''
        :returns: the name of this Symbol.
        :rtype: str
        '''
        return self._name

    def find_symbol_table(self, node: Node) -> SymbolTable:
        '''
        Searches back up the tree for the SymbolTable that contains this Symbol.

        :param node: the PIR node from which to search.
        :type node: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :returns: the SymbolTable containing this Symbol or None
        :rtype: :py:type:`HartreeParticleDSL.Particle_IR.symbols.symboltable.SymbolTable` \
                or None

        :raises TypeError: if the supplied node is not a PIR Node.
        '''
        # pylint: disable=fixme
        if not isinstance(node, Node):
            raise TypeError(f"find symbol table expected to be passed an instance of "
                            f"Particle_IR.nodes.Node but got {type(node)}.")

        #TODO If/When needed
        raise NotImplementedError()
