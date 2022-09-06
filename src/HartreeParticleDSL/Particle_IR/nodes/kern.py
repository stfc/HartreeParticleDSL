from __future__ import annotations

from abc import ABCMeta
from typing import List
from HartreeParticleDSL.Particle_IR.nodes.node import Node
from HartreeParticleDSL.Particle_IR.nodes.body import Body

class Kern(Node, metaclass=ABCMeta):
    '''
    Abstract class representing a general Kernel object.
    
    :param children: List of Nodes to be contained in this Kern region.
    :type children: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node.
    '''

    def __init__(self, children: List[Node]=None) -> None:
        super().__init__(children=children)

        # If no children were specified, we need to add the Body of this Kern
        if not children or len(children)==0:
            self.addchild(Body())

        # Specialisations of this class will need argument lists
        self._arguments = []

        # Add a symbol table
        from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable
        self._symbol_table = SymbolTable(kern=self)

    @property
    def symbol_table(self) -> SymbolTable:
        '''
        :returns: The symbol table for this Kernel.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol_table.SymbolTable`
        '''
        return self._symbol_table

    @property
    def body(self) -> Body:
        '''
        :returns: The body of this Kernel.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.body.Body`
        '''
        return self.children[0]

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given position and node are valid as a child
        of this node.

        Kern nodes can contain a single Statement node as its child.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if position == 0 and isinstance(child, Body):
            return True
        return False

    def node_str(self) -> str:
        '''
        :returns: a text description of this assignment
        :rtype: str
        '''
        arg_strings = []
        for x in self._arguments:
            arg_strings.append(x.node_str())
        arg_string = ", ".join(arg_strings)
        nodestr = type(self).__name__ + f"[{arg_string}: {self.children[0].node_str()}]"
        return nodestr
