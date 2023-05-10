'''
This module contains the Body class.
'''

from __future__ import annotations
from typing import Union, List

from HartreeParticleDSL.Particle_IR.nodes.node import Node
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement

import psyclone.psyir.nodes.schedule as psysched

class Body(Node, psysched.Schedule):
    '''
    Class to represent a Body of a code region. Can contain any
    number of Statement nodes.

    :param children: List of Nodes to be contained in this Body region.
    :type children: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node` \
            or None.
    '''
    # pylint: disable=undefined-variable

    def __init__(self, children: Union[None,List[Node]]=None) -> None:
        super().__init__(children=children)

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given position and node are valid as a child
        of this node.

        Body nodes can have any number of Statement nodes as children.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if isinstance(child, Statement):
            return True
        return False

    def node_str(self) -> str:
        '''
        :returns: a text description of this assignment
        :rtype: str
        '''
        nodestr = "Body[\n"
        for i in self.children:
            nodestr += f"    {i.node_str()}\n"
        nodestr +="] End Body"
        return nodestr
