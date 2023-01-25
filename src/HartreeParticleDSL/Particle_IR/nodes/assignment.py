'''
This module contaisn the Assignment class.
'''

from __future__ import annotations

from typing import List, Union

from HartreeParticleDSL.Particle_IR.nodes.node import DataNode, Node
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement

class Assignment(Statement):
    '''
    Class to represent an Assignment. Assignments have 2 datanode
    children, which can be accessed as lhs and rhs.

    :param children: List of Nodes to be contained in this Assignment region.
    :type children: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node` \
            or None.
    '''
    def __init__(self, children: Union[List[Node], None]=None) -> None:
        super().__init__(children=children)

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given position and node are valid as a child
        of this node.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if (position in (0, 1)) and isinstance(child, DataNode):
            return True
        return False

    @property
    def lhs(self) -> DataNode:
        '''
        :returns: The lhs of this assignment (child 0).
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        '''
        return self.children[0]

    @property
    def rhs(self) -> DataNode:
        '''
        :returns: the rhs of this assignment (child 1).
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        '''
        return self.children[1]

    @staticmethod
    def create(lhs: DataNode, rhs: DataNode) -> Assignment:
        '''
        Creates an assignment with the given lhs and rhs nodes.

        :param lhs: The lhs of the assignment to be created.
        :type lhs: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        :param rhs: The rhs of the assignment to be created.
        :type rhs: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode

        :returns: An assignment node representing the lhs and rhs provided.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.assignment.Assignment`
        '''
        assign = Assignment()
        assign.addchild(lhs)
        assign.addchild(rhs)
        return assign

    def node_str(self) -> str:
        '''
        :returns: a text description of this assignment
        :rtype: str
        '''
        nodestr = f"Assignment[{self.lhs.node_str()}, {self.rhs.node_str()}]"
        return nodestr
