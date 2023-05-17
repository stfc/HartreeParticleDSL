'''
This module contaisn the Assignment class.
'''

from __future__ import annotations

from typing import List, Union

from HartreeParticleDSL.Particle_IR.nodes.statement import Statement

import psyclone.psyir.nodes.assignment as psyassign

class Assignment(psyassign.Assignment, Statement):
    _text_name = "Assignment"

    def node_str(self) -> str:
        '''
        :returns: a text description of this assignment
        :rtype: str
        '''
        nodestr = f"Assignment[{self.lhs.node_str()}, {self.rhs.node_str()}]"
        return nodestr

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
