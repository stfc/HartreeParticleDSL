'''
This module contains the abstract Operation class, and its direct descendants,
BinaryOperation and UnaryOperation.
'''

from __future__ import annotations

from abc import ABCMeta
from enum import Enum
from typing import List

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from psyclone.psyir.nodes import DataNode, Node

from psyclone.psyir.nodes import Operation
import psyclone.psyir.nodes.operation as psyOP

#class BinaryOperation(Operation):
#    '''
#    Class representing a BinaryOperation. Has exactly two children to represent
#    the two sides of the operation.
#    '''
#    _text_name = "BinaryOperation"
#    class BinaryOp(Enum):
#        '''Enumeration of Binary Operations supported in Particle DSL.'''
#        # Arithmetic operations
#        ADDITION=1
#        SUBTRACTION=2
#        MULTIPLY=3
#        DIVISION=4
#        # Comparison operations
#        LESS_THAN=5
#        LESS_THAN_EQUAL=6
#        GREATER_THAN=7
#        GREATER_THAN_EQUAL=8
#        EQUALITY=9
#        #Logical operations
#        LOG_AND=10
#        LOG_OR=11
#
#    Operator = BinaryOp
#
#    @staticmethod
#    def create(operator: BinaryOp, children: List[DataNode]) -> BinaryOperation:
#        '''
#        Create a binary operation representing `children[0] operation children[1]`.
#
#        :param operator: The binary operation represented by the created node.
#        :type operator: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperation.BinaryOp`
#        :param children: The children represented by this BinaryOperation.
#        :type children: list of py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
#
#        :raises IRGenerationError: if the children input doesn't contain exactly 2 children \
#                                   DataNodes.
#
#        :returns: The new binary operation node.
#        :rtype: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperation
#        '''
#        oper = BinaryOperation(operator=operator)
#        if len(children) != 2:
#            raise IRGenerationError(f"Attempting to create a BinaryOperation with "
#                                    f"wrong number of children. Was provided "
#                                    f"{len(children)} but can only accept 2.")
#        if not isinstance(children[0], DataNode):
#            raise IRGenerationError(f"Attempting to create a BinaryOperation "
#                                    f"but first provided child is {type(children[0])} "
#                                    f"instead of a DataNode.")
#        if not isinstance(children[1], DataNode):
#            raise IRGenerationError(f"Attempting to create a BinaryOperation "
#                                    f"but second provided child is {type(children[1])} "
#                                    f"instead of a DataNode.")
#
#        oper.addchild(children[0])
#        oper.addchild(children[1])
#        return oper
#
#    @staticmethod
#    def _validate_child(position: int, child: Node) -> bool:
#        '''
#        Determines whether a given child and index are valid for this node.
#
#        :param int position: the position to be validated.
#        :param child: a child to be validated.
#        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
#
#        :return: whether the given child and position are valid for this node.
#        :rtype: bool
#        '''
#        if (position in (0, 1)) and isinstance(child, DataNode):
#            return True
#        return False
#
#    def node_str(self, colour=True) -> str:
#        '''
#        :returns: a string representation of this node.
#        :rtype: str
#        '''
#        return self.coloured_name(colour) + f"[{self.operator}: ({self.children[0]}, {self.children[1]})]"

