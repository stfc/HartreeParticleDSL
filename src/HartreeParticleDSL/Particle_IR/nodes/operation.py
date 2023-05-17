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

import psyclone.psyir.nodes.operation as psyop

class Operation(psyop.Operation, metaclass=ABCMeta):
    '''
    Base class for Operations in Particle IR.

    :param operator: the operator describing this operation. This is specialised \
                     in subclasses of Operation.
    :type operator: object

    :raises TypeError: if the supplied operator is not a valid type for this operation.
    '''
    Operator = object
    _text_name = "Operation"

    def __init__(self, operator: object):
        super().__init__(operator=operator)

        if not isinstance(operator, self.Operator):
            raise TypeError(f"{type(self).__name__} expects an operator of type "
                            f"{self.Operator} but received type "
                            f"{type(operator)}.")
        self._operator = operator
        self._argument_names = []

    def __str__(self):
        return self.node_str(False)

    @property
    def operator(self) -> object:
        '''
        :returns: the operator for this Operation.
        :rtype: Any
        '''
        return self._operator

class BinaryOperation(Operation):
    '''
    Class representing a BinaryOperation. Has exactly two children to represent
    the two sides of the operation.
    '''
    _text_name = "BinaryOperation"
    class BinaryOp(Enum):
        '''Enumeration of Binary Operations supported in Particle DSL.'''
        # Arithmetic operations
        ADDITION=1
        SUBTRACTION=2
        MULTIPLY=3
        DIVISION=4
        # Comparison operations
        LESS_THAN=5
        LESS_THAN_EQUAL=6
        GREATER_THAN=7
        GREATER_THAN_EQUAL=8
        EQUALITY=9
        #Logical operations
        LOG_AND=10
        LOG_OR=11

    Operator = BinaryOp

    @staticmethod
    def create(operator: BinaryOp, children: List[DataNode]) -> BinaryOperation:
        '''
        Create a binary operation representing `children[0] operation children[1]`.

        :param operator: The binary operation represented by the created node.
        :type operator: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperation.BinaryOp`
        :param children: The children represented by this BinaryOperation.
        :type children: list of py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`

        :raises IRGenerationError: if the children input doesn't contain exactly 2 children \
                                   DataNodes.

        :returns: The new binary operation node.
        :rtype: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperation
        '''
        oper = BinaryOperation(operator=operator)
        if len(children) != 2:
            raise IRGenerationError(f"Attempting to create a BinaryOperation with "
                                    f"wrong number of children. Was provided "
                                    f"{len(children)} but can only accept 2.")
        if not isinstance(children[0], DataNode):
            raise IRGenerationError(f"Attempting to create a BinaryOperation "
                                    f"but first provided child is {type(children[0])} "
                                    f"instead of a DataNode.")
        if not isinstance(children[1], DataNode):
            raise IRGenerationError(f"Attempting to create a BinaryOperation "
                                    f"but second provided child is {type(children[1])} "
                                    f"instead of a DataNode.")

        oper.addchild(children[0])
        oper.addchild(children[1])
        return oper

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if (position in (0, 1)) and isinstance(child, DataNode):
            return True
        return False

    def node_str(self, colour=True) -> str:
        '''
        :returns: a string representation of this node.
        :rtype: str
        '''
        return self.coloured_name(colour) + f"[{self.operator}: ({self.children[0]}, {self.children[1]})]"

class UnaryOperation(Operation):
    '''
    Class representing a UnaryOperation. Has exactly one child.
    '''
    _text_name = "UnaryOperation"
    class UnaryOp(Enum):
        '''Enumeration of Unary Operations supported in Particle DSL.'''
        UNARYSUB=1
        LOG_NOT=2

    Operator = UnaryOp

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if position == 0 and isinstance(child, DataNode):
            return True
        return False

    @staticmethod
    def create(operator: UnaryOp, child: DataNode) -> UnaryOperation:
        '''
        Create a unary operation representing `operation children[0]`.

        :param operator: The unary operation represented by the created node.
        :type operator: py:class:`HartreeParticleDSL.Particle_IR.nodes.UnaryOpeartion.UnaryOp`
        :param children: The child of this UnaryOperation.
        :type children: py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`

        :raises IRGenerationError: if the child is not a DataNode.

        :returns: The new unary operation node.
        :rtype: py:class:`HartreeParticleDSL.Particle_IR.nodes.UnaryOperation
        '''
        oper = UnaryOperation(operator=operator)
        if not isinstance(child, DataNode):
            raise IRGenerationError(f"Attempting to create a UnaryOperation "
                                    f"but provided child is {type(child)} "
                                    f"instead of a DataNode.")

        oper.addchild(child)
        return oper

    def node_str(self, colour=True) -> str:
        '''
        :returns: a string representation of this node.
        :rtype: str
        '''
        return self.coloured_name(colour) + f"[{self.operator}: {self.children[0]}]"

# pylint: disable=fixme
#TODO Do we need an N-ary class for general operations?
