from abc import ABCMeta
import Enum

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.node import DataNode, Node

class Operation(DataNode, metaclass=ABCMeta):
    '''
    Base class for Operations in Particle IR.

    :param operator: the operator describing this operation. This is specialised \
                     in subclasses of Operation.
    :type operator: object

    :raises TypeError: if the supplied operator is not a valid type for this operation.
    '''
    Operator = object

    def __init__(self, operator: object):
        super().__init__()

        if not isinstance(operator, self.Operator):
            raise TypeError(f"{type(self)} expects an operator of type "
                            f"{self.Operator} but received type "
                            f"{type(operator)}.")
        self._operator = operator

    @property
    def operator(self) -> object:
        return self._operator

class BinaryOperator(Operation):
    '''
    Class representing a BinaryOperation. Has exactly two children to represent
    the two sides of the operation.
    '''
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
    def create(operator: BinaryOp, children: List[DataNode]) -> BinaryOperator:
        '''
        Create a binary operation representing `children[0] operation children[1]`.

        :param operator: The binary operation represented by the created node.
        :type operator: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperator.BinaryOp`
        :param children: The children represented by this BinaryOperation.
        :type children: list of py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`

        :raises IRGenerationError: if the children input doesn't contain exactly 2 children DataNodes.
    
        :returns: The new binary operation node.
        :rtype: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperator
        '''
        op = BinaryOperator(operator=operator)
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

        op.addchild(children[0])
        op.addchild(children[1])
        return op

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
        if (position == 0 or position == 1) and isinstance(child, DataNode):
            return True
        return False

class UnaryOperation
    '''
    Class representing a UnaryOperation. Has exactly one child.
    '''
    class UnaryOp(Enum):
        '''Enumeration of Unary Operations supported in Particle DSL.'''
        UNARYSUB=1
        LOGICAL_NOT=2


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
    def create(operator: UnaryOp, child: DataNode) -> BinaryOperator:
        '''
        Create a binary operation representing `children[0] operation children[1]`.

        :param operator: The binary operation represented by the created node.
        :type operator: py:class:`HartreeParticleDSL.Particle_IR.nodes.UnaryOpeartion.UnaryOp`
        :param children: The child of this UnaryOperation.
        :type children: py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`

        :raises IRGenerationError: if the child is not a DataNode.
    
        :returns: The new binary operation node.
        :rtype: py:class:`HartreeParticleDSL.Particle_IR.nodes.BinaryOperator
        '''
        op = BinaryOperator(operator=operator)
        if not isinstance(child, DataNode):
            raise IRGenerationError(f"Attempting to create a UnaryOperation "
                                    f"but provided child is {type(child)} "
                                    f"instead of a DataNode.")

        op.addchild(child)
        return op

#TODO Do we need an N-ary class for general operations?
