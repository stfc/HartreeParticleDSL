'''
This module contains the IfElseBlock class.
'''

from __future__ import annotations
from typing import List

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.node import Node
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement

class IfElseBlock(Statement):
    '''
    Class representing an If else block in the tree. The first child
    is the if condition, the second is the if body and the final is
    the else body.

    Backends can output their else if conditional if the child of
    the else body is another IfElseBlock, however they don't have to.
    '''

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.

        An IfElseBlock should have 1 Node child and 2 Body nodes.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if position == 0 and isinstance(child, Node):
            return True
        if (position in (1, 2)) and isinstance(child, Body):
            return True
        return False

    @staticmethod
    def create(condition: Node, ifbody: List[Node], elsebody: List[Node]) -> IfElseBlock:
        '''
        Create an IfElseBlock node representing the provided condition,
        ifbody and elsebody.

        :param condition: the if condition to evaluate for this IfElseBlock.
        :type condition: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node`
        :param ifbody: the nodes to execute if the condition evaluates to true.
        :type ifbody: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node`
        :param elsebody: the nodes to execute if the condition evaluates to false.
        :type elsebody: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node`

        :raises TypeError: If the condition is the wrong type.

        :returns: The IfElseBlock node representing this input.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.ifelse.IfElseBlock`
        '''
        if not isinstance(condition, Node):
            raise TypeError(f"The condition input needs to be a Node, but "
                            f"found {type(condition)}.")
        ifblock = IfElseBlock()
        body1 = Body()
        body2 = Body()
        for child in ifbody:
            body1.addchild(child)
        for child in elsebody:
            body2.addchild(child)
        ifblock.addchild(condition)
        ifblock.addchild(body1)
        ifblock.addchild(body2)

        return ifblock

    def _check_completeness(self) -> None:
        '''
        Checks this IFElseBlock is fully formed as expected.

        :raises IRGenerationError: if this node doesn't have the correct \
                                   number of children
        '''
        if len(self.children) != 3:
            raise IRGenerationError(f"IfElseBlock must have 3 children but "
                                    f"found only {len(self.children)}.")
        if not isinstance(self.children[0], Node):
            raise IRGenerationError(f"IfElseBlock first child must be a Node "
                                    f"but found {type(self.children[0])}.")
        if not isinstance(self.children[1], Body):
            raise IRGenerationError(f"IfElseBlock second child must be a Body "
                                    f"but found {type(self.children[1])}.")
        if not isinstance(self.children[2], Body):
            raise IRGenerationError(f"IfElseBlock third child must be a Body "
                                    f"but found {type(self.children[2])}.")

    @property
    def condition(self) -> Node:
        '''
        :returns: the condition node associated with this IfBlock.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node`
        '''
        self._check_completeness()
        return self.children[0]

    @property
    def ifbody(self) -> Body:
        '''
        :returns: the body assoiciated with the true path for this IfBlock.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.body.Body`
        '''
        self._check_completeness()
        return self.children[1]

    @property
    def elsebody(self) -> Body:
        '''
        :returns: the body assoiciated with the false path for this IfBlock.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.body.Body`
        '''
        self._check_completeness()
        return self.children[2]

    def is_else_if(self) -> bool:
        '''
        Computes if the else body is an elseif block. To be an elseif, the
        first child of the elsebody needs to be an Ifblock.

        :returns: true if the else is an elseif, otherwise false.
        :rtype: bool
        '''
        if ( len(self.elsebody.children) > 0 and
             isinstance(self.elsebody.children[0], IfElseBlock)):
            return True
        return False

    def node_str(self) -> str:
        '''
        :returns: a string representation of this node.
        :rtype: str
        '''
        val = f"IfElseBlock[{self.condition}:\n"
        val = val + f"{self.ifbody}\n"
        val = val + f"{self.elsebody}\n"
        val = val + "]"
        return val
