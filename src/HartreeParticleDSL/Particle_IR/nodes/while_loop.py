from __future__ import annotations
from typing import List

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.node import DataNode, Node
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol

class While(Statement):
    '''
    Class to represent a While loop in the Particle IR tree.
    '''
    def __init__(self) -> None:
        super().__init__()

    @staticmethod
    def create(condition: DataNode, children: List[DataNode]) -> While:
        '''
        Creates and returns a While node.

        :param condition: The condition used for this while loop.
        :type condition: :py:class:`HartreeParticleDSL.Particle_IR.symbols.node.DataNode`
        :param children: The children nodes used in this loop.
        :type children: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`

        :raises IRGenerationError: if any of the inputs are invalid.

        :returns: a new While node representing the inputs.
        :rtype: :py:class::py:class:`HartreeParticleDSL.Particle_IR.nodes.while.While`
        '''
        if not isinstance(condition, DataNode):
            raise IRGenerationError(f"While condition needs to be a DataNode but "
                                    f"got {type(condition)}.")

        while_loop = While()
        body = Body()
        while_loop.addchild(condition)
        while_loop.addchild(body)
        for child in children:
            body.addchild(child)

        return while_loop
        

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.

        While loops need a DataNode as first child and Body as second child.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if (position == 0) and isinstance(child, DataNode):
            return True
        if position == 1 and isinstance(child, Body):
            return True
        return False

    def _check_completeness(self):
        '''
        Checks whether this while loop is fully formed, i.e. has 2 children, 
        condition and body.

        :raises IRGeneration: If this node is not fully formed as expected.
        '''
        if len(self.children) != 2:
            raise IRGenerationError("While node should have 2 children "
                                    f"but found {len(self.children)}")

    @property
    def condition(self) -> Statement:
        '''
        :returns: the condition of this While.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.statement.Statement`
        '''
        self._check_completeness()
        return self.children[0]


    @property
    def body(self) -> Body:
        '''
        :returns: the Body representing the body of this Loop.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Body`
        '''
        self._check_completeness()
        return self.children[1]

    def node_str(self) -> str:
        '''
        :returns: a text description of this loop
        :rtype: str
        '''
        my_str = (f"While[{self.children[0].node_str()}: {self.children[1].node_str()}\n]")
        return my_str

