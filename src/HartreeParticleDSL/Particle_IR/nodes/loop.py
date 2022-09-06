from __future__ import annotations
from typing import List

from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.nodes.node import DataNode, Node
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol

class Loop(Statement):
    '''
    Class to represent a Loop in the Particle IR tree.

    :param variable: The loop variable used in this loop.
    :type variable: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`

    :raises IRGenerationError: if the variable input is invalid.
    '''
    def __init__(self, variable: Symbol) -> None:
        super().__init__()
        if not isinstance(variable, Symbol):
            raise TypeError(f"Loop variable needs to be a symbol but "
                                    f"got {type(variable)}.")
        self._variable = variable

    @staticmethod
    def create(variable: Symbol, start: DataNode, stop: DataNode, step: DataNode, children: List[DataNode]) -> Loop:
        '''
        Creates and returns a Loop node.

        :param variable: The loop variable used in this loop.
        :type variable: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
        :param start: The start node used in this loop.
        :type start: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        :param stop: The stop node used in this loop.
        :type stop: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        :param step: The step node used in this loop.
        :type step: :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        :param children: The children nodes used in this loop.
        :type children: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`

        :raises IRGenerationError: if any of the inputs are invalid.

        :returns: a new Loop node representing the inputs.
        :rtype: :py:class::py:class:`HartreeParticleDSL.Particle_IR.nodes.loop.Loop`
        '''
        if not isinstance(variable, Symbol):
            raise IRGenerationError(f"Loop variable needs to be a symbol but "
                                    f"got {type(variable)}.")
        if not isinstance(start, DataNode):
            raise IRGenerationError(f"Loop start needs to be a DataNode but "
                                    f"got {type(start)}.")
        if not isinstance(stop, DataNode):
            raise IRGenerationError(f"Loop stop needs to be a DataNode but "
                                    f"got {type(stop)}.")
        if not isinstance(step, DataNode):
            raise IRGenerationError(f"Loop step needs to be a DataNode but "
                                    f"got {type(step)}.")

        loop = Loop(variable=variable)
        body = Body()
        loop.addchild(start)
        loop.addchild(stop)
        loop.addchild(step)
        loop.addchild(body)
        for child in children:
            body.addchild(child)

        return loop
        

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.

        Loops are "range" style loops, i.e. i in range(0, 100, 10)
        First 3 children are therefore the DataNodes for this range.
        The fourth child is the Body node for this loop

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if (position >= 0 and position <= 2) and isinstance(child, DataNode):
            return True
        if position == 3 and isinstance(child, Body):
            return True
        return False

    def _check_completeness(self):
        '''
        Checks whether this loop is fully formed, i.e. has 4 children, 
        start, stop, step and body.

        :raises IRGeneration: If this node is not fully formed as expected.
        '''
        if len(self.children) != 4:
            raise IRGenerationError("Loop node should have 4 children "
                                    f"but found {len(self.children)}")

    @property
    def start_expr(self) -> DataNode:
        '''
        :returns: the DataNode representing the start in this Loop.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
        '''
        self._check_completeness()
        return self.children[0]

    @property
    def stop_expr(self) -> DataNode:
        '''
        :returns: the DataNode representing the stop in this Loop.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
        '''
        self._check_completeness()
        return self.children[1]

    @property
    def step_expr(self) -> DataNode:
        '''
        :returns: the DataNode representing the step in this Loop.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
        '''
        self._check_completeness()
        return self.children[2]

    def node_str(self) -> str:
        '''
        :returns: a text description of this loop
        :rtype: str
        '''
        my_str = (f"Loop[ ({self.start_expr.node_str()}, {self.stop_expr.node_str()} "
                 f"{self.step_expr.node_str()}): {self.children[3].node_str()}\n]")
        return my_str

