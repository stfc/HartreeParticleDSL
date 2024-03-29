'''
This module contains the Call class.
'''

from __future__ import annotations
from typing import List

from HartreeParticleDSL.Particle_IR.nodes.node import DataNode
from HartreeParticleDSL.Particle_IR.nodes.statement import Statement

class Call(Statement, DataNode):
    '''
    Contains a generic call in the tree.

    :param str name: the function name called by this Call.
    '''
    # pylint: disable=undefined-variable

    def __init__(self, func_name: str) -> None:
        super().__init__()
        self._func_name = func_name

    @property
    def func_name(self) -> str:
        '''
        :returns: The function name of this call.
        :rtype: str
        '''
        return self._func_name

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
        if isinstance(child, DataNode):
            return True
        return False

    @staticmethod
    def create(func_name: str, arguments: List[DataNode]) -> Call:
        '''
        Creates a Call node for the supplied function with the supplied
        arguments.

        :param func_name: The name of the function call represented by this \
                          Call.
        :type func_name: str
        :param arguments: The list of arguments of this Call.
        :type arguments: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`

        :raises TypeError: If the func_name argument is not a str.

        :returns: The new Call node to represent the inputs.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.call.Call`
        '''
        if not isinstance(func_name, str):
            raise TypeError(f"Expected func_name to be a str but got "
                            f"{type(func_name)}.")

        rval = Call(func_name=func_name)
        for arg in arguments:
            rval.addchild(arg)
        return rval

    def node_str(self) -> str:
        '''
        :returns: a text description of this node.
        :rtype: str
        '''
        rstring = f"Call[{self.func_name}: ("
        arg_strs = []
        for arg in self.children:
            arg_strs.append(str(arg))
        rstring = rstring + ", ".join(arg_strs) + ")]"

        return rstring
