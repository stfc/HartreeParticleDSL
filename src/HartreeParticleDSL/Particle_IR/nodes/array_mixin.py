'''
This module contains the abstract ArrayMixin class.
'''

from __future__ import annotations
from abc import ABCMeta
from typing import List

from HartreeParticleDSL.Particle_IR.nodes.node import DataNode

class ArrayMixin(metaclass=ABCMeta):
    '''
    Abstract class used to handle all common functionality for Arrays.
    This should always preceed other inheritence.
    '''
    # pylint: disable=undefined-variable

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
        # pylint: disable=unused-argument
        return isinstance(child, DataNode)

    @property
    def is_array(self) -> bool:
        '''
        :returns: whether this instance indicated an array access.
        :rtype: bool
        '''
        return True

    @property
    def indices(self) -> List[DataNode]:
        '''
        Returns the list of nodes that make up the index accesses for this
        array.

        :returns: the nodes representing the array-index expressions.
        :rtype: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
        '''
        return self.children
