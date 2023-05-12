'''
This module contains the ChildrenList class and the base Node class.
'''

from __future__ import annotations
from typing import Union, Tuple, List, Any
from abc import ABCMeta
from collections.abc import Callable
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

import psyclone.psyir.nodes.node as psynode
import psyclone.psyir.nodes.datanode as psydatanode

class ChildrenList(psynode.ChildrenList):
    pass

class Node(psynode.Node):

    def same_parent(self, node_2: Node) -> bool:
        '''
        :param node_2: The node to check if has the same parent.
        :type node_2: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :returns: True if node_2 has the same parent as this node, False otherwise.
        :rtype: bool
        '''
        if self.parent is None or node_2.parent is None:
            return False
        return self.parent is node_2.parent

    def validate_constraints(self) -> None:
        '''
        Validates this node in the context of the PIR tree.
        This checks constraints that can only be checked once the
        tree is complete.

        By default this routine does nothing, and must be overridden if
        required.
        '''

class DataNode(Node, psydatanode.DataNode, metaclass=ABCMeta):
    '''
    Abstract root node for any Node involving data accesses, e.g. Reference
    '''

