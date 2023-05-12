'''
This module contains the Member class.
'''

from HartreeParticleDSL.Particle_IR.nodes.node import Node

import psyclone.psyir.nodes.member as psyMem

class Member(Node, psyMem.Member):
    '''
    Node representing a member of a structure. This is the leaf member node.

    :param str name: the name of the member of the structure that is \
                            being referenced.
    '''
    _text_name = "Member"

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self._name = name

    def node_str(self) -> str:
        '''
        :returns: a string representation of this node.
        :rtype: str
        '''
        return f"Member[{self.name}]"
