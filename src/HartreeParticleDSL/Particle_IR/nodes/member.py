'''
This module contains the Member class.
'''

from HartreeParticleDSL.Particle_IR.nodes.node import Node

class Member(Node):
    '''
    Node representing a member of a structure. This is the leaf member node.

    :param str name: the name of the member of the structure that is \
                            being referenced.
    '''

    def __init__(self, name: str) -> None:
        super().__init__()
        self._name = name

    @property
    def name(self) -> str:
        '''
        :returns: the name of this member.
        :rtype: str
        '''
        return self._name

    @property
    def is_array(self) -> bool:
        '''
        :returns: if this member is an array.
        :rtype: bool
        '''
        return False

    def node_str(self) -> str:
        '''
        :returns: a string representation of this node.
        :rtype: str
        '''
        return f"Member[{self.name}]"
