from __future__ import annotations

from HartreeParticleDSL.Particle_IR.nodes.member import Member

class StructureMember(Member):
    '''
    Node representing a member of a structure that is also a structure.

    :param str name: The name of this structure member.
    :param sub_member: The member this structure contains as an access.
    :type sub_member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
    '''
    def __init__(self, name: str, sub_member: Member) -> None:
        super().__init__(name)
        self.addchild(sub_member)

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
        return (position == 0 and isinstance(child, Member))

    @property
    def member(self) -> Member:
        '''
        :returns: the member of the structure being accessed.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
        '''
        return self.children[0]

    @member.setter
    def member(self, sub_member: Member) -> None:
        '''
        Sets the member attribute of this StructureMember.

        :param sub_member: The member this structure contains as an access.
        :type sub_member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
        '''
        self.children[0] = sub_member

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"StructureMember[{self.name}: {self.member.node_str()}]"
        return nodestr
