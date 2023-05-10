'''
This module contains the StructureReference class.
'''

from __future__ import annotations

from typing import Union

from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

class StructureReference(Reference):
    '''
    Contains a StructureReference in the IR tree.

    :param symbol: The structure type symbol that this reference points to.
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`
    :param member: The structure access member used in this reference.
    :type member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
    '''
    # pylint: disable=undefined-variable
    def __init__(self, symbol: StructureSymbol, member: Union[None,Member]=None) -> None:
        super().__init__(symbol=symbol)
        self.symbol = symbol
        self.member = member

    @staticmethod
    def _validate_child(position, child):
        '''
        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`psyclone.psyir.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool

        '''
        if position == 0:
            return isinstance(child, Member)
        return False

    @property
    def symbol(self) -> StructureSymbol:
        return super().symbol

    @symbol.setter
    def symbol(self, symbol: StructureSymbol) -> None:
        '''
        Sets the symbol this StructureReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: \
                :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, StructureSymbol):
            raise TypeError("Attempted to make a StructureReference to a non-StructureSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol

    @property
    def member(self) -> Union[None,Member]:
        '''
        :returns: the member of the structure accessed in this reference.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member` or None
        '''
        return self.children[0]

    @member.setter
    def member(self, member: Member) -> None:
        '''
        Sets the member this StructureReference accesses.

        :param member: the member of the structure to be accessed in this reference.
        :type member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`

        :raises TypeError: if the member is not the correct type.
        '''
        if member is not None and not isinstance(member, Member):
            raise TypeError("Attempted to make a StructureReference with a non-Member access. "
                            f"Got {type(member)} as input.")
        if len(self.children) > 0:
            self.children[0] = member
        else:
            self.addchild(member)

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"StructureReference[{self.symbol.name}: {self.member.node_str()}]"
        return nodestr
