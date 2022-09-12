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
    def __init__(self, symbol: StructureSymbol, member: Union[NoneType,Member]=None) -> None:
        super().__init__()
        self.symbol = symbol
        self.member = member

    @property
    def symbol(self) -> StructureSymbol:
        return super().symbol

    @symbol.setter
    def symbol(self, symbol: StructureSymbol) -> None:
        '''
        Sets the symbol this StructureReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, StructureSymbol):
            raise TypeError("Attempted to make a StructureReference to a non-StructureSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol

    @property
    def member(self) -> Union[NoneType,Member]:
        '''
        :returns: the member of the structure accessed in this reference.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member` or None
        '''
        return self._member

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
        self._member = member

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"StructureReference[{self.symbol.name}: {self.member.node_str()}]"
        return nodestr
