'''
This module contains the ConfigReference class.
'''

from __future__ import annotations
from HartreeParticleDSL.Particle_IR.nodes.structure_reference import StructureReference
from HartreeParticleDSL.Particle_IR.nodes.member import Member
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol

class ConfigReference(StructureReference):
    '''
    Contains a ConfigReference in the IR tree.

    :param symbol: The structure type symbol that this reference points to.
    :type symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`
    :param member: The structure access member used in this reference.
    :type member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
    '''
    # pylint: disable=undefined-variable
    def __init__(self, symbol: StructureSymbol, member: Member) -> None:
        super().__init__(symbol, member)

    @property
    def symbol(self) -> ScalarTypeSymbol:
        return super().symbol

    @symbol.setter
    def symbol(self, symbol: StructureSymbol) -> None:
        '''
        Sets the symbol this ConfigReference refers to.

        :param symbol: The symbol to make this symbol refer to.
        :type symbol: \
                :py:class:`HartreeParticleDSL.Particle_IR.symbols.structuresymbol.StructureSymbol`

        :raises TypeError: if the symbol is not the correct type.
        '''
        if not isinstance(symbol, StructureSymbol):
            raise TypeError("Attempted to make a ConfigReference to a non-StructureSymbol. "
                            f"Got {type(symbol)} as input.")
        self._symbol = symbol

    @property
    def member(self) -> Member:
        '''
        :returns: the member of the structure accessed in this reference.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`
        '''
        return self._member

    @member.setter
    def member(self, member: Member) -> None:
        '''
        Sets the member this ConfigReference accesses.

        :param member: the member of the structure to be accessed in this reference.
        :type member: :py:class:`HartreeParticleDSL.Particle_IR.nodes.member.Member`

        :raises TypeError: if the member is not the correct type.
        '''
        if not isinstance(member, Member):
            raise TypeError("Attempted to make a ConfigReference with a non-Member access. "
                            f"Got {type(member)} as input.")
        self._member = member

    def node_str(self) -> str:
        '''
        :returns: a text description of this node
        :rtype: str
        '''
        nodestr = f"ConfigReference[{self.symbol.name}: {self.member.node_str()}]"
        return nodestr
