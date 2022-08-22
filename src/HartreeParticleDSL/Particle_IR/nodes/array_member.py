from HartreeParticleDSL.Particle_IR.nodes.node import DataNode
from HartreeParticleDSL.Particle_IR.nodes.array_mixin import ArrayMixin
from HartreeParticleDSL.Particle_IR.nodes.member import Member

class ArrayMember(ArrayMixin, Member):
    '''
    Node representing an access to the elemnts of an array inside a structure.
    Must have one or more children which give the array-index expressions
    for the array access.

    :param str name: the name of the member of the structure being referenced.
    :param indices: the list of indices that make up the array acces.
    :type indices: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
    '''

    def __init__(self, name, indices):
        super().__init__(name=name)

        for index in indices:
            self.addchild(index)