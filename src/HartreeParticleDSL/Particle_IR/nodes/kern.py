from abc import ABCMeta
from HartreeParticleDSL.Particle_IR.nodes.node import Node
from HartreeParticleDSL.Particle_IR.nodes.body import Body

class Kern(Node, metaclass=ABCMeta):
    '''
    Abstract class representing a general Kernel object.
    
    :param children: List of Nodes to be contained in this Kern region.
    :type children: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node.
    '''

    def __init__(self, children=None: List[Node]) -> None:
        super.__init__(children=children)

        # If no children were specified, we need to add the Body of this Kern
        if len(children=0):
            self.addchild(Body())

        # Specialisations of this class will need argument lists
        self._arguments = []

    @classmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given position and node are valid as a child
        of this node.

        Kern nodes can contain a single Statement node as its child.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        if position == 0 and isinstance(child, Body):
            return True
        return False
