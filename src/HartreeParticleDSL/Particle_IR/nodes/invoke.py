from typing import List, Union

from HartreeParticleDSL.Particle_IR.nodes.statement import Statement
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.kernels import Main

def Invoke(Statement):
    '''
    Dictates an Invoke call in Particle_IR.

    An Invoke is used for any Kern objects that are not Main.
    '''

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
        if isinstance(child, Kern) and not isinstance(child, Main):
            return True
        return False

    @staticmethod
    def create(kernels: Union[Kern, List[Kern]]) -> Invoke:
        '''
        Creates an Invoke node for the Kern or list of Kern objects
        passed to this function.

        :param kernels: kernel or list of kernels used in this invoke.
        :type kernels: (list of) :py:class:`HartreeParticleDSL.Particle_IR.nodes.kern.Kern`

        :returns: the Invoke object representing the supplied kernels.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.invoke.Invoke`
        '''
        rval = Invoke()
        if isinstance(kernels, list):
            for child in kernels:
                rval.addchild(child)
        else:
            rval.addchild(kernels)

        return rval


