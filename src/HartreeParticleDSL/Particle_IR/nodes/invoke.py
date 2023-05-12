'''
This module contains the Invoke class.
'''

from __future__ import annotations
from typing import List, Union

from HartreeParticleDSL.Particle_IR.nodes.statement import Statement
from HartreeParticleDSL.Particle_IR.nodes.literal import Literal

class Invoke(Statement):
    '''
    Dictates an Invoke call in Particle_IR.

    An Invoke is used for any Kern objects that are not Main.
    '''
    _text_name = "Invoke"
    # pylint: disable=undefined-variable

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
        if isinstance(child, Literal):
            return True
        return False

    @staticmethod
    def create(kernels: Union[Literal, List[Literal]]) -> Invoke:
        '''
        Creates an Invoke node for the kernels represented by the literals
        passed to this function.

        :param kernels: kernel or list of kernels used in this invoke.
        :type kernels: (list of) :py:class:`HartreeParticleDSL.Particle_IR.nodes.literal.Literal`

        :returns: the Invoke object representing the supplied kernel calls.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.invoke.Invoke`
        '''
        rval = Invoke()
        if isinstance(kernels, list):
            for child in kernels:
                rval.addchild(child)
        else:
            rval.addchild(kernels)

        return rval
