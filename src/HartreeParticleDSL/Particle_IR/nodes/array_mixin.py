'''
This module contains the abstract ArrayMixin class.
'''

from __future__ import annotations
from abc import ABCMeta
from typing import List

import psyclone.psyir.nodes.array_mixin as psyArrmix

class ArrayMixin(psyArrmix.ArrayMixin):
    _text_name = "ArrayMixin"

    @property
    def indices(self) -> List[DataNode]:
        '''
        Returns the list of nodes that make up the index accesses for this
        array.

        :returns: the nodes representing the array-index expressions.
        :rtype: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.DataNode`
        '''
        return self.children
