'''
This module contains the abstract Reference class.
'''

from abc import ABCMeta
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol

import psyclone.psyir.nodes.reference as psyRef

#class Reference(psyRef.Reference, metaclass=ABCMeta):
#    '''
#    Contains a Reference to a variable in the ParticleIR tree.
#    '''
#    _text_name = "Reference"
