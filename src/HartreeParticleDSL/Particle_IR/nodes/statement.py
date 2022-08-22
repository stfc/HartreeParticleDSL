from abc import ABCMeta

from HartreeParticleDSL.Particle_IR.nodes.node import Node

class Statement(Node, metaclass=ABCMeta):
    '''
    Abstract node representing a Statement.
    '''

