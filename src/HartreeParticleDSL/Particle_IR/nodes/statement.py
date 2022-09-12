from abc import ABCMeta

from HartreeParticleDSL.Particle_IR.nodes.node import Node

class Statement(Node, metaclass=ABCMeta):
    '''
    Abstract node representing a Statement.
    '''

class EmptyStatement(Statement):
    '''
    Node used to representing a part of the AST tree that
    is not included in the Particle_IR output, e.g.
    a create variable call without value which can be done
    through the symbol table output.
    '''

