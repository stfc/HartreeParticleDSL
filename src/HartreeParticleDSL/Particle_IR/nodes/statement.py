'''
This module contains the abstract Statement class, as well as some
child classes, EmptyStatement, Break and Return.
'''

from abc import ABCMeta

from HartreeParticleDSL.Particle_IR.nodes.node import Node, DataNode

import psyclone.psyir.nodes.statement as psyStatement

class Statement(Node, psyStatement.Statement, metaclass=ABCMeta):
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
    _text_name = "EmptyStatement"

class Break(Statement):
    '''
    Node used to represent a break statement.
    '''
    _text_name = "Break"

class Return(Statement):
    '''
    Node used to represent a Return statement.
    '''
    _text_name = "Return"

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
        if position == 0 and isinstance(child, DataNode):
            return True
        return False
