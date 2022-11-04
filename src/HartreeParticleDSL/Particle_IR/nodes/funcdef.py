'''
This module contains the FuncDef class.
'''

from __future__ import annotations

from typing import List, Union
from HartreeParticleDSL.Particle_IR.nodes.node import Node, DataNode
from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable

class FuncDef(Node):
    '''
    Class representing a general (non-kernel) Function Definition.

    :param str name: The name of this function definition.
    :param children: List of Nodes to be contained in this Kern region.
    :type children: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.Node. \
            or None.
    '''
    # pylint: disable=undefined-variable

    def __init__(self, name: str,  children: Union[None,List[Node]]=None) -> None:
        super().__init__(children=children)

        self.name = name
        # If no children were specified, we need to add the Body of this Kern
        if not children or len(children)==0:
            self.addchild(Body())

        # Specialisations of this class will need argument lists
        self._arguments = []

        # Add a symbol table
        self._symbol_table = SymbolTable(kern=self)

    @property
    def symbol_table(self) -> SymbolTable:
        '''
        :returns: The symbol table for this FuncDef.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol_table.SymbolTable`
        '''
        return self._symbol_table

    @property
    def body(self) -> Body:
        '''
        :returns: The body of this FuncDef.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.body.Body`
        '''
        return self.children[0]

    @property
    def arguments(self) -> List[DataNode]:
        '''
        :returns: The argument list of this FuncDef.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[DataNode]) -> None:
        '''
        Sets the arguments of this Function Definition.

        :param arguments: list of arguments for this pairwise kernel.
        :type arguments: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`

        :raises TypeError: If any provided argument is not a DataNode.
        '''
        self._arguments = []
        for arg in arguments:
            if not isinstance(arg, DataNode):
                raise TypeError("Each argument must be a DataNode, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
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

    def node_str(self) -> str:
        '''
        :returns: a text description of this FuncDef
        :rtype: str
        '''
        arg_strings = []
        for arg in self._arguments:
            arg_strings.append(arg.node_str())
        arg_string = ", ".join(arg_strings)
        nodestr = type(self).__name__ + f"[{arg_string}: {self.children[0].node_str()}]"
        return nodestr

    @property
    def name(self) -> str:
        '''
        :returns: The name of this perpart kernel.
        :rtype: str
        '''
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        '''
        Sets the name of this PerPartKernel

        :param str name: The name of this perpart kernel.

        :raises TypeError: If the provided name is not a str.
        '''
        if not isinstance(name, str):
            raise TypeError("Expected FuncDef name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @staticmethod
    def create(name: str, arguments: List[DataNode], body: List[Statement]) -> FuncDef:
        '''
        Creates a FuncDef containing the input arguments and kernel body.

        :param str name: The name of this FuncDef.
        :param arguments: The list of arguments to use for this FuncDef.
        :type arguments: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.node.DataNode`
        :param body: The list of Statements that make up this FuncDef
        :type body: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

        :returns: The new FuncDef created
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.funcdef.FuncDef
        '''
        kernel = FuncDef(name)
        kernel.arguments = arguments

        for node in body:
            kernel.body.addchild(node)

        return kernel
