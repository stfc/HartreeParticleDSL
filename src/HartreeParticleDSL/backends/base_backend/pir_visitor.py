from __future__ import annotations
import inspect

from HartreeParticleDSL.Particle_IR.nodes.node import Node

class PIR_Visitor():
    '''
    Generic PIR visitor.

    :param str indent_string: The indentation to use in this Visitor. Default \
                              "    ".
    :param int initial_indent_depth: Specifies the initial indentation, default \
                                     is 0.

    :raises TypeError: if any of the inputs are the wrong type.
    '''

    def __init__(self, indent_string: str="    ", initial_indent_depth: int=0) -> None:

        if not isinstance(indent_string, str):
            raise TypeError("Expected indent_string to be a string but got "
                            f"{type(indent_string)}")
        if not isinstance(initial_indent_depth, int):
            raise TypeError("Expected initial_indent_depth to be an int but "
                            f"got {type(initial_indent_depth)}")

        self._indent = indent_string
        self._depth = initial_indent_depth

    @property
    def _nindent(self) -> str:
        '''
        :returns: the current indentation string.
        :rtype: str

        '''
        return self._depth * self._indent

    def indent(self) -> None:
        '''
        Increase the current indentation.
        '''
        self._depth = self._depth + 1

    def dedent(self) -> None:
        '''
        Decrease the current indentation.
        '''
        self._depth = self._depth - 1

    def __call__(self, node: Node) -> str:
        '''
        This is called when an instance of this visitor class is called as a
        function. It walks through the tree and visits them according to the
        visitor functions defined here.

        :param node: A Particle IR Node.
        :type node: py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :returns: Code representaiton of this Particle IR tree.
        :rtype: str

        :raises TypeError: if provided node is not a Particle IR Node.
        '''
        if not isinstance(node, Node):
            raise TypeError("Expected a Particle IR Node as input but got "
                            f"{type(node)}.")

        return self._visit(node)

    def _visit(self, node: Node) -> str:
        '''
        Implements the Particle IR visitors. This is implemented using the
        class hierarchy names as the candidate method names, starting from
        the class name of the object and then recursing to its parents until
        all parents are exhausted. Names are converted to lower case but
        otherwise not modified.

        This code is heavily adapted from PSyclone.

        :param node: A Particle IR Node.
        :type node: py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :returns: Code representaiton of this Particle IR tree.
        '''

        # Make a list of the node's ancestor classes (including
        # itself) in method resolution order (mro), apart from the
        # base "object" class.
        possible_method_names = ["visit_"+curr_class.__name__.lower()+"_node"
                                 for curr_class in inspect.getmro(type(node))]
        possible_method_names.remove("visit_object_node")

        # Try to call methods using the class names in the order of
        # the class hierarchy (starting from the current class name).
        for method_name in possible_method_names:
            try:
                # pylint: disable=eval-used
                node_result = eval(f"self.{method_name}(node)")

                return node_result

            except AttributeError as excinfo:
                if f"attribute '{method_name}'" in str(excinfo):
                    # This attribute error is because the method that
                    # was tried does not exist.
                    pass
                else:
                    # The method does exist but it has raised an
                    # attribute error so re-raise it here.
                    raise AttributeError(excinfo) from excinfo

        # We haven't found a handler for this node.
        raise VisitorError(
            f"Unsupported node '{type(node).__name__}' found: method names "
            f"attempted were {possible_method_names}.")

    # Implement visitors
    def visit_emptystatement_node(self, node: EmptyStatement) -> str:
        '''
        Upon finding an Empty Statement node, we just return an empty string
        '''
        return ""
