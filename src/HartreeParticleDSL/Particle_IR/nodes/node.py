from __future__ import annotations
from abc import ABCMeta
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
import typing
from typing import Union, Tuple, List
from collections.abc import Callable

class ChildrenList(list):
    '''
    Customized list to keep track of children nodes. It handles validation
    of inserted children according to the node it is created for.

    Since this is a subclass of the standard list, all operations are
    conserved and make use fo the validation.

    :param node: reference to the node where the list belongs.
    :type node: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
    :param validation_function: the callback function to use to validation children.
    :type validation_function: \
            function(int, :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`)
    '''

    def __init__(self, node: Node, validation_function: Callable[int, Node]):
        super().__init__()
        self._node_reference = node
        self._validation_function = validation_function

    def _validate_item(self, index: int, item: Any) -> None:
        '''
        Validates the provided index and item.

        :param int index: position where the item is inserted.
        :param item: object to be validated in the specified index.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :raises IRGenerationError: If the item cannot be a child at the index given.
        '''
        if not self._validation_function(index, item):
            raise IRGenerationError(f"Item '{type(item)}' can't be child {index}"
                                    f" of '{type(self._node_reference)}'.")

    def _check_is_orphan(self, item: Node) -> None:
        '''
        Checks that the provided item has no parent specified.

        :param item: object that needs to be validated.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :raises IRGenerationError: If the given item is not an orphan.
        '''
        if item.parent is not None:
            raise IRGenerationError(f"Item '{item}' can't be added as child of "
                f"'{self._node_reference}' because it is not an orphan. It already "
                f"has a '{item.parent}' as a parent.")

    def _set_parent_link(self, node: Node) -> None:
        '''
        Set parent connection of the given node to the node that this ChildrenList belongs to

        :param node: node for which the parent connection needs to be updated.
        :type node: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        '''
        node._parent = self._node_reference

    @staticmethod
    def _del_parent_link(node: Node) -> None:
        '''
        Delete parent connection of the given node.

        :param node: node for which the parent connection needs to be updated.
        :type node: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        node._parent = None

    def append(self, item: Node) -> None:
        '''
        Extends list append method to include children node validation.

        :param item: item to be appended to the list.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        self._validate_item(len(self), item)
        self._check_is_orphan(item)
        super().append(item)
        self._set_parent_link(item)

    def __setitem__(self, index: int, item: Node) -> None:
        '''
        Extends list __setitem__ method with children node validation.

        :param int index: position where to insert the item.
        :param item: item to be inserted to the list.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        '''
        self._validate_item(index, item)
        self._check_is_orphan(item)
        self._del_parent_link(self[index])
        super().__setitem__(index, item)
        self._set_parent_link(item)

    def insert(self, index: int, item: Node) -> None:
        '''
        Extends list insert method with children node validation.

        :param int index: position where to insert the item.
        :param item: item to be inserted to the list.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        positiveindex = index if index >= 0 else len(self) - index
        self._validate_item(positiveindex, item)
        self._check_is_orphan(item)
        # Check that all displaced items will still be in valid positions
        for position in range(positiveindex, len(self)):
            self._validate_item(position+1, self[position])
        super().insert(index, item)
        self._set_parent_link(item)

    def extend(self, items: list[Node]) -> None:
        '''
        Extends list extend method with children node validation.
        
        :param items: list of items to be inserted to the list.
        :type items: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        for index, item in enumerate(items):
            self._validate_item(len(self) + index, item)
            self._check_is_orphan(item)
        super().extend(items)
        for item in items:
            self._set_parent_link(item)

    def __delitem__(self, index: int) -> Node:
        '''
        Extends list __delitem__ method with children node validation.

        :param int index: position where to remove the item.

        '''
        positiveindex = index if index >= 0 else len(self) - index
        for position in range(positiveindex + 1, len(self)):
            self._validate_item(position-1, self[position])
        self._del_parent_link(self[index])
        super().__delitem__(index)

    def remove(self, item: Node) -> None:
        '''
        Extends list remove method with children node validation.

        :param item: item to be deleted from the list.
        :type item: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        for position in range(self.index(item)+1, len(self)):
            self._validate_item(position - 1, self[position])
        self._del_parent_link(item)
        super().remove(item)

    def pop(self, index: int=-1) -> Node:
        '''
        Extends list pop method with children node validation.

        :param int index: position of the index that is popped out. Default
                          value is -1.

        :returns: the last value of the given index value from the list.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        positiveindex = index if index >= 0 else len(self) - index
        for position in range(positiveindex + 1, len(self)):
            self._validate_item(position-1, self[position])
        self._del_parent_link(self[index])
        return super().pop(index)

    def reverse(self):
        ''' Extends list reverse method with children node validation.'''
        for index, item in enumerate(self):
            self._validate_item(len(self) - index - 1, item)

        super().reverse()

class Node:

    START_DEPTH = 0

    START_POSITION = 0

    def __init__(self, children: List[Node]=None):
        self._children = ChildrenList(self, self._validate_child)
        if children:
            self._children.extend(children)
        self._parent=None

    @staticmethod
    def _validate_child(position: int, child: Node) -> bool:
        '''
        Determines whether a given child and index are valid for this node.
        The generic implementation always returns False, and should be specialised
        in other Nodes.

        :param int position: the position to be validated.
        :param child: a child to be validated.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :return: whether the given child and position are valid for this node.
        :rtype: bool
        '''
        return False

    @property
    def children(self) -> list[Node]:
        '''
        :returns: the children of this Node.
        :rtype: List[:py:class::py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`]
        '''
        return self._children

    @property
    def parent(self) -> Node:
        '''
        :returns: the parent node.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node` or NoneType
        '''
        return self._parent

    def node_str(self) -> str:
        '''
        :returns: a text description of this node. Will typically be \
                overriden by sub-class.
        :rtype: str
        '''
        return type(self).__name__ + "[]"

    def __str__(self) -> str:
        '''
        :returns: the string representation of this node.
        :rtype: str
        '''
        return self.node_str()

    @property
    def depth(self) -> int:
        '''
        Returns this Node's depth in the tree.
        
        :returns: depth of the Node in the tree.
        :rtype: int
        '''

        my_depth = self.START_DEPTH
        node = self
        while node is not None:
            node = node.parent
            my_depth += 1
        return my_depth

    def view(self, depth:int=0, indent: str= "    ", _index: Union[int, None]=None) -> str:
        '''
        Output a human readable description of the current node and all of its descendents as a string.

        :param int depth: depth of the tree hierarch for output text. Default is 0.
        :param str indent: the indent to apply as the depth increases. Defaults to 4 spaces.
        :param int _index: the position of this node wrt its siblins or None. Defaults to None.

        :returns: a string representation of this node and its children.
        :rtype: str

        :raises TypeError: if one of the arguments is the wrong type.
        :raises ValueError: if the depth argument is negative.
        '''

        if not isinstance(depth, int):
            raise TypeError(f"depth argument should be an int but found {type(depth)}.")
        if depth < 0:
            raise ValueError(f"depth argument should be a positive integer but found {depth}.")
        if not isinstance(indent, str):
            raise TypeError(f"indent argument should be a string but found {type(indent)}.")

        full_indent: str = depth*indent
        description = self.node_str()
        if _index is None:
            result = f"{full_indent}{description}\n"
        else:
            result = f"{full_indent}{_index}: {description}\n"
        children_result_list = []
        for idx, node in enumerate(self._children):
            children_result_list.append( node.view(depth=depth+1, _index=idx, indent=indent))
        result = result + "".join(children_result_list)
        return result

    def addchild(self, child: Node, index: Union[None, int]=None) -> None:
        '''
        Adds the supplied node as a child of this node( at position index if supplied).
        The suppled node must not have an existing parent.

        :param child: the node to add as a child of this one.
        :type child: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        :param index: optional position at which to insert new child. Default is
                      None, which means to append the new child to the list of children.
        :type index: Optional[int]
        '''
        if index is not None:
            self._children.insert(index, child)
        else:
            self._children.append(child)

    @property
    def position(self) -> int:
        '''
        Find a Node's position relative to its parent Node (starting with
        0 if it does not have a parent node).

        :returns: relative position of a Node to its parent.
        :rtype: int
        '''
        if self.parent is None:
            return self.START_POSITION
        for index, child in enumerate(self.parent.children):
            if child is self:
                return index

    @property
    def abs_position(self) -> int:
        '''
        Find a Node's absolute position in the tree containing it, starting
        with 0 if it is the root. Is always compute dynamically as its position my change.

        :returns: absolute position of a Node in the tree.
        :rtype: int
        '''
        if self.root is self:
            return self.START_POSITION
        found, position = self._find_position(self.root.children,
                                              self.START_POSITION)

        assert found

        return position

    def _find_position(self, children: list[Node], position: Union[None, int]=None) -> Tuple[bool,int]:
        '''
        Recurse through the tree depth first returning position of self
        if found.

        :param children: list of Nodes that are children of a Node.
        :type children: list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        :param int position: stat counting from this position. Defaults to START_POSITION.

        :returns: the relative position of this node in the tree and if it was found in
                  the tree.
        :rtype: bool, int

        '''

        if position is None:
            position=self.START_POSITION
        for child in children:
            position += 1
            if child is self:
                return True, position
            if child.children:
                found, position = self._find_position(child.children, position)
                if found:
                    return True, position
        return False, position

    @property
    def root(self) -> Node:
        '''
        :returns: the root node of the PIR tree.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        if self.parent is None:
            return self
        node = self.parent
        while node.parent is not None:
            node = node.parent
        return node

    def sameParent(self, node_2: Node) -> bool:
        '''
        :param node_2: The node to check if has the same parent.
        :type node_2: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`

        :returns: True if node_2 has the same parent as this node, False otherwise.
        :rtype: bool
        '''
        if self.parent is None or node_2.parent is None:
            return False
        return self.parent is node_2.parent

    def walk(self, t_type: Union[Tuple[type[Any]], type[Any]], stop_type: Union[None, Tuple[type[Any]], type[Any]]=None) -> List[Node]:
        '''
        Recurse through the PIR tree and return all objects that are an instance of
        t_type, which is either a single class or a tuple of classes.

        The recursion is stopped in a section of the tree if an instance of 'stop_type' is found.

        :param t_type: the class(es) for which the instance are collected.
        :type t_type: type | Tuple[type,...]
        :param stop_type: class(es at which recurtion is halted (optional).
        :type stop_type: Options[type | Tuple[type,...]]

        :returns: list with all nodes that are isntances of t_type starting at \
                  and including this node.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''

        local_list = []
        if isinstance(self, t_type):
            local_list.append(self)
        if stop_type and isinstance(self, stop_type):
            return local_list
        for child in self.children:
            local_list += child.walk(t_type, stop_type)
        return local_list

    def ancestor(self, t_type: Union[Tuple[type[Any]], type[Any]], excluding: Union[None, Tuple[type[Any]], type[Any]]=None,
            include_self: bool=False, limit: Union[None, Node]=None) -> Union[Node, None]:
        '''
        Search back up the tree and find an ancestor that is an instance of the supplied type.
        Return it if found, else return None.
        Node types to ignore may be provided via the `excluding` argument. If `include_self` is True
        then the current node is included in the search.
        If `limit` is provided then the search ceases if the supplied node is encountered.

        :param t_type: class(es) to search for.
        :type t_type: type | Tuple[type,...]
        :param excluding: class(es) to ignore or None.
        :type excluding: Optional[type | Tuple[type, ...]]
        :param bool include_self: whether or not to include this node in the search.
        :param limit: an optional node at which to stop the search.
        :type limit: Optional[:py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`]

        :returns: The first ancestor node that is an instance of the requested classes
                  or None if not found.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node` or None.
        '''

        if include_self:
            myparent = self
        else:
            myparent = self.parent

        if excluding is not None:
            if isinstance(excluding, type):
                excludes = (excluding, )
            elif isinstance(excluding, tuple):
                excludes = excluding
            else:
                raise TypeError(f"The 'excluding' argument to ancestor() must be a type or tuple of types."
                                f" Was supplied {type(excluding)}.")

        if limit and not isinstance(limit, Node):
            raise TypeError(f"The 'limit' argument to ancestor() must be an instance of Node, but got "
                            f"{type(limit)}.")

        if limit:
            while myparent is not None and myparent is not limit:
                if isinstance(myparent, t_type):
                    if not (excluding and isinstance(myparent, excludes)):
                        return myparent
                myparent = myparent.parent
        else:
            while myparent is not None:
                if isinstance(myparent, t_type):
                    if not (excluding and isinstance(myparent, excludes)):
                        return myparent
                myparent = myparent.parent

        return None 

    def pop_all_children(self) -> List[Node]:
        '''
        Remove all children of thise node and retrurn them as a list.

        :returns: all of the children of this node as orphans.
        :rtype list of :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        children = []
        while self.children:
            children.insert(0, self.children.pop())
        return children

    def detach(self) -> Node:
        '''
        Detach this node from the tree it belongs to and return it.
        
        :returns: this node as an orphan.
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.Node`
        '''
        if self.parent:
            index = self.position
            self.parent.children.pop(index)
        return self

    #TODO copy??

    def validate_constraints(self) -> None:
        '''
        Validates this node in the context of the PIR tree.
        This checks constraints that can only be checked once the
        tree is complete.

        By default this routine does nothing, and must be overridden if
        required.
        '''
        pass

class DataNode(Node, metaclass=ABCMeta):
    '''
    Abstract root node for any Node involving data accesses, e.g. Reference
    '''
