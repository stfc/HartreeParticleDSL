import pytest

from HartreeParticleDSL.Particle_IR.nodes.node import Node, ChildrenList
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

'''
Module for testing the children list and Node implementations.
Unit testing quality is quite bad here as Node doesn't support a lot
of the functionality itself, so a fake child class was implemented
'''


def test_children_list():
    # Create a dummy validation function
    def validation_function(item, node):
        return True
    # Create a failing validation function
    def validation2(item, node):
        return False
    node1 = Node()
    node2 = Node()
    node3 = Node()
    a = ChildrenList(node1, validation_function)
    b = ChildrenList(node2, validation2)

    a._validate_item(2, node2)

    with pytest.raises(IRGenerationError) as excinfo:
        b._validate_item(1, node3)
    assert ("Item '<class 'HartreeParticleDSL.Particle_IR.nodes.node.Node'>' can't "
            "be child 1 of '<class 'HartreeParticleDSL.Particle_IR.nodes.node.Node'>'."
            in str(excinfo.value))

    a.insert(0,node3)

    a._check_is_orphan(node1)

    with pytest.raises(IRGenerationError) as excinfo:
        a._check_is_orphan(node3)
    assert ("Item 'Node[]' can't be added as child of 'Node[]' "
            "because it is not an orphan. It already has a 'Node[]'"
            " as a parent." in str(excinfo.value))

    a.pop(0)
    assert node3.parent is None

    a._set_parent_link(node3)
    assert node3.parent == node1

    a._del_parent_link(node3)
    assert node3.parent is None

    a.append(node3)
    assert len(a) == 1
    assert node3.parent is node1
    
    a.append(node2)
    assert len(a) == 2
    assert a[0] is node3
    assert a[1] is node2
    assert node2.parent is node1

    a.pop(0)
    a.pop(0)
    assert len(a) == 0
    assert node2.parent is None
    assert node3.parent is None

    a.append(node2)
    a.__setitem__(0, node3)
    assert len(a) == 1
    assert a[0] is node3
    assert node2.parent is None
    a.pop(0)

    a.append(node2)
    a.insert(0, node3)
    assert len(a) == 2
    assert a[0] is node3
    assert a[1] is node2
    a.pop(0)
    a.pop(0)

    x = [node2, node3]
    a.extend(x)
    assert len(a) == 2
    assert a[0] is node2
    assert a[1] is node3

    del(a[1])
    assert len(a) == 1
    assert node3.parent is None
    assert a[0] is node2
    a.append(node3)
    del(a[0])
    assert a[0] is node3
    del(a[0])
    a.append(node2)
    
    a.remove(node2)
    assert len(a) == 0
    assert node2.parent is None
    a.append(node2)
    a.append(node3)
    a.remove(node2)
    assert a[0] is node3
    a.remove(node3)

    a.extend(x)
    a.reverse()
    assert a[0] is node3
    assert a[1] is node2


def test_node_init():
    node = Node()
    assert node.START_DEPTH == 0
    assert node.START_POSITION == 0

    assert len(node._children) == 0
    assert node._parent is None

def test_node_validate_child():
    assert Node._validate_child(0, "a") is False

def test_node_init_children1():
    class fakeNode(Node):
        @staticmethod
        def _validate_child(position, child):
            return True

    n1 = fakeNode()
    n2 = fakeNode()

    mynode = fakeNode(children=[n1, n2])
    assert len(mynode._children) == 2
    assert mynode._children[0] is n1
    assert mynode._children[1] is n2

    assert mynode.children is mynode._children

def test_node_nodestr():
    node = Node()
    assert node.node_str() == "Node[]"
    assert str(node) == "Node[]"

def test_node_init_children2():
    class fakeNode(Node):
        @staticmethod
        def _validate_child(position, child):
            return True

    class fakeNode2(Node):
        @staticmethod
        def _validate_child(position, child):
            return True

    n1 = fakeNode()
    n2 = fakeNode()

    mynode = fakeNode(children=[n1, n2])
    assert mynode.depth is 1
    assert n1.depth is 2

    correct = '''fakeNode[]
    fakeNode[]
    fakeNode[]
'''
    assert mynode.view() == correct

    with pytest.raises(TypeError) as excinfo:
        mynode.view(depth="x")
    assert ("depth argument should be an int but found str." in
            str(excinfo.value))
    with pytest.raises(ValueError) as excinfo:
        mynode.view(depth=-1)
    assert ("depth argument should be a positive integer but found -1." in
            str(excinfo.value))
    with pytest.raises(TypeError) as excinfo:
        mynode.view(indent=3)
    assert ("indent argument should be a str but found int." in
            str(excinfo.value))

    mynode.children.pop(1)
    assert n2.parent is None
    mynode.addchild(n2)
    assert n2.parent is mynode
    assert mynode.children[0] is n1
    mynode.children.pop(1)
    mynode.addchild(n2, 0)
    assert mynode.children[0] is n2

    assert mynode.position is 0
    assert n2.position is 0
    assert n1.position is 1

    n3 = fakeNode()
    n4 = fakeNode()
    n5 = fakeNode()
    n2.addchild(n3)
    n1.addchild(n4)
    n1.addchild(n5)

    assert mynode.abs_position is 0
    assert n2.abs_position is 1
    assert n3.abs_position is 2
    assert n1.abs_position is 3
    assert n4.abs_position is 4
    assert n5.abs_position is 5

    assert n5.root is mynode

    assert n5.same_parent(n4) is True
    assert n5.same_parent(n3) is False

    assert len(mynode.walk(Node)) is 6
    x = fakeNode2()
    n2.addchild(x)
    z = fakeNode()
    x.addchild(z)
    x.addchild(fakeNode())
    assert len(mynode.walk(Node, stop_type=fakeNode2)) is 7

    assert n5.ancestor(Node) is n1
    assert n5.ancestor(Node,include_self=True) is n5
    assert n5.ancestor(Node,excluding=fakeNode) is None
    assert n5.ancestor(Node,excluding=(fakeNode, Node)) is None
    assert n5.ancestor(int) is None
    assert z.ancestor(Node,limit=n2) is x
    assert z.ancestor(fakeNode, limit=mynode) is n2

    with pytest.raises(TypeError) as excinfo:
        n5.ancestor(Node,excluding=123)
    assert ("The 'excluding' argument to ancestor() must be a type or a tuple of types but got: 'int'" in str(excinfo.value))

    with pytest.raises(TypeError) as excinfo:
        n5.ancestor(Node,limit=123)
    assert("The 'limit' argument to ancestor() must be an instance of Node but got 'int'" in str(excinfo.value))

    li = n1.pop_all_children()
    assert n4 is li[0]
    assert n5 is li[1]

    assert n1.detach() is n1
    assert n1.parent is None

    n1.validate_constraints()
    
    li = [n5]
    assert n5._find_position(li, None) == (True, 1)
    assert n5.same_parent(None) == False
