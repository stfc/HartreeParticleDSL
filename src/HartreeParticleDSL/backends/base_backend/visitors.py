import ast
from abc import ABCMeta
import re

class baseVisitor(ast.NodeVisitor, metaclass=ABCMeta):

    def __init__(self, parent, indent=4):
        self._parent = parent
        super().__init__()

    def check_position(self, node):
        rval = ""
        if type(node.value) is ast.Attribute:
            rval = rval + self.visit(node.value)
            rval = rval + f".{node.attr}"
        a = re.match("core_part.position.[xyz]$", rval[rval.find("core_part.position"):])
        if a is not None:
            end = self._parent.get_particle_access(a.group(0)[-1])
            return rval[0:rval.find("core_part.position")] + end
        return None

class main_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class pairwise_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class per_part_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass
