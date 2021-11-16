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
            rval = rval + str(self.visit(node.value))
            rval = rval + f".{node.attr}"
        # Search for the core_part.position.[...] in the string.
        # Upon locating it, we pull out the contents of the [ ]
        a = re.match("core_part.position.[a-zA-Z]*$", rval[rval.find("core_part.position"):])
        if a is not None:
            print(a.group(0), a.group(0)[(a.group(0).find("position.") +9 ):])
            end = self._parent.get_particle_position(a.group(0)[(a.group(0).find("position.") +9 ):])
            return rval[0:rval.find("core_part.position")] + end
        return None

class main_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class pairwise_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class per_part_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass
