import ast
from abc import ABCMeta

class baseVisitor(ast.NodeVisitor, metaclass=ABCMeta):

    def __init__(self, indent=4):
        pass

class main_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class pairwise_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass

class per_part_baseVisitor(baseVisitor, metaclass=ABCMeta):
    pass
