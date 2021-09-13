from abc import ABCMeta, abstractmethod

class IO_Module(metaclass=ABCMeta):
    ''' Base implementation of an IO Module. '''

    @abstractmethod
    def __init__(self):
        pass
