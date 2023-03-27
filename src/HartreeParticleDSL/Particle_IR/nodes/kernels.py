'''
This module contains all of the Kernel class types.
'''
# pylint: disable=undefined-variable

from __future__ import annotations

from typing import List

from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.reference import Reference
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

class PairwiseKernel(Kern):
    '''
    Class to represent a Pairwise Kernel.

    :param str name: the name of this Kernel.
    '''

    def __init__(self, name: str) -> None:
        super().__init__()
        self._name = name

    @property
    def name(self) -> str:
        '''
        :returns: The name of this pairwise kernel.
        :rtype: str
        '''
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        '''
        Sets the name of this PairwiseKernel

        :param str name: The name of this pairwise kernel.

        :raises TypeError: If the provided name is not a str.
        '''
        if not isinstance(name, str):
            raise TypeError("Expected PairwiseKernel name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @property
    def arguments(self) -> List[Reference]:
        '''
        :returns: The argument list of this pairwise kernel.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[Reference]) -> None:
        '''
        Sets the arguments of this pairwise Kernel. There is a minimum of three
        arguments required, to represent the two particles and the config.

        :param arguments: list of arguments for this pairwise kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises TypeError: If any provided argument is not a Reference.
        :raises IRGenerationError: If there are not at least three arguments.
        '''
        self._arguments = []
        if len(arguments) < 3:
            raise IRGenerationError("Pairwise kernel requires at least three arguments"
                                    f", but only got {len(arguments)}.")
        for arg in arguments:
            if not isinstance(arg, Reference):
                raise TypeError("Each argument must be a Reference, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
    def create(name: str, arguments: List[Reference], kernel_body: List[Statement]) \
            -> PairwiseKernel:
        '''
        Creates a PairwiseKernel containing the input arguments and kernel body.

        :param str name: The name of this PairwiseKernel.
        :param arguments: The list of arguments to use for this Pairwise Kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        :param kernel_body: The list of Statements that make up this Pairwise Kernel
        :type kernel_body: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

        :returns: The new Pairwise Kernel created
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kernels.PairwiseKernel
        '''
        kernel = PairwiseKernel(name)
        kernel.arguments = arguments

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel

class PerPartKernel(Kern):
    '''
    Class to represent a Per Part Kernel.
    '''
    def __init__(self, name: str) -> None:
        super().__init__()
        self._name = name

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
            raise TypeError("Expected PerPartKernel name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @property
    def arguments(self) -> List[Reference]:
        '''
        :returns: The argument list of this perpart kernel.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[Reference]) -> None:
        '''
        Sets the arguments of this perpart Kernel. There is a minimum of two
        arguments required, to represent the particle and the config.

        :param arguments: list of arguments for this perpart kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises TypeError: If any provided argument is not a Reference.
        :raises IRGenerationError: If there are not at least two arguments.
        '''
        self._arguments = []
        if len(arguments) < 2:
            raise IRGenerationError("Perpart kernel requires at least two arguments"
                                    f", but only got {len(arguments)}.")
        for arg in arguments:
            if not isinstance(arg, Reference):
                raise TypeError("Each argument must be a Reference, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
    def create(name: str, arguments: List[Reference], kernel_body: List[Statement]) \
            -> PerPartKernel:
        '''
        Creates a PerPartKernel containing the input arguments and kernel body.

        :param str name: The name of this PerPartKernel.
        :param arguments: The list of arguments to use for this PerPart Kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        :param kernel_body: The list of Statements that make up this PerPart Kernel
        :type kernel_body: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

        :returns: The new PerPart Kernel created
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kernels.PerPartKernel
        '''
        kernel = PerPartKernel(name)
        kernel.arguments = arguments

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel

class MainKernel(Kern):
    '''
    Class to represent the main function of the Particle method.

    :param str name: the name of this Kernel.
    '''

    def __init__(self, name: str) -> None:
        super().__init__()
        self._name = name

    @property
    def name(self) -> str:
        '''
        :returns: The name of this main kernel.
        :rtype: str
        '''
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        '''
        Sets the name of this MainKernel

        :param str name: The name of this main kernel.

        :raises TypeError: If the provided name is not a str.
        '''
        if not isinstance(name, str):
            raise TypeError("Expected MainKernel name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @property
    def arguments(self) -> List[Reference]:
        '''
        :returns: The argument list of this main kernel.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[Reference]) -> None:
        '''
        Sets the arguments of this main Kernel. Must be empty

        :param arguments: list of arguments for this perpart kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises IRGenerationError: If there are any provided arguments.
        '''
        self._arguments = []
        if len(arguments) > 0:
            raise IRGenerationError("Main kernel requires zero arguments"
                                    f", but got {len(arguments)}.")

    @staticmethod
    def create(name: str, kernel_body: List[Statement]) -> MainKernel:
        '''
        Creates a MainKernel containing the input arguments and kernel body.

        :param str name: The name of this MainKernel.
        :param kernel_body: The list of Statements that make up this Main Kernel
        :type kernel_body: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

        :returns: The new Main Kernel created
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kernels.MainKernel
        '''
        kernel = MainKernel(name)

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel


class SourceBoundaryKernel(Kern):
    def __init__(self, name: str, source_count: int) -> None:
        super().__init__()
        self._name = name
        self._source_count = source_count

    @property
    def name(self) -> str:
        '''
        '''
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        '''
        '''
        if not isinstance(name, str):
            raise TypeError("Expected SourceBoundaryKernel name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @property
    def source_count(self) -> int:
        '''
        '''
        return self._source_count

    @source_count.setter
    def source_count(self, source_count: int) -> None:
        if not isinstance(source_count, int):
            raise TypeError("Expected SourceBoundaryKernel source_count to be an "
                            f"int but got {type(source_count)}.")
        self._source_count = source_count

    @property
    def arguments(self) -> List[Reference]:
        '''
        :returns: The argument list of this perpart kernel.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[Reference]) -> None:
        '''

        :param arguments: list of arguments for this perpart kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises TypeError: If any provided argument is not a Reference.
        :raises IRGenerationError: If there are not at least two arguments.
        '''
        self._arguments = []
        if len(arguments) < 2:
            raise IRGenerationError("Source boundary kernel requires at least two arguments"
                                    f", but only got {len(arguments)}.")
        for arg in arguments:
            if not isinstance(arg, Reference):
                raise TypeError("Each argument must be a Reference, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
    def create(name: str, source_count: int, arguments: List[Reference],
               kernel_body: List[Statement]) -> SourceBoundaryKernel:
        '''
        '''
        kernel = SourceBoundaryKernel(name, source_count)
        kernel.arguments = arguments

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel

class SinkBoundaryKernel(Kern):
    def __init__(self, name: str) -> None:
        super().__init__()
        self._name = name

    @property
    def name(self) -> str:
        '''
        '''
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        '''
        '''
        if not isinstance(name, str):
            raise TypeError("Expected SinkBoundaryKernel name to be a str but got "
                            f"{type(name)}.")
        self._name = name

    @property
    def arguments(self) -> List[Reference]:
        '''
        :returns: The argument list of this perpart kernel.
        :rtype: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        '''
        return self._arguments

    @arguments.setter
    def arguments(self, arguments: List[Reference]) -> None:
        '''

        :param arguments: list of arguments for this perpart kernel.
        :type arguments: List of \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises TypeError: If any provided argument is not a Reference.
        :raises IRGenerationError: If there are not at least two arguments.
        '''
        self._arguments = []
        if len(arguments) < 2:
            raise IRGenerationError("Sink boundary kernel requires at least two arguments"
                                    f", but only got {len(arguments)}.")
        for arg in arguments:
            if not isinstance(arg, Reference):
                raise TypeError("Each argument must be a Reference, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
    def create(name: str, arguments: List[Reference], kernel_body: List[Statement]) \
            -> SinkBoundaryKernel:
        '''
        '''
        kernel = SinkBoundaryKernel(name)
        kernel.arguments = arguments

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel
