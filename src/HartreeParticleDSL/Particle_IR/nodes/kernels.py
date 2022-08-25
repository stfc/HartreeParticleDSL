from HartreeParticleDSL.Particle_IR.nodes.body import Body
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.nodes.reference import Reference

class PairwiseKernel(Kern):
    '''
    Class to represent a Pairwise Kernel.

    :param str name: the name of this Kernel.
    '''

    def __init__(self, name: str) -> None:
        super().__init__()
        self.name = name

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
    def name(self) -> str:
        '''
        :returns: The name of this pairwise kernel.
        :rtype: str
        '''
        return self._name

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
        Sets the arguments of this pairwise Kernel. There is a minimum of two
        arguments required, to represent the particle and the config.

        :param arguments: list of arguments for this pairwise kernel.
        :type arguments: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`

        :raises TypeError: If any provided argument is not a Reference.
        :raises IRGenerationError: If there are not at least two arguments.
        '''
        self._arguments = []
        if len(arguments) < 2:
            raise IRGenerationError("Pariwise kernel requires at least two arguments"
                                    f", but only got {len(arguments)}.")
        for arg in arguments:
            if not isinstance(arg, Reference):
                raise TypeError("Each argument must be a Reference, but found "
                                f"{type(arg)}.")
            self._arguments.append(arg)

    @staticmethod
    def create(self, name: str, arguments: List[Reference], kernel_body: List[Statement]) -> PairwiseKernel:
        '''
        Creates a PairwiseKernel containing the input arguments and kernel body.

        :param str name: The name of this PairwiseKernel.
        :param arguments: The list of arguments to use for this Pairwise Kernel.
        :type arguments: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.reference.Reference`
        :param kernel_body: The list of Statements that make up this Pairwise Kernel
        :type kernel_body: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

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
    # TODO Init with argument list

class MainKernel(Kern):
    '''
    Class to represent the main function of the Particle method.

    :param str name: the name of this Kernel.
    '''

    def __init__(self, name: str) -> None:
        super().__init__()
        self.name = name

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
    def name(self) -> str:
        '''
        :returns: The name of this main kernel.
        :rtype: str
        '''
        return self._name

    @staticmethod
    def create(self, name: str, kernel_body: List[Statement]) -> MainKernel:
        '''
        Creates a MainKernel containing the input arguments and kernel body.

        :param str name: The name of this MainKernel.
        :param kernel_body: The list of Statements that make up this Main Kernel
        :type kernel_body: List of :py:class:`HartreeParticleDSL.Particle_IR.nodes.statement.Statement`

        :returns: The new Main Kernel created
        :rtype: :py:class:`HartreeParticleDSL.Particle_IR.nodes.kernels.MainKernel
        '''
        kernel = MainKernel(name)

        for node in kernel_body:
            kernel.body.addchild(node)

        return kernel
