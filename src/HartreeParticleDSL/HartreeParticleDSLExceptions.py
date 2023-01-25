class SingletonInstanceError(Exception):
    '''
    Error used to note when a user attempts to make the protected
    _HartreeParticleDSL object which is forbidden.
    '''
    pass

class RepeatedNameError(Exception):
    '''
    Inbuilt Error thrown when attemping to declare an object
    with a name that already exists in old-style backends.
    '''
    pass

class IllegalLoopError(Exception):
    '''
    Inbuilt Error thrown when attemping to create an unsupported
    type of Loop (usually this means a non-range loop) in
    old-style backends.
    '''
    pass

class UnsupportedCodeError(Exception):
    '''
    Inbuilt Error thrown when attempting to parse code that isn't
    supported yet in an old-style backend.
    '''
    pass

class IllegalArgumentCountError(Exception):
    '''
    Inbuilt Error thrown when a Kernel has the wrong number of
    arguments in old-style backends.
    '''
    pass

class UnsupportedTypeError(Exception):
    '''
    Inbuilt Error thrown when attempting to create a variable of
    an unknown type in an old-style backend.
    '''
    pass

class InvalidNameError(Exception):
    '''
    Inbuilt Error thrown when attempting to create a variable with
    a name not supported by an old-style backend.
    '''
    pass

class InternalError(Exception):
    '''
    An Internal Error thrown in HartreeParticleDSL backends.
    '''
    pass

class NoBackendError(Exception):
    '''
    Currently unused error.
    '''
    pass

class IRGenerationError(Exception):
    '''
    An general error thrown by Particle_IR when illegal IR trees
    are discovered. Error messages explain the actual error found.
    '''
    pass
