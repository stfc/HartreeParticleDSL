class SingletonInstanceError(Exception):
    pass

class RepeatedNameError(Exception):
    pass

class IllegalLoopError(Exception):
    pass

class UnsupportedCodeError(Exception):
    pass

class IllegalArgumentCountError(Exception):
    pass

class UnsupportedTypeError(Exception):
    pass

class InvalidNameError(Exception):
    pass

class InternalError(Exception):
    pass

class NoBackendError(Exception):
    pass

class IRGenerationError(Exception):
    pass
