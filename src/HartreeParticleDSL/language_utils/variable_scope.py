from HartreeParticleDSL.HartreeParticleDSLExceptions import RepeatedNameError
import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL

class variable_access():
    '''
    A class to handle variable accesses. Used whenever variable accesses
    occur in the AST to handle later creation of strings.

    :param variable: The _variable object accessed.
    :type variable: _variable
    :param child: The child variable accesses of this variable access. \
                  For example, if some language code is expected to be \
                  `obj.element` then the variable access to `obj` will have \
                  the variable acces to `element` as its child. Default is None.
    :type child: variable_access or None
    :param array_index: The variable access(es) or constant value(s) used to access this
                        variable access if it is an array access.
    :type array_index: List of variable_access or str
    '''
    def __init__(self, variable, child=None, array_index=None):
        assert variable is not None
        assert (child is None) or (type(child) is variable_access)
        self._variable = variable
        if array_index is None:
            self._array_indices = []
        else:
            self._array_indices = array_index
        self._child = child
        self._is_child = False

    def __str__(self):
        '''
        :returns: The string representation of this variable access. This uses \
                  the currently set backend to create the string.
        '''
        return HartreeParticleDSL.get_backend().access_to_string(self)


    def add_array_index(self, array_index):
        '''
        Add an array index to this variable access
        :param array_index: The array index to add to this variable_access
        :type array_index: variable_access or str

        :raises SyntaxError: If the array_index is not a variable_access or str
        '''
        if not isinstance(array_index, (variable_access, str)):
            raise SyntaxError("The supplied array_index was neither a "
                              f"variable_access or str, received {type(array_index)}")
        assert array_index is not self
        self._array_indices.append(array_index)


    @property
    def array_indices(self):
        '''
        Returns the array indices used for this variable access

        :returns: The array indices used for this variable access.
        :rtype: List of variable_access or str
        '''
        return self._array_indices

    @property
    def variable(self):
        '''
        Returns the variable object associated with this variable access.

        :returns: The variable object used for this variable access.
        :rtype: _variable
        '''
        return self._variable

    @property
    def child(self):
        '''
        Returns the child variable access of this variable access.

        :returns: The child variable access of this variable access.
        :rtype: variable_access or None
        '''
        return self._child

    @child.setter
    def child(self, child):
        '''
        Add a child to this variable access.
        :param child: The child variable access to add to this variable_access
        :type child: variable_access or str

        :raises SyntaxError: If there is already a child.
        '''
        if self._child is not None:
            raise SyntaxError("Variable access already contains a child")
        child.is_child = True
        self._child = child

    @property
    def is_child(self):
        '''
        :returns: Whether this is a child variable.
        :rtype: bool
        '''
        return self._is_child

    @is_child.setter
    def is_child(self, is_child):
        '''
        Sets whether this variable is a child or not.
        :param bool is_child: Whether this variable is a child
        '''
        self._is_child = is_child

class variable():
    '''
    A class to handle each variable. Used by variable_scope to handle
    creating variables

    :param str var_name: The variable's name
    :param str var_type: The variable's type.
    :param bool is_pointer: Determine whether this i
    '''
    def __init__(self, var_name, var_type, is_pointer):
        self._var_name = var_name
        self._var_type = var_type
        self._is_pointer = is_pointer

    @property
    def var_name(self):
        return self._var_name

    @property
    def var_type(self):
        return self._var_type

    @property
    def is_pointer(self):
        return self._is_pointer

class variable_scope():
    '''
    A class to handle variable scoping. Used by backends to keep track of
    variable scope and types to determine whether things are pointers, etc.
    '''

    def __init__(self):
        self._scope = {}

    def add_variable(self, var_name, var_type, is_pointer):
        '''
        Add a variable to the current scope.

        :param str var_name: The variable's name.
        :param str var_type: The variable's type.
        :param bool is_pointer: Whether the variable is a pointer \
                                or not.

        :raises RepeatedNameError: If var_name is already in the scope.
        '''
        if var_name in self._scope:
            raise RepeatedNameError(f"{var_name} is already defined in the "
                                     "current scope")
        var = variable(var_name, var_type, is_pointer)
        self._scope[var_name] = var

    def get_variable(self, var_name):
        '''
        Check if a variable is in the current scope
        
        :param str var_name: The variable's name.

        :returns: The variable object corresponding to the variable name. \
                  If not found then returns None.
        :rtype: variable or None
        '''
        var = self._scope.get(var_name, None)
        return var

    def remove_variable(self, var_name):
        '''
        Removes a variable from the current scope.

        :param str var_name: The variable's name to be removed.

        :raises InvalidNameError: If var_name is not in the scope.
        '''
        if var_name not in self._scope:
            raise InvalidNameError(f"{var_name} is not defined in the current "
                                    "scope.")
        self._scope.pop(var_name, None)
