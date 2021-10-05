from abc import ABCMeta, abstractmethod

class FDPS_IO_Mixin(metaclass=ABCMeta):

    def gen_code_fdps(self, part_type):
        '''
        Generates and returns the FDPS code required for this IO module.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code FDPS code required for this IO module.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "gen_code_fdps")

    def call_input_fdps(self, part_count, filename):
        '''
        Returns the call required to use this IO module for input.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_input_fdps")

    def call_output_fdps(self, part_count, filename):
        '''
        Returns the call required to use this IO module for output.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_output_fdps")

    def get_includes_fdps(self):
        '''
        Returns the includes required to use this IO module for FDPS.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "get_includes_fdps")
