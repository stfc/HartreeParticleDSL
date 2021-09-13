from abc import ABCMeta, abstractmethod

class C_AOS_IO_Mixin(metaclass=ABCMeta):

    @abstractmethod
    def gen_code_c(self, part_type):
        '''
        Generates and returns the C code required for this IO module.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code C_AOS code required for this IO module.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "gen_code_c")


    @abstractmethod
    def call_input_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for input.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_input_c")

    @abstractmethod
    def call_output_c(self, part_count, filename):
        '''
        Returns the C call required to use this IO module for output.

        :raises NotImplementedError: Abstract method that must be 
                                     overriden by children

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_output_c")
    
    @abstractmethod
    def get_includes_c(self):
        '''
        Returns the C includes required to use this IO module for C_AOS.

        :raises NotImplementedError: Abstract method that must be 
                                     overriden by children

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "get_includes_c")
