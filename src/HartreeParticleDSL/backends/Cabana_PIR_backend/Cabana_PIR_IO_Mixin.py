from abc import ABCMeta, abstractmethod

class Cabana_PIR_IO_Mixin(metaclass=ABCMeta):

    def gen_code_cabana_pir(self, part_type):
        '''
        Generates and returns the Cabana_PIR code required for this IO module.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code Cabana_PIR code required for this IO module.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "gen_code_cabana_pir")

    def call_get_box_size_pir(self, part_count, filename, current_indent=4):
        '''
        Returns the code required for the Cabana_PIR backend to use this IO module for reading
        in the box dimensions from the HDF5 file.

        :param int part_count: Unused particle count value.
        :param str filename: Filename to read the data from.
        :param int current_indent: The code current indentation level. Default is 4.

        :returns: The function call used to get the box size from the file.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_get_box_size_pir")

    def call_input_cabana_pir(self, part_count, filename, current_indent=4):
        '''
        Returns the call required to use this IO module for input.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code required to use this IO module for input.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_input_cabana_pir")

    def call_output_cabana_pir(self, part_count, filename, variable, current_indent=4):
        '''
        Returns the call required to use this IO module for output.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The code required to use this IO module for output.
        :rtype: str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "call_output_cabana_pir")

    def get_includes_cabana_pir(self):
        '''
        Returns the includes required to use this IO module for Cabana_PIR.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The includes for this IO module.
        :rtype: List of str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "get_includes_cabana_pir")

    def get_header_includes_cabana_pir(self):
        '''
        Returns the includes required for the header to use this IO module for
        Cabana PIR.

        :raises NotImplementedError: Abstract method that must be
                                     overriden by children

        :returns: The header includes for this IO module.
        :rtype: List of str
        '''
        raise NotImplementedError(f"{self.__class__.__name__} does not "
                                  "implement required function "
                                  "get_header_includes_cabana_pir")
