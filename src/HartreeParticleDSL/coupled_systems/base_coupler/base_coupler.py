from abc import ABCMeta

class base_coupler(metaclass=ABCMeta):
    
    def has_preferred_decomposition(self):
        return False

    def get_preferred_decomposition(self, field_str, current_indent, indent):
        raise NotImplementedError()

    def get_extra_symbols(self, function_list):
        return []

    def copy_files(self):
        raise NotImplementedError()

    def compilation_files(self):
        return []

    def get_required_packages(self):
        return []
