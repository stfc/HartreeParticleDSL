from __future__ import annotations
from collections import OrderedDict
import inspect
from typing import Union, Dict

import HartreeParticleDSL.HartreeParticleDSL as HartreeParticleDSL
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from HartreeParticleDSL.Particle_IR.datatypes.datatype import DataType
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol


class SymbolTable():
    '''
    Symbol table object for Particle IR.
    Symbols must be unique by name.
    Scopes are only handled at a Kernel level or a global level.

    :param kern: The Kernel or Main or HartreeParticleDSL that this symbol \
                 table belongs to.
    :type kern: :py:class:`HartreeParticleDSL.HartreeParticleDSL.HartreeParticleDSL` or \
                :py:class:`HartreeParticleDSL.Particle_IR.nodes.Kern`

    :raises TypeError: if the kern argument is not an accepted type.
    '''

    def __init__(self, kern: Union[HartreeParticleDSL, Kern]) -> None:
        self._symbols = OrderedDict()

        # Ordered list of the arguments to keep track of inputs.
        self._argument_list = []

        # Store if we're a global symbol table.
        self._is_global = False

        if not isinstance(kern, (HartreeParticleDSL._HartreeParticleDSL, Kern)):
            raise TypeError("Argument 'kern' should be of type "
                            "HartreeParticleDSL or Kern, but instead got "
                            f"{type(kern)}.")
        
        if isinstance(kern, HartreeParticleDSL._HartreeParticleDSL):
            self._is_global = True

        self._kern = kern

    def is_empty(self) -> bool:
        '''
        :returns: True if the symbol table is empty, and False otherwise.
        :rtype: bool
        '''
        return len(self._symbols) == 0

    def get_symbols(self) -> Dict[str, Symbol]:
        '''
        Return the symbols from this symbol table. Note this is a copy of the data
        structure containing this information, not the original.

        :returns: ordered dictionary of symbols indexed by symbol name.
        :rtype: OrderedDict[str] = :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
        '''
        all_symbols = OrderedDict()
        for symbol_name, symbol in self._symbols.items():
            all_symbols[symbol_name] = symbol

        return all_symbols

    def add(self, new_symbol: Symbol):
        ''' Add a new symbol to the symbol table if the symbol name is not already
        in use.
    
        :param new_symbol: the symbol to add to the symbol table.
        :type new_symbol: :py:class:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
    
        :raises TypeError: if the new_symbol argument is not a symbol.
        :raises IRGenerationError: if the symbol table already contains a symbol of this name.
        '''
        if not isinstance(new_symbol, Symbol):
            raise IRGenerationError(f"Symbol {new_symbol} is not a symbol, but "
                                f"{type(new_symbol)}.")
    
        if new_symbol.name in self._symbols:
            raise IRGenerationError(f"Tried to add a new symbol {new_symbol.name} "
                                    "but it was already present in the symbol table.")
    
        self._symbols[new_symbol.name] = new_symbol
    
    def new_symbol(self, name: str, datatype: DataType, symbol_type: type[Symbol]) ->  Symbol:
        '''
        Create a new symbol of symbol type `symbol_type`, with name `name` and datatype `datatype`.
        The visibility will match the visibility of this symbol table.
    
        :param name: Name of the symbol. This name must be unique in this symbol table.
        :type name: str
        :param datatype: The data type of the symbol.
        :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.DataType`
        :param symbol_type: The Particle IR type of the symbol to create.
        :type symbol_type: :py:type:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`
    
        :raises TypeError: if the type_symbol argument is not the type of a Symbol \
                           object class or one of its subclasses.
        :raises TypeError: if the name argument is not a str.
        '''
        if not (isinstance(symbol_type, type) and
                Symbol in inspect.getmro(symbol_type)):
            raise TypeError(f"The symbol_type parameter should be a type class of "
                            f"a sublcass of Symbol, but found {type(symbol_type)} instead.")
        if not isinstance(name, str):
            raise TypeError(f"The name parameter should be a str, but found {type(name)}.")
    
        visibility = Symbol.Visibility.LOCAL
        if self._is_global:
            visibility = Symbol.Visibility.GLOBAL
    
        symbol = symbol_type(name=name, datatype=datatype, visibility=visibility)
        self.add(symbol)

        return symbol

    def lookup(self, name: str) -> Union[Symbol, NoneType]:
        '''
        Search for a Symbol named name in the symbol table. Returns the symbol
        if it exists, otherwise returns None.

        :param name: Name of the symbol to lookup.
        :type name: str

        :returns: The symbol found in the symbol table.
        :rtype: Symbol or None.

        :raises TypeError: If the name argument is not a string.
        '''
        if not isinstance(name, str):
            raise TypeError("Expected the name argument to be a string, but got "
                            f"{type(name)}.")

        # TODO Should we look in the global symbol table if we don't find locally?
        return self._symbols.get(name)            

    def find_or_create(self, name: str, datatype: DataType, symbol_type: type[Symbol]) ->  Symbol:
        '''
        Search for a symbol in the symbol table and return it, else create it if it doesn't
        yet exist.

        :param name: Name of the symbol.
        :type name: str
        :param datatype: The data type of the symbol.
        :type datatype: :py:class:`HartreeParticleDSL.Particle_IR.datatypes.datatype.DataType`
        :param symbol_type: The Particle IR type of the symbol to create.
        :type symbol_type: :py:type:`HartreeParticleDSL.Particle_IR.symbols.symbol.Symbol`

        :raises IRGenerationError: if the symbol found doesn't match the datatype requested
        '''

        sym = self.lookup(name)
        if sym is not None:
            if sym.datatype != datatype:
                raise IRGenerationError(f"Found a symbol with specified name {name}, but "
                                        f"it had datatype {sym.datatype} which doesn't "
                                        f"match requested datatype {datatype}.")
            return sym
        return self.new_symbol(name=name, datatype=datatype, symbol_type=symbol_type)


    # Future potential requirements
    # remove
    # merge
    # argument_list property
